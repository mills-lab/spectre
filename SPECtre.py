#!/usr/bin/python

# DESCRIPTION: Takes as input a BAM alignment file of ribosome profiling reads,
# gene or transcript abundance estimates from Cufflinks (mRNA-Seq or ribo-Seq)
# and a GTF annotation file (UCSC genes.gtf or Ensembl have been tested. From
# these inputs, the normalized read coverage over each transcript is computed
# and the spectral coherence relative to an enriched signal over a single
# reading frame is calculated (eg. 1:0:0, repeating). Additional implementations
# to calculate the FLOSS read distribution and FLOSS score (Ingolia, 2014), and
# ORFScore (Bazzini, 2014) are included. This version is able to handle one
# sample (BAM) at a time, however an additional script (intersect_experiments.py
# is included to merge multiple experiments into a single analysis.

# If you use SPECtre as part of your analysis, please cite the following:
#

# This script was written and tested using Python v2.7.8 on an OS X v10.9.4
# workstation and on a server running RHEL r6.4.

# DEPENDENCIES:
# bx-python:	https://bitbucket.org/james_taylor/bx-python/wiki/Home
# pyfasta:		https://pypi.python.org/pypi/pyfasta/
# rpy:			http://rpy.sourceforge.net/
# ROCR (in R)	https://rocr.bioinf.mpi-sb.mpg.de/

import sys
import os
import re
import math
import logging
import argparse
import itertools
import collections

from bx.intervals.intersection import IntervalTree

from functools import partial
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

from rpy2.robjects import r
from rpy2.robjects.packages import importr

#############
# UTILITIES #
#############
def hash():
	return collections.defaultdict(hash)

def flatten(d, delimiter=":"):
	def expand(key, value):
		if isinstance(value, dict):
			return [(delimiter.join([key, k]), v) for k, v in flatten(value, delimiter).items()]
		else:
			return [(key, value)]
	return dict([item for k, v in d.items() for item in expand(k, v)])

def add_custom_offsets(custom, default):
	custom_offsets = custom.split(",")
	for offset in custom_offsets:
		read, pos = offset.split(":")
		if "custom" in default:
			if int(read) in default["custom"]:
				pass
			else:
				default["custom"][int(read)] = int(pos)
		else:
			default["custom"] = {int(read): int(pos)}
	return default

def check_chrom_prefix(bam_file):
	header = os.popen("samtools view -H " + bam_file)
	chroms = [line for line in header if "@SQ" in line]
	has_prefix = False
	for chrom in chroms:
		has_prefix = True if "chr" in chroms else False
	return has_prefix

def convert_chromosome(name):
	return name if "chr" in name else "chr" + name

def decode(flag):
	return "-" if int(flag) & 16 else "+"

def transcript_length(coords):
	return sum([(end+1 - start) for start, end in coords])

def expand_transcript_coordinates(coords):
	coordinates = list()
	# Check if the supplied coordinates are of type str():
	if "," in coords:
		coords = coords.split(",")
	if isinstance(coords, str):
		start, end = coords.split("-")
		region = range(int(start), int(end)+1)
		for pos in region:
			coordinates.append(pos)
	elif isinstance(coords, list):
		for coord in coords:
			start, end = coord.split("-") if isinstance(coord, str) else coord
			region = range(int(start), int(end)+1)
			for pos in region:
				coordinates.append(pos)
	else:
		start, end = coords
		region = range(int(start), int(end)+1)
		for pos in region:
			coordinates.append(pos)
	return sorted(coordinates)

def regroup_coordinates(coords):
	first = last = coords[0]
	for n in coords[1:]:
		if n-1 == last:
			last = n
		else:
			yield first, last
			first = last = n
	yield first, last

def calculate_offset_position(read_position, read_strand, read_length, read_cigar, offsets):
	# This function takes as input the 5' position of a read, the read length and strand, and calculates the offset
	# position of the read relative to the A or P-site based on the CIGAR string and pre-defined position offsets based
	# on previously pusblished methods, or user-provided offsets.
	def calculate_offset(offsets, length):
		if length in offsets:
			return offsets[length]
		else:
			# All Bazzini offset positions are the same, Ingolia and other non-standard read offsets should be about
			# half of the read length:
			return offsets.values()[0] if len(set(offsets.values())) == 1 else length / 2

	def extract_coordinates_from_cigar(pos, cigar):
		# Return a list of expanded coordinates based on the read position and its CIGAR string:
		coordinates = list()
		ops = re.findall("[0-9]+[DHIMNPSX=]{1}", cigar)
		for op in ops:
			increment, modifier = int(op[:-1]), op[-1]
			if not modifier == "N":
				increments = range(increment)
				for x in xrange(len(increments)):
					coordinates.append(pos + increments[x])
			else:
				pos = coordinates[-1] + increment
		return coordinates

	return extract_coordinates_from_cigar(read_position, read_cigar)[calculate_offset(offsets, read_length)-1] if read_strand == "+" else extract_coordinates_from_cigar(read_position, read_cigar)[-calculate_offset(offsets, read_length)]

class SAM(object):
	# Parse the default fields from a BAM/SAM alignment record:
	def __init__(self, sam):
		self.sam = sam
	def qname(self):
		return self.sam.strip().split("\t")[0]
	def flag(self):
		return self.sam.strip().split("\t")[1]
	def rname(self):
		return self.sam.strip().split("\t")[2]
	def pos(self):
		return self.sam.strip().split("\t")[3]
	def mapq(self):
		return self.sam.strip().split("\t")[4]
	def cigar(self):
		return self.sam.strip().split("\t")[5]
	def rnext(self):
		return self.sam.strip().split("\t")[6]
	def pnext(self):
		return self.sam.strip().split("\t")[7]
	def tlen(self):
		return self.sam.strip().split("\t")[8]
	def seq(self):
		return self.sam.strip().split("\t")[9]
	def qual(self):
		return self.sam.strip().split("\t")[10]
	def optional(self):
		return self.sam.strip().split("t")[11:]

class Checks(object):
	# Series of QC checks for minimum length and coverage:
	def __init__(self, region, frame_reads, strand, length, buffers):
		self.region = region
		self.frame_reads = frame_reads
		self.strand = strand
		self.length = length
		self.buffers = buffers

	def coverage(self):
		max_left_trim, max_right_trim = max([i[0] for i in self.buffers]), max([j[-1] for j in self.buffers])
		if self.region == "NA":
			return False
		else:
			if self.strand == "+":
				if max_left_trim == 0:
					if max_right_trim == 0:
						return True if math.fsum(self.region) > 0.0 else False
					else:
						return True if math.fsum(self.region[:-max_right_trim]) > 0.0 else False
				else:
					if max_right_trim == 0:
						return True if math.fsum(self.region[max_left_trim-1:]) > 0.0 else False
					else:
						return True if math.fsum(self.region[max_left_trim-1:-max_right_trim]) > 0.0 else False
			else:
				if max_left_trim == 0:
					if max_right_trim == 0:
						return True if math.fsum(self.region) > 0.0 else False
					else:
						return True if math.fsum(self.region[max_right_trim-1:]) > 0.0 else False
				else:
					if max_right_trim == 0:
						return True if math.fsum(self.region[:-max_left_trim]) > 0.0 else False
					else:
						return True if math.fsum(self.region[max_right_trim-1:-max_left_trim]) > 0.0 else False

	def trimming(self):
		# Require that each post-trimmed transcript region be at least the length of the sliding window:
		max_left_trim, max_right_trim = max([i[0] for i in self.buffers]), max([j[-1] for j in self.buffers])
		return True if len(self.region) - max_left_trim - max_right_trim >= int(self.length) else False

	def frames(self):
		if self.frame_reads in (None, "NA", ""):
			return False
		elif sum(self.frame_reads) == 0:
			return False
		else:
			return True

##############
# ANNOTATION #
##############
def extract_fpkms(cufflinks):
	# This function takes as input a Cufflinks isoforms.fpkm_tracking output file and
	# extracts the transcript-level expression in FPKM. Pleast note that, as written,
	# this function should support a generic tab-delimited file such that the name of
	# the transcript is in the first column, and the column holding the expression
	# values is titled as "FPKM".
	logger.info("extract_fpkms(): Parsing transcript FPKMs from file: " + cufflinks + " into memory... [STARTED]")
	fpkms = dict()
	position = 0
	for linenum, line in enumerate(open(cufflinks)):
		if linenum == 0:
			for i, j in enumerate(line.strip().split("\t")):
				if j.upper() == "FPKM":
					position = i
		else:
			transcript, fpkm = line.strip().split("\t")[0], float(line.strip().split("\t")[position])
			if transcript not in fpkms:
				fpkms[transcript] = fpkm
	logger.info("extract_fpkms(): Parsing transcript FPKMs from file: " + cufflinks + " into memory... [COMPLETE]")
	return fpkms

def parse_gtf(gtf_file, fpkms, window_length, buffers, prefix, sanitize):
	# This function takes as input a user-supplied GTF transcript annotation file,
	# and extracts the CDS and UTR coordinates of protein-coding genes, and start
	# and end coordinates for non-coding transcripts then loads then into a dict().
	# In addition, this function will parse the start and end coordinates for each
	# transcript into an IntervalTree(). Default transcript annotation types may be
	# found in main(), and can be changed to suit user requirements.
	def append_stop_coordinates(gtf):
		# Since Ensembl annotation GTFs have separate entries for stop coordinates, these
		# coordinates must be manually appended to the end of each protein-coding transcript.
		def extend_coordinates(coords, strand):
			if strand == "+":
				coords[-1] = (coords[-1][0], coords[-1][-1]+3)
			else:
				coords[0] = (coords[0][0]-3, coords[0][-1])
			return coords
		logger.info("parse_gtf/append_stop_coordinates(): Appending stop coordinates to parsed GTF.. [STARTED]")
		flat = [(transcript, coordinates) for transcript, coordinates in flatten(gtf).items() if "protein_coding" in transcript if "CDS" in transcript]
		for transcript, coordinates in flat:
			gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
			gtf[gene_type][chrom][strand][gene_id][transcript_id][feature] = extend_coordinates(coordinates, strand)
		logger.info("parse_gtf/append_stop_coordinates(): Appending stop coordinates to parsed GTF.. [COMPLETE]")
		return gtf

	def partition_utr_coordinates(gtf):
		# Since Ensembl does not differentiate between the 5'UTR and 3'UTR, they must be
		# annotated based on their position relative to the CDS and strand:
		logger.info("parse_gtf/partition_utr_coordinates(): Partitioning 5' and 3' UTR entries from loaded GTF.. [STARTED]")
		flat = [(transcript, coordinates) for transcript, coordinates in flatten(gtf).items() if "protein_coding" in transcript if "UTR" in transcript]
		for transcript, coordinates in flat:
			gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
			five_prime, three_prime = list(), list()
			try:
				if strand == "+":
					cds_start = gtf[gene_type][chrom][strand][gene_id][transcript_id]["CDS"][0][0]
					[five_prime.append((start, end)) if end <= cds_start else three_prime.append((start, end)) for start, end in coordinates]
				else:
					cds_start = gtf[gene_type][chrom][strand][gene_id][transcript_id]["CDS"][-1][-1]
					[five_prime.append((start, end)) if start >= cds_start else three_prime.append((start, end)) for start, end in coordinates]
				# Remove the previous combined entry for UTRs of a transcript:
				del gtf[gene_type][chrom][strand][gene_id][transcript_id]["UTR"]
				# Replace with separate key entries for 5'UTR and 3'UTR:
				gtf[gene_type][chrom][strand][gene_id][transcript_id]["UTR5"] = five_prime
				gtf[gene_type][chrom][strand][gene_id][transcript_id]["UTR3"] = three_prime
			except KeyError:
				pass
		logger.info("parse_gtf/partition_utr_coordinates(): Partitioning 5' and 3' UTR entries from loaded GTF.. [COMPLETE]")
		return gtf

	def interval_tree(gtf):
		logger.info("parse_gtf/interval_tree(): Parsing transcript intervals from GTF into an IntervalTree()... [STARTED]")
		# Load the start and end coordinates for each transcript into an IntervalTree()
		# object. A 5kb buffer is added to each end per Ingolia, et al. (2014) in order
		# to help filter out proximal transcripts. This buffer can be changed or removed
		# according to user requirements, but is not recommended.
		intervals = dict()
		flat = [(transcript, coordinates) for transcript, coordinates in flatten(gtf).items() if "CDS" in transcript or "exon" in transcript]
		for transcript, coordinates in flat:
			gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
			tree = None
			seq_id = chrom + "|" + strand
			start = coordinates[0][0] - 5000
			end = coordinates[-1][-1] + 5000
			if seq_id in intervals:
				tree = intervals[seq_id]
			else:
				tree = IntervalTree()
				intervals[seq_id] = tree
			tree.add(start, end, (gene_id, transcript_id))
		logger.info("parse_gtf/interval_tree(): Parsing transcript intervals from GTF into an IntervalTree()... [COMPLETE]")
		return intervals

	# Further refinement required...
	def remove_overlapping_transcripts(transcript_gtf, transcript_intervals):
		cleaned_gtf = hash()
		logger.info("parse_gtf/remove_overlapping_transcripts(): Sanitizing transcript database [STARTED].")
		for gene_type in gtf:
			for chrom in gtf[gene_type]:
				for strand in gtf[gene_type][chrom]:
					for gene in gtf[gene_type][chrom][strand]:
						for transcript in gtf[gene_type][chrom][strand][gene]:
							# Use the CDS or transcript sequence as the test region:
							for feature in gtf[gene_type][chrom][strand][gene][transcript]:
								if feature in ("CDS", "exon"):
									# Add a 5000 nt buffer to each transcript start and end coordinate:
									seq_id = chrom + "|" + strand
									start, end = gtf[gene_type][chrom][strand][gene][transcript][feature][0][0] - 5000, gtf[gene_type][chrom][strand][gene][transcript][feature][-1][-1] + 5000
									overlaps = list()
									# Test for overlap with other transcripts:
									if seq_id in transcript_tracker:
										overlaps.extend(transcript_tracker[seq_id].find(start, end))
									# Transcript is allowed to overlap with itself, therefore more than two (2) overlaps
									# should be flagged for removal:
									if len(overlaps) > 1:
										if gene_type not in cleaned_gtf:
											if chrom not in cleaned_gtf[gene_type]:
												if strand not in cleaned_gtf[gene_type][chrom]:
													if gene not in cleaned_gtf[gene_type][chrom][strand]:
														if transcript not in cleaned_gtf[gene_type][chrom][strand][gene]:
															cleaned_gtf[gene_type][chrom][strand][gene][transcript] = gtf[gene_type][chrom][strand][gene][transcript]
		logger.info("parse_gtf/remove_overlapping_transcripts(): Sanitizing transcript database [COMPLETE].")
		return cleaned_gtf

	def filter_transcripts(gtf, transcript_fpkms, transcript_intervals, window_length, buffers, sanitize):
		def check_expression(transcript, transcript_fpkms):
			# Requires that the transcript be expressed:
			if transcript in transcript_fpkms:
				return True if transcript_fpkms[transcript] > 0 else False
			else:
				return False

		def check_length(coordinates, transcript_buffers, window_length):
			# Requires that the CDS or exonic regions of a transcript be of a minimum length to
			# accomodate the sliding SPECtre windows:
			def transcript_length(coordinates):
				length = 0
				for start, end in coordinates:
					length += (end+1) - start
				return length

			max_left_trim = max([left for left, right in transcript_buffers])
			max_right_trim = max([right for left, right in transcript_buffers])
			return True if transcript_length(coordinates) - max_left_trim - max_right_trim >= window_length else False

		logger.info("parse_gtf/filter_transcripts(): Scrubbing transcripts based on minimum length and expression [STARTED].")
		filtered_gtf = hash()
		flat = [(transcript, coordinates) for transcript, coordinates in flatten(gtf).items() if "CDS" in transcript or "exon" in transcript]
		for transcript, coordinates in flat:
			gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
			if check_expression(transcript_id, fpkms) == True:
				# Removed check_length() check because Coherence() should return the full length coherence for short transcripts:
				#if check_length(coordinates, buffers, window_length) == True:
				try:
					filtered_gtf[gene_type][chrom][strand][gene_id][transcript_id] = gtf[gene_type][chrom][strand][gene_id][transcript_id]
				except KeyError:
					pass
		logger.info("parse_gtf/filter_transcripts(): Scrubbing transcripts based on minimum length and expression [COMPLETE].")
		return remove_overlapping_transcripts(filtered_gtf, transcript_intervals) if sanitize == True else filtered_gtf

	# Load transcripts from the provided GTF into a hash():
	logger.info("parse_gtf(): Parsing transcript coordinates from GTF: " + gtf_file + " into memory... [STARTED]")
	transcripts = hash()
	for line in open(gtf_file):
		if not line.startswith("#"):
			seq_name, source, feature, start, end, score, strand, frame, attributes = line.strip().split("\t")
			gene_type = re.findall('gene_[a-z]{0,3}type\s\"[\w\.]+\"', attributes)[0].split('"')[1]
			# Parse annotated protein-coding CDS into the GTF dictionary:
			if (gene_type == "protein_coding" and feature in ("CDS", "UTR")) or (gene_type != "protein_coding" and feature == "exon"):
				gene, transcript = re.findall('gene_id\s\"[\w\_\:\:\.]+\"', attributes)[0].split('"')[1], re.findall('transcript_id\s\"[\w\W]+\"', attributes)[0].split('"')[1]
				if isinstance(transcripts[gene_type][seq_name][strand][gene][transcript][feature], list):
					transcripts[gene_type][seq_name][strand][gene][transcript][feature].append((int(start), int(end)))
				else:
					transcripts[gene_type][seq_name][strand][gene][transcript][feature] = [(int(start), int(end))]
	logger.info("parse_gtf(): Parsing transcript coordinates from GTF: " + gtf_file + " into memory... [COMPLETE]")
	# Instantiate the transcript IntervalTree():
	intervals = interval_tree(transcripts)
	# Filter transcripts GTF based on minimum length, expression, and sanitize of overlapping transcripts, as requested:
	return filter_transcripts(partition_utr_coordinates(append_stop_coordinates(transcripts)), fpkms, intervals, window_length, buffers, sanitize)

#############
# ALIGNMENT #
#############
def load_alignments_to_tree(bam_file):
	# This function loads alignments from the input BAM file into an IntervalTree for faster access. Timer() results
	# are shown below for an abnormally large region, but illustrates the potential savings in time and computation
	# over thousands of transcripts or ORFs:
	#
	# t_tree = Timer(lambda: reads["X|+"].find(1,10000000))
	# print t_tree.timeit(number=1)
	# 2.8133392334e-05
	# t_popen = Timer(lambda: os.popen("samtools view " + bam + " X:1-10000000"))
	# print t_popen.timeit(number=1)
	# 0.0262911319733
	intervals = dict()
	logger.info("load_alignments_to_tree(): Loading alignments from: " + bam_file + " into IntervalTree()... [STARTED].")
	reads = os.popen("samtools view " + bam_file)
	for read in reads:
		tree = None
		sam = SAM(read)
		seq_id = sam.rname() + "|" + decode(sam.flag())
		if decode(sam.flag()) == "+":
			start, end = int(sam.pos()), int(sam.pos())+1
		else:
			start, end = int(sam.pos())-1, int(sam.pos())
		if seq_id in intervals:
			tree = intervals[seq_id]
		else:
			tree = IntervalTree()
			intervals[seq_id] = tree
		tree.add(start, end, sam)
	logger.info("load_alignments_to_tree(): Loading alignments from: " + bam_file + " into IntervalTree()... [COMPLETE].")
	return intervals

############################
# READ COVERAGE EXTRACTION #
############################
class Coverage(object):
	# Based on the annotation (gene_type, chrom, strand, gene_id, transcript_id, and feature) and coordinates input to
	# this class(), extract the A-site and P-site read coverage (both raw and normalized), the number of A-site reads
	# by read length, and the number of P-site reads by reading frame.
	def __init__(self, bam_tree, asite_buffers, psite_buffers, methods, offsets, transcript_coordinates):
		self.bam_tree = bam_tree
		self.asite_buffers = asite_buffers
		self.psite_buffers = psite_buffers
		self.methods = methods
		self.offsets = offsets
		self.transcript_coordinates = transcript_coordinates

	def annotation(self):
		transcript, coordinates = self.transcript_coordinates
		return transcript

	def asite_length_distribution(self):
		# This function takes as input a tuple of transcript annotation information and its coordinates, and outputs 
		# the A-site adjusted reads within those coordinates organized by read length. Other primary inputs include a
		# BAM file, and a pre-defined set of transcript boundary buffers. The read length distribution is output as a
		# dict() with the read lengths (from 26nt to 34nt) and their respective number of reads aligned to those
		# transcript coordinates.
		def check_length(coords, buffers, feature):
			asite_coordinates = expand_transcript_coordinates(coords)
			return False if len(asite_coordinates) <= sum(buffers[feature]) else True

		if "FLOSS" in self.methods:
			annotation, coordinates = self.transcript_coordinates
			# The annotation_coordinates variable is defined as a tuple containing (annotation, coordinates), where the
			# annotation variable is further broken down into a colon-delimited string of: gene_type, chromosome, 
			# strand, gene_id, transcript_id, and transcript_feature (eg. "CDS", "exon", or "UTR").
			gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
			# UTRs are stored as "UTR5" or "UTR3", and exon boundaries are identical to those defined for CDS regions:
			feature = "UTR" if feature in ("UTR5", "UTR3") else feature
			feature = "CDS" if feature == "exon" else feature
			# Initialize the read length distribution for this transcript:
			read_lengths = {26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0}
			asite_coordinates = expand_transcript_coordinates(coordinates)
			if check_length(coordinates, self.asite_buffers, feature) == True:
				if strand == "+":
					asite_coordinates = asite_coordinates[self.asite_buffers[feature][0]:] if self.asite_buffers[feature][-1] == 0 else asite_coordinates[self.asite_buffers[feature][0]:-self.asite_buffers[feature][-1]]
				else:
					asite_coordinates = asite_coordinates[self.asite_buffers[feature][-1]:] if self.asite_buffers[feature][0] == 0 else asite_coordinates[self.asite_buffers[feature][-1]:-self.asite_buffers[feature][0]]
				buffered_coordinates = regroup_coordinates(asite_coordinates)
				for start, end in buffered_coordinates:
					# Coordinates are buffered by 30nt on both ends to account for offset position of reads:
					reads = self.bam_tree[chrom + "|" + strand].find(start-30, end+30)
					for read in reads:
						offset_position = calculate_offset_position(int(read.pos()), strand, len(read.seq()), read.cigar(), self.offsets["ingolia"])
						if offset_position in range(start, end+1):
							if len(read.seq()) in read_lengths:
								read_lengths[len(read.seq())] += 1
							elif len(read.seq()) <= 26:
								read_lengths[26] += 1
							elif len(read.seq()) >= 34:
								read_lengths[34] += 1
							else:
								pass
				return read_lengths
			else:
				return "NA"
		else:
			return "NA"

	def psite_frame_distribution(self):
		# This function takes as input a BAM alignment file, P-site regional boundary buffers, and an
		# annotation_coordinates() object that defines a transcripts and its constituent coordinates.
		# The transcript is delimited into: gene_type, chromosome, strand, gene, transcript, and feature
		# (CDS, UTR, or Exon). Based on the reads that overlap the buffered transcript coordinates, this
		# function will output the read distribution over those coordinates by reading frame (+0/+1/+2).
		def frame_reads(coverage, read_strand, feature, psite_buffers):
			frame_coverage = list()
			if strand == "+":
				frame_coverage = coverage[psite_buffers[feature][0]:] if psite_buffers[feature][-1] == 0 else coverage[psite_buffers[feature][0]:-psite_buffers[feature][-1]]
			else:
				frame_coverage = coverage[psite_buffers[feature][-1]:] if psite_buffers[feature][0] == 0 else coverage[psite_buffers[feature][-1]:-psite_buffers[feature][0]]
			if sum(frame_coverage) == 0:
				return "NA"
			else:
				masked_region = list()
				for x in xrange(len(frame_coverage)):
					masked_region.append(0) if frame_coverage[x]/math.fsum(frame_coverage) >= 0.7 else masked_region.append(frame_coverage[x])
				return [sum(masked_region[i::3]) for i in (0,1,2)] if strand == "+" else [sum(masked_region[::-1][i::3]) for i in (0,1,2)]	

		if "ORFscore" in self.methods:
			# Transcript is defined as a tuple of (annotation, coordinates):
			annotation, coordinates = self.transcript_coordinates
			# Annotation is further delimited into: gene_type, chromosome, strand, gene, transcript, and feature:
			gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
			feature = "UTR" if feature in ("UTR5", "UTR3") else feature
			feature = "CDS" if feature == "exon" else feature
			transcript_coverage = [0]*transcript_length(coordinates)
			psite_coordinates = expand_transcript_coordinates(coordinates)
			for start, end in coordinates:
				reads = self.bam_tree[chrom + "|" + strand].find(start-30, end+30)
				for read in reads:
					offset_position = calculate_offset_position(int(read.pos()), strand, len(read.seq()), read.cigar(), self.offsets["bazzini"])
					if offset_position in psite_coordinates:
						transcript_coverage[psite_coordinates.index(offset_position)] += 1
			return "NA" if sum(transcript_coverage) == 0 else frame_reads(transcript_coverage, strand, feature, self.psite_buffers)

	def psite_normalized_coverage(self):
		# This function takes as input a BAM alignment file, P-site regional boundary buffers, and an
		# annotation_coordinates() object that defines a transcript and its constituent coordinates.
		# The transcript is de-limited into: gene_type, chromosome, strand, gene, transcrip and feature
		# (CDS, UTR, or Exon). Based on the reads that overlap the buffered transcript coordinates, this
		# function will output the normalized read coverage over those coordinates.
		if "SPECtre" in self.methods:
			annotation, coordinates = self.transcript_coordinates
			gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
			# Instantiate the normalized read coverage for output, and the raw coverage as a list() of length L as
			# defined by the transcript coordinates:
			normalized_coverage = list()
			transcript_coverage = [0]*transcript_length(coordinates)
			# Extract the coordinates into a list():
			transcript_positions = expand_transcript_coordinates(coordinates)
			for start, end in coordinates:
				# Coordinates are buffered by 30nt on both ends to account for offset position of reads:
				reads = self.bam_tree[chrom + "|" + strand].find(start-30, end+30)
				for read in reads:
					if "custom" in self.offsets:
						offset_position = calculate_offset_position(int(read.pos()), strand, len(read.seq()), read.cigar(), self.offsets["custom"])
					else:
						offset_position = calculate_offset_position(int(read.pos()), strand, len(read.seq()), read.cigar(), self.offsets["bazzini"])
					if offset_position in transcript_positions:
						transcript_position = transcript_positions.index(offset_position)
						transcript_coverage[transcript_position] += 1
			if sum(transcript_coverage) == 0:
				return "NA"
			else:
				# Normalize the transcript coverage over each position to the position with the highest read coverage:
				for x in xrange(len(transcript_coverage)):
					normalized_coverage.append(transcript_coverage[x]/float(max(transcript_coverage)))
				return normalized_coverage if strand == "+" else normalized_coverage[::-1]
		else:
			return "NA"

	# This will be deprecrated in future versions, as SPECtre will use P-site adjusted positions.
	def psite_raw_coverage(self):
		# This function takes as input a BAM alignment file, P-site regional boundary buffers, and an
		# annotation_coordinates() object that defines a transcript and its constituent coordinates.
		# The transcript is de-limited into: gene_type, chromosome, strand, gene, transcrip and feature
		# (CDS, UTR, or Exon). Based on the reads that overlap the buffered transcript coordinates, this
		# function will output the raw read coverage over those coordinates.
		if "SPECtre" in self.methods:
			annotation, coordinates = self.transcript_coordinates
			gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
			# Instantiate the raw coverage as a list() of the same length as defined by the transcript coordinates.
			transcript_coverage = [0]*transcript_length(coordinates)
			# Extract the coordinates into a list():
			transcript_positions = expand_transcript_coordinates(coordinates)
			for start, end in coordinates:
				# Coordinates are buffered by 30nt on both ends to account for offset position of reads:
				reads = self.bam_tree[chrom + "|" + strand].find(start-30, end+30)
				for read in reads:
					if "custom" in self.offsets:
						offset_position = calculate_offset_position(int(read.pos()), strand, len(read.seq()), read.cigar(), self.offsets["custom"])
					else:
						offset_position = calculate_offset_position(int(read.pos()), strand, len(read.seq()), read.cigar(), self.offsets["bazzini"])
					if offset_position in transcript_positions:
						transcript_coverage[transcript_positions.index(offset_position)] += 1
			if sum(transcript_coverage) == 0:
				return "NA"
			else:
				return transcript_coverage if strand == "+" else transcript_coverage[::-1]
		else:
			return "NA"

################################
# TRANSCRIPT SCORING FUNCTIONS #
################################
class Coherence(object):

	def __init__(self, psite_buffers, window_length, step_size, spectre_analysis, methods, transcript_coverage):
		self.psite_buffers = psite_buffers
		self.window_length = window_length
		self.step_size = step_size
		self.spectre_analysis = spectre_analysis
		self.methods = methods
		self.transcript_coverage = transcript_coverage

	def annotation(self):
		transcript, psite_coverage = self.transcript_coverage
		return transcript

	def spectre_signal(self):
		annotation, psite_coverage = self.transcript_coverage
		gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
		check = Checks(psite_coverage, "NA", strand, self.window_length, self.psite_buffers.values())
		# Both check.coverage() and check.trimming() must pass:
		if check.coverage() == True and check.trimming() == True:
			# Generate the reference coverage signal based on the length of the normalized region:
			reference_signal = ([4/6.0,1/6.0,1/6.0]*int(math.ceil(len(psite_coverage)/3.0)))[0:len(psite_coverage)]
			coherences = list()
			for i in range(0, len(psite_coverage))[::int(self.step_size)]:
				j = i + self.window_length
				if (math.fsum(psite_coverage[i:j]) == 0) or (len(psite_coverage[i:j]) < self.window_length):
					coherences.append(0.0)
				else:
					r('window.region <- c(%s)' % ",".join(str(n) for n in psite_coverage[i:j]))
					r('window.coding <- c(%s)' % ",".join(str(n) for n in reference_signal[i:j]))
					r('test.spec <- spec.pgram(data.frame(window.region, window.coding), spans=c(3,3), plot=FALSE)')
					coherences.append(r('test.spec$coh[which(abs(test.spec$freq-1/3)==min(abs(test.spec$freq-1/3)))]')[0])
			return coherences
		else:
			return "NA" if check.coverage == True else "NA"

	def spectre_score(self):
		def mean(signal):
			return math.fsum(signal) / float(len(signal))
		def median(signal):
			sorted_signal = sorted(signal)
			midpoint = (len(sorted_signal)-1) // 2
			return sorted_signal[midpoint] if len(sorted_signal) % 2 else math.fsum([sorted_signal[midpoint], sorted_signal[midpoint+1]]) / 2.0
		def maximum(signal):
			return sorted(signal)[-1]
		def nonzero_mean(signal):
			nonzero_scores = [n for n in signal if n > 0]
			return math.fsum(nonzero_scores) / float(len(nonzero_scores))
		def nonzero_median(signal):
			nonzero_scores = sorted([n for n in signal if n > 0])
			midpoint = (len(nonzero_scores)-1) // 2
			return nonzero_scores[midpoint] if len(nonzero_scores) % 2 else math.fsum([nonzero_scores[midpoint], nonzero_scores[midpoint+1]]) / 2.0

		annotation, psite_coverage = self.transcript_coverage
		gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
		check = Checks(psite_coverage, "NA", strand, self.window_length, self.psite_buffers.values())
		# Both check.coverage() and check.trimming() must pass:
		if check.coverage() == True and check.trimming() == True:
			if self.spectre_signal() == "NA":
				# Cannot calculate a SPECtre score:
				return "NA"
			else:
				if self.spectre_analysis == "mean":
					return mean(self.spectre_signal())
				elif self.spectre_analysis == "median":
					return median(self.spectre_signal())
				elif self.spectre_analysis == "maximum":
					return maximum(self.spectre_signal())
				elif self.spectre_analysis == "nonzero_mean":
					return nonzero_mean(self.spectre_signal())
				elif self.spectre_analysis == "nonzero_median":
					return nonzero_median(self.spectre_signal())
				else:
					# Default to calculate the mean of the spectre_signal():
					return mean(self.spectre_signal())
		else:
			return "NA" if check.coverage() else "NA" 

	def spectre_phase(self):
		annotation, psite_coverage = self.transcript_coverage
		gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
		check = Checks(psite_coverage, "NA", strand, self.window_length, self.psite_buffers.values())
		# Both check.coverage() and check.trimming() must pass:
		if check.coverage() == True and check.trimming() == True:
			if self.spectre_signal() == "NA":
				return "NA"
			else:
				# Generate the reference coverage signal based on the length of the normalized region:
				reference_signal = ([4/6.0,1/6.0,1/6.0]*int(math.ceil(len(psite_coverage)/3.0)))[0:len(psite_coverage)]
				r('test.region <- c(%s)' % ",".join(str(n) for n in psite_coverage))
				r('test.coding <- c(%s)' % ",".join(str(n) for n in reference_signal))
				r('spec.coding <- spec.pgram(data.frame(test.region, test.coding), spans=c(3,3), plot=FALSE)')
				frame = r('round(test.spec$phase[which(abs(test.spec$freq-1/3)==min(abs(test.spec$freq-1/3)))])')[0]
				return "+0" if frame == 0 else "+1" if frame == 2 else "+2" if frame == -2 else "NA"
		else:
			return "NA"

	def coherence_signal(self):
		annotation, psite_coverage = self.transcript_coverage
		gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
		return psite_coverage if "Full" in self.methods else "NA"

	def coherence_score(self):
		annotation, psite_coverage = self.transcript_coverage
		gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
		check = Checks(psite_coverage, "NA", strand, self.window_length, self.psite_buffers.values())
		# As long as there is sufficient read coverage, the full coherence should be calculated for transcripts
		# of all lengths:
		if "Full" in self.methods:
			if check.coverage() == True:
				# Load the normalized and reference signals into R and calculate the spectral coherence
				# over the full length of the normalized region:
				reference_signal = ([4/6.0, 1/6.0, 1/6.0]*int(math.ceil(len(psite_coverage)/3.0)))[0:len(psite_coverage)]
				r('test.region <- c(%s)' %",".join(str(n) for n in psite_coverage))
				r('test.coding <- c(%s)' %",".join(str(n) for n in reference_signal))
				r('spec.coding <- spec.pgram(data.frame(test.region, test.coding), spans=c(3,3), plot=FALSE)')
				return r('spec.coding$coh[which(abs(spec.coding$freq-1/3)==min(abs(spec.coding$freq-1/3)))]')[0]
			else:
				return "NA"
		else:
			return "NA"

class FLOSS(object):
	# Calculate the FLOSS distribution for the input transcript or ORF based on the read length distribution:
	def __init__(self, methods, transcript_reads):
		self.methods = methods
		self.transcript_reads = transcript_reads

	def annotation(self):
		transcript, asite_reads = self.transcript_reads
		return transcript

	def distribution(self):
		annotation, asite_reads = self.transcript_reads
		if "FLOSS" in self.methods:
			if asite_reads == "NA":
				return "NA"
			else:
				total_reads = float(sum(asite_reads.values()))
				if total_reads == 0:
					return "NA"
				else:
					for read_length in asite_reads:
						asite_reads[read_length] /= total_reads
					return asite_reads
		else:
			return "NA"

class ORF(object):
	# Calculate the ORFscore based on the enrichment of reads in the first frame of a transcript or ORF relative to the
	# number of reads in the other two frames:
	def __init__(self, psite_buffers, methods, transcript_reads):
		self.psite_buffers = psite_buffers
		self.methods = methods
		self.transcript_reads = transcript_reads

	def annotation(self):
		transcript, psite_reads = self.transcript_reads
		return transcript

	def score(self):
		def frame_score(reads, mean_reads):
			return math.pow((reads - mean_reads), 2) / mean_reads

		annotation, psite_reads = self.transcript_reads
		gene_type, chrom, strand, gene_id, transcript_id, feature = annotation.split(":")
		check = Checks("NA", psite_reads, strand, 0, self.psite_buffers.values())
		if "ORFscore" in self.methods:
			if check.frames() == True:
				frames_mean = sum(psite_reads) / float(len(psite_reads))
				score = math.log(math.fsum([frame_score(psite_reads[0], frames_mean), frame_score(psite_reads[1], frames_mean), frame_score(psite_reads[2], frames_mean), 1]), 2)
				return score if (psite_reads[0] > psite_reads[1]) and (psite_reads[0] > psite_reads[2]) else -score
			else:
				return "NA"
		else:
			return "NA"

#######################################################
# FINAL TRANSCRIPT SCORE CALCULATIONS AND AGGREGATION #
#######################################################
def calculate_transcript_scores(gtf, fpkms, fpkm_cutoff, asite_buffers, psite_buffers, bam_tree, window_length, step_size, spectre_analysis, methods, offsets, threads):

	'''
	For each transcript, this function is to calculate (if so designated) its SPECtre metrics
	(score, and windowed coherence signal), FLOSS metrics (score, and read distribution), and
	ORFscore metrics (score, and distribution of reads over each frame). If directed, the
	SPECtre metrics over a given window length(s) will also be calculated.
	'''

	def calculate_reference_distribution(distributions):
		reference_distribution = {26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0}
		reference_transcripts = len(distributions)
		for dist in distributions:
			for read_length in dist:
				if read_length in reference_distribution:
					reference_distribution[read_length] += dist[read_length]
		for read_length in reference_distribution:
			reference_distribution[read_length] /= reference_transcripts
		return reference_distribution

	def calculate_floss_score(reference_distribution, transcript_score):
		annotation, distribution = transcript_score
		gene_type, chrom, strand, gene_id, transcript_id, feature = annotation.split(":")
		if distribution == "NA":
			return "NA"
		else:
			floss = 0
			for read_length in distribution:
				if read_length in reference_distribution:
					floss += abs(distribution[read_length] - reference_distribution[read_length])
			return floss / 2.0

	def build_translation_distributions(fpkms, fpkm_cutoff, transcript_scores):
		translated, not_translated = list(), list()
		for annotation, score in transcript_scores:
			 gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
			 if not score == "NA":
				 if transcript in fpkms:
				 	translated.append(score) if fpkms[transcript] >= fpkm_cutoff else not_translated.append(score)
		return translated, not_translated

	def calculate_posterior_probability(coding_scores, noncoding_scores, score):
		if score == "NA" or score == 0.0:
			return 0.0
		elif isinstance(score, float) or isinstance(score, int):
			if score > max(coding_scores):
				return 1.0
			else:
				prob_score_coding = len([n for n in coding_scores if n >= score])/float(len(coding_scores))
				prob_coding = len(coding_scores)/float(len(coding_scores)+len(noncoding_scores))
				prob_score = len([n for n in coding_scores+noncoding_scores if n >= score])/float(len(coding_scores)+len(noncoding_scores))
			return (prob_score_coding * prob_coding) / prob_score
		else:
			return 0.0

	def posterior_probability(coding_scores, noncoding_scores, transcript_score):
		annotation, score = transcript_score
		gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
		if isinstance(score, int) or isinstance(score, float):
			return annotation, calculate_posterior_probability(coding_scores, noncoding_scores, score)
		elif isinstance(score, list):
			windowed = list()
			for coh in score:
				windowed.append(calculate_posterior_probability(coding_scores, noncoding_scores, coh))
			return annotation, windowed
		else:
			return annotation, "NA"

	def populate_hash(transcripts, metric, transcript_metrics):
		for annotation, score in transcript_metrics:
			gene_type, chrom, strand, gene, transcript, feature = annotation.split(":")
			transcripts[gene_type][chrom][strand][gene][transcript][feature][metric] = score
		return transcripts

	logger.info("calculate_transcript_scores(): Calculating transcript-level metrics [STARTED].")
	# Instantiate the Pool objects for multiprocessing:
	pool = ThreadPool(threads)
	# Parse the transcript annotations and coordinates into separate lists:
	transcripts, coordinates = zip(*flatten(gtf).iteritems())

	logger.info("calculate_transcript_scores(): Calculating A/P-site read distribution and coverage... [STARTED].")
	##############################################################
	# EXTRACT THE A/P-SITE COVERAGE OVER EACH TRANSCRIPT OR ORF: #
	##############################################################
	# Instantiate the partial parameters to the Coverage() class:
	Coverage_func = partial(Coverage, bam_tree, asite_buffers, psite_buffers, methods, offsets)
	# Calculate the A-site coverage and read distribution, and P-site reading frame coverage, where the results are a
	# Coverage() object for each transcript or ORF:
	coverages = pool.map(Coverage_func, flatten(gtf).iteritems())
	# Partition coverage metrics to their respective containers:
	transcript_asite_distributions = dict(zip([transcript.annotation() for transcript in coverages], [transcript.asite_length_distribution() for transcript in coverages]))
	transcript_psite_distributions = dict(zip([transcript.annotation() for transcript in coverages], [transcript.psite_frame_distribution() for transcript in coverages]))
	transcript_psite_normalized_coverages = dict(zip([transcript.annotation() for transcript in coverages], [transcript.psite_normalized_coverage() for transcript in coverages]))
	transcript_psite_raw_coverages = dict(zip([transcript.annotation() for transcript in coverages], [transcript.psite_raw_coverage() for transcript in coverages]))
	logger.info("calculate_transcript_scores(): Calculating A/P-site read distribution and coverage... [COMPLETE].")

	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence, FLOSS and ORFscore... [STARTED].")
	##################################################################################
	# CALCULATE THE SPECTRE, FULL COHERENCE, FLOSS AND ORFSCORE FOR EACH TRANSCRIPT: #
	##################################################################################
	# Calculate SPECtre and Coherence scores for each transcript:
	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence scores... [STARTED].")
	Coherence_func = partial(Coherence, asite_buffers, window_length, step_size, spectre_analysis, methods)
	coherences = pool.map(Coherence_func, transcript_psite_normalized_coverages.iteritems())
	# Partition coherence scores to their respective containers:
	transcript_spectre_scores = dict(zip([scores.annotation() for scores in coherences], [scores.spectre_score() for scores in coherences]))
	transcript_spectre_signals = dict(zip([scores.annotation() for scores in coherences], [scores.spectre_signal() for scores in coherences]))
	transcript_coherence_scores = dict(zip([scores.annotation() for scores in coherences], [scores.coherence_score() for scores in coherences]))
	transcript_coherence_signals = dict(zip([scores.annotation() for scores in coherences], [scores.coherence_signal() for scores in coherences]))
	transcript_spectre_phases = dict(zip([scores.annotation() for scores in coherences], [scores.spectre_phase() for scores in coherences]))
	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence scores... [COMPLETE].")

	# Calculate FLOSS read distributions for each transcript:
	logger.info("calculate_transcript_scores(): Calculating FLOSS read length distributions... [STARTED].")
	FLOSS_func = partial(FLOSS, methods)
	floss_distributions = pool.map(FLOSS_func, transcript_asite_distributions.iteritems())
	# Parition FLOSS distributions to their respective containers:
	transcript_floss_distributions = dict(zip([floss.annotation() for floss in floss_distributions], [floss.distribution() for floss in floss_distributions]))
	logger.info("calculate_transcript_scores(): Calculating FLOSS read length distributions... [COMPLETE].")

	# Build FLOSS reference distribution:
	logger.info("calculate_transcript_scores(): Calculating FLOSS reference read distribution... [STARTED].")
	protein_coding_distributions = [distribution for transcript, distribution in transcript_floss_distributions.items() if "protein_coding" in transcript.lower() if "cds" in transcript.lower()]
	reference_distribution = calculate_reference_distribution(protein_coding_distributions)
	logger.info("calculate_transcript_scores(): Calculating FLOSS reference read distribution... [COMPLETE].")

	logger.info("calculate_transcript_scores(): Calculating FLOSS from read length distributions... [STARTED].")
	FLOSS_metric = partial(calculate_floss_score, reference_distribution)
	floss_scores = pool.map(FLOSS_metric, transcript_floss_distributions.iteritems())
	transcript_floss_scores = dict(zip(transcripts, floss_scores))
	logger.info("calculate_transcript_scores(): Calculating FLOSS from read length distributions... [COMPLETE].")

	# Calculate ORFscore for each transcript:
	logger.info("calculate_transcript_scores(): Calculating ORFscore from frame distributions... [STARTED].")
	ORF_func = partial(ORF, psite_buffers, methods)
	orf_scores = pool.map(ORF_func, transcript_psite_distributions.iteritems())
	# Parition ORF scores to their respective containers:
	transcript_orf_scores = dict(zip([orf.annotation() for orf in orf_scores], [orf.score() for orf in orf_scores]))
	logger.info("calculate_transcript_scores(): Calculating ORFscore from frame distributions... [COMPLETE].")

	#####################################################################################
	# BUILD DISTRIBUTIONS BASED ON TRANSLATIONAL STATUS FOR PROTEIN-CODING TRANSCRIPTS: #
	#####################################################################################
	logger.info("calculate_transcript_scores(): Building distributions for posteriors calculation... [STARTED].")
	# For transcripts with SPECtre scores:
	logger.info("calculate_transcript_scores(): Buliding distributions for SPECtre posteriors... [STARTED].")
	protein_coding_spectre_scores = [(transcript, score) for transcript, score in transcript_spectre_scores.items() if "protein_coding" in transcript.lower() if "cds" in transcript.lower()]
	translated_spectre_scores, untranslated_spectre_scores = build_translation_distributions(fpkms, fpkm_cutoff, protein_coding_spectre_scores)
	logger.info("calculate_transcript_scores(): Buliding distributions for SPECtre posteriors... [COMPLETE].")

	# For transcripts with Coherence scores:
	logger.info("calculate_transcript_scores(): Building distributions for Coherence posteriors... [STARTED].")
	protein_coding_coherence_scores = [(transcript, score) for transcript, score in transcript_coherence_scores.items() if "protein_coding" in transcript.lower() if "cds" in transcript.lower()]
	translated_coherence_scores, untranslated_coherence_scores = build_translation_distributions(fpkms, fpkm_cutoff, protein_coding_coherence_scores)
	logger.info("calculate_transcript_scores(): Building distributions for Coherence posteriors... [STARTED].")
	logger.info("calculate_transcript_scores(): Building distributions for posteriors calculation... [COMPLETE].")

	##############################################################################
	# CALCULATE THE SPECTRE AND COHERENCE SCORE POSTERIORS FOR EACH TRANSCRIPT : #
	##############################################################################
	logger.info("calculate_transcript_scores(): Calculating transcript posteriors... [STARTED].")
	# Build partial functions:
	Coherence_post = partial(posterior_probability, translated_coherence_scores, untranslated_coherence_scores)
	Windowed_post = partial(posterior_probability, translated_spectre_scores, untranslated_spectre_scores)
	SPECtre_post = partial(posterior_probability, translated_spectre_scores, untranslated_spectre_scores)
	# Calculate the posterior probabilities for each transcript:
	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence transcript posteriors... [STARTED].")
	coherence_posteriors = pool.map(Coherence_post, transcript_coherence_scores.iteritems())
	windowed_posteriors = pool.map(Windowed_post, transcript_spectre_signals.iteritems())
	spectre_posteriors = pool.map(SPECtre_post, transcript_spectre_scores.iteritems())
	# Parition SPECtre/Coherence posteriors to their respective containers:
	transcript_coherence_posteriors = dict(zip([transcript for transcript, posterior in coherence_posteriors], [posterior for transcript, posterior in coherence_posteriors]))
	transcript_windowed_posteriors = dict(zip([transcript for transcript, posterior in windowed_posteriors], [posterior for transcript, posterior in windowed_posteriors]))
	transcript_spectre_posteriors = dict(zip([transcript for transcript, posterior in spectre_posteriors], [posterior for transcript, posterior in spectre_posteriors]))
	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence transcript posteriors... [COMPLETE].")
	logger.info("calculate_transcript_scores(): Calculating transcript-level metrics [COMPLETE].")

	#####################################################################################
	# OUTPUT THE TRANSCRIPT COVERAGES, SCORES AND POSTERIORS TO A COMPOSITE DICTIONARY: #
	#####################################################################################
	logger.info("calculate_transcript_scores(): Output transcript metrics to hash()... [STARTED].")
	# Instantiate the composite transcript hash:
	metrics = hash()

	# Populate the hash() with the A-site and P-site read coverages:
	logger.info("calculate_transcript_scores(): Output transcript A/P-site read coverages... [STARTED].")
	metrics = populate_hash(metrics, "A_distributions", transcript_asite_distributions.items())
	metrics = populate_hash(metrics, "P_distributions", transcript_psite_distributions.items())
	metrics = populate_hash(metrics, "P_normalized", transcript_psite_normalized_coverages.items())
	metrics = populate_hash(metrics, "P_raw", transcript_psite_raw_coverages.items())
	logger.info("calculate_transcript_scores(): Output transcript A/P-site read coverages... [COMPLETE].")

	# Populate the hash() with the SPECtre, Coherence, FLOSS, and ORFscore metrics for each transcript:
	logger.info("calculate_transcript_scores(): Output SPECtre scores, signals and posteriors... [STARTED].")
	metrics = populate_hash(metrics, "SPEC_score", transcript_spectre_scores.items())
	metrics = populate_hash(metrics, "SPEC_signal", transcript_spectre_signals.items())
	metrics = populate_hash(metrics, "SPEC_score_posterior", transcript_spectre_posteriors.items())
	metrics = populate_hash(metrics, "SPEC_signal_posterior", transcript_windowed_posteriors.items())
	metrics = populate_hash(metrics, "SPEC_phase", transcript_spectre_phases.items())
	logger.info("calculate_transcript_scores(): Output SPECtre scores, signals and posteriors... [COMPLETE].")

	if "Full" in methods:
		logger.info("calculate_transcript_scores(): Output Coherence scores, signals and posteriors... [STARTED].")
		metrics = populate_hash(metrics, "FULL_score", transcript_coherence_scores.items())
		metrics = populate_hash(metrics, "FULL_signal", transcript_coherence_signals.items())
		metrics = populate_hash(metrics, "FULL_score_posterior", transcript_coherence_posteriors.items())
		logger.info("calculate_transcript_scores(): Output Coherence scores, signals and posteriors... [COMPLETE].")
	if "FLOSS" in methods:
		logger.info("calculate_transcript_scores(): Output FLOSS scores and distributions... [STARTED].")
		metrics = populate_hash(metrics, "FLOSS_score", transcript_floss_scores.items())
		metrics = populate_hash(metrics, "FLOSS_distribution", transcript_floss_distributions.items())
		logger.info("calculate_transcript_scores(): Output FLOSS scores and distributions... [COMPLETE].")
	if "ORFscore" in methods:
		logger.info("calculate_transcript_scores(): Output ORFscores... [STARTED].")
		metrics = populate_hash(metrics, "ORF_score", transcript_orf_scores.items())
		logger.info("calculate_transcript_scores(): Output ORFscores... [COMPLETE].")
	logger.info("calculate_transcript_scores(): Output transcript metrics to hash()... [COMPLETE].")
	return metrics, reference_distribution

class ExperimentMetrics(object):

	'''
	Once all transcripts have been scored according to the metrics requested, calculate
	experiment-level metrics including thresholds for protein-coding versus non-coding
	SPECtre scores, translated versus untranslated SPECtre scores, and analytical AUCs
	using the ROCR package in R.
	'''

	def __init__(self, stats, fpkms, analyses, cutoff, fdr):
		self.stats = stats
		self.fpkms = fpkms
		self.analyses = analyses
		self.cutoff = cutoff
		self.fdr = fdr

	@staticmethod
	def build_score_distributions(transcript_stats, transcript_fpkms, fpkm_cutoff, method):
		# Further refinement required after checking input format...
		translated, not_translated = list(), list()
		# Get the scores by transcript and method:
		metric = "_".join([method.upper(), "score"])
		for gene_type in transcript_stats:
			if gene_type == "protein_coding":
				for chrom in transcript_stats[gene_type]:
					for strand in transcript_stats[gene_type][chrom]:
						for gene in transcript_stats[gene_type][chrom][strand]:
							for transcript in transcript_stats[gene_type][chrom][strand][gene]:
								for feature in transcript_stats[gene_type][chrom][strand][gene][transcript]:
									if feature == "CDS":
										if metric in transcript_stats[gene_type][chrom][strand][gene][transcript][feature]:
											if not transcript_stats[gene_type][chrom][strand][gene][transcript][feature][metric] == None:
												if not transcript_stats[gene_type][chrom][strand][gene][transcript][feature][metric] == "NA":
													score = transcript_stats[gene_type][chrom][strand][gene][transcript][feature][metric]
													if transcript in transcript_fpkms:
														if transcript_fpkms[transcript] >= fpkm_cutoff:
															translated.append(score)
														else:
															not_translated.append(score)
		return translated, not_translated

	def translation_threshold(self):
		logger.info("ExperimentMetrics.translation_threshold(): Calculating experiment-level translation threshold [STARTED].")
		translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "SPEC")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
		r('scores <- rbind(active, inactive)')
		logger.info("ExperimentMetrics.translation_threshold(): Calculating experiment-level translation threshold [COMPLETE].")
		return str(r('quantile(scores$SPEC[scores$biotype=="not_translated"], probs=%s)' %(1-float(self.fdr)))[0])

	def spectre_auc(self):
		logger.info("ExperimentMetrics.spectre_auc(): Calculating experiment-level SPECtre AUC [STARTED].")
		if "SPECtre" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "SPEC")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.spectre_auc(): Calculating experiment-level SPECtre AUC [COMPLETE].")
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

	def floss_auc(self):
		logger.info("ExperimentMetrics.floss_auc(): Calculating experiment-level FLOSS AUC [STARTED].")
		if "FLOSS" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "FLOSS")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.floss_auc(): Calculating experiment-level FLOSS AUC [COMPLETE].")
			return str(r('performance(prediction(-scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

	def orfscore_auc(self):
		logger.info("ExperimentMetrics.orfscore_auc(): Calculating experiment-level ORFscore AUC [STARTED].")
		if "ORFscore" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "ORF")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.orfscore_auc(): Calculating experiment-level ORFscore AUC [COMPLETE].")
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

	def full_auc(self):
		logger.info("ExperimentMetrics.full_auc(): Calculating experiment-level Spectral Coherence AUC [STARTED].")
		if "Full" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "Full")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.full_auc(): Calculating experiment-level Spectral Coherence AUC [COMPLETE].")
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

def print_metrics(output_file, transcript_stats, experiment_stats, reference_distribution, gtf, fpkms, analyses, parameters, verbose_check):

	def format_coordinates(coords):
		if len(coords) == 0:
			return "NA"
		else:
			return ",".join([str(start) + "-" + str(end) for start, end in coords])

	def return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, feature):
		try:
			return format_coordinates(gtf[gene_type][chrom][strand][gene][transcript][feature])
		except KeyError:
			return "NA"

	def return_metric(stats, gene_type, chrom, strand, gene, transcript, feature, metric):
		if gene_type == "protein_coding":
			try:
				if isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], list):
					return ",".join([str(x) for x in stats[gene_type][chrom][strand][gene][transcript][feature][metric]])
				elif isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], dict):
					if len(stats[gene_type][chrom][strand][gene][transcript][feature][metric]) == 0:
						return "NA"
					else:
						return zip(*stats[gene_type][chrom][strand][gene][transcript][feature][metric].items())
				else:
					return stats[gene_type][chrom][strand][gene][transcript][feature][metric]
			except KeyError:
				return "NA"
		else:
			if feature == "CDS":
				try:
					if isinstance(stats[gene_type][chrom][strand][gene][transcript]["exon"][metric], list):
						return ",".join([str(x) for x in stats[gene_type][chrom][strand][gene][transcript]["exon"][metric]])
					elif isinstance(stats[gene_type][chrom][strand][gene][transcript]["exon"][metric], dict):
						if len(stats[gene_type][chrom][strand][gene][transcript]["exon"][metric]) == 0:
							return "NA"
						else:	
							return zip(*stats[gene_type][chrom][strand][gene][transcript]["exon"][metric].items())
					else:
						return stats[gene_type][chrom][strand][gene][transcript]["exon"][metric]
				except KeyError:
					return "NA"
			else:
				try:
					if isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], list):
						return ",".join([str(x) for x in stats[gene_type][chrom][strand][gene][transcript][feature][metric]])
					elif isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], dict):
						if len(stats[gene_type][chrom][strand][gene][transcript][feature][metric]) == 0:
							return "NA"
						else:
							return zip(*stats[gene_type][chrom][strand][gene][transcript][feature][metric].items())
					else:
						return stats[gene_type][chrom][strand][gene][transcript][feature][metric]
				except KeyError:
					return "NA"

	def format_fpkm(transcript_fpkms, transcript):
		if transcript in transcript_fpkms:
			return str(transcript_fpkms[transcript])
		else:
			return "NA"

	def write_parameters(parameters, analyses):
		parameters_string = "\n# Input = " + str(parameters.input) + " \n# Output = " + str(parameters.output) + "\n# FPKM = " + str(parameters.fpkm) + "\n# GTF = " + str(parameters.gtf)
		parameters_string += "\n# PARAMETERS\n# Analyses = " + str(analyses) + "\n# Window Length = " + str(parameters.len) + "\n# FPKM Cutoff = " + str(parameters.min) + "\n# FDR = " + str(parameters.fdr) + "\n# Test = " + str(parameters.type)
		return parameters_string

	def return_experiment_metric(experiment_stats, analysis):
		try:
			if analysis == "Full":
				return str(experiment_stats.full_auc())
			elif analysis == "FLOSS":
				return str(experiment_stats.floss_auc())
			elif analysis == "ORFscore":
				return str(experiment_stats.orfscore_auc())
		except:
			return "NA"

	def write_experiment_metrics(experiment_stats):
		metric_string = "\n# Translation Threshold = " + str(experiment_stats.translation_threshold()) + "\n# SPECtre AUC = " + str(experiment_stats.spectre_auc()) + "\n# Full AUC = " + return_experiment_metric(experiment_stats, "Full") + "\n# FLOSS AUC = " + return_experiment_metric(experiment_stats, "FLOSS") + "\n# ORFscore AUC = " + return_experiment_metric(experiment_stats, "ORFscore")
		return metric_string

	# Instantiate necessary R packages:
	rocr = importr("ROCR")
	graphics = importr("grDevices")

	logger.info("print_metrics/main(): Output results to file... [STARTED]")
	# Initialize the default header text:
	header = "\nid\tchr\tstrand\tgene_id\ttranscript_id\tgene_type\tribo_fpkm\tcoordinates_5UTR\tcoordinates_CDS\tcoordinates_3UTR"
	for feature in ("UTR5", "CDS", "UTR3"):
		for analysis in analyses:
			if analysis == "SPECtre":
				header += "\t" + "\t".join(["SPEC_score_" + feature, "SPEC_score_posterior_" + feature, "SPEC_signal_" + feature, "SPEC_signal_posterior_" + feature])
			if analysis == "Full":
				header += "\t" + "\t".join(["FULL_score_" + feature, "FULL_score_posterior_" + feature, "FULL_signal_" + feature])
			if analysis == "FLOSS":
				header += "\t" + "\t".join(["FLOSS_score_" + feature, "FLOSS_distribution_" + feature])
			if analysis == "ORFscore":
				header += "\t" + "\t".join(["ORF_score_" + feature, "ORF_reads_" + feature])

	output_file.write("# FILE I/O:" + write_parameters(parameters, analyses))
	output_file.write("\n# METRICS:" + write_experiment_metrics(experiment_stats))
	output_file.write("\n" + str(zip(*reference_distribution.items())))
	output_file.write(header)

	count = 1
	for gene_type in sorted(transcript_stats):
		for chrom in sorted(transcript_stats[gene_type]):
			for strand in transcript_stats[gene_type][chrom]:
				for gene in sorted(transcript_stats[gene_type][chrom][strand]):
					for transcript in sorted(transcript_stats[gene_type][chrom][strand][gene]):
						if gene_type == "protein_coding":
							line = "\n" + "\t".join(str(field) for field in [count, chrom, strand, gene, transcript, gene_type, format_fpkm(fpkms, transcript),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR5"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "CDS"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR3")])
						else:
							line = "\n" + "\t".join(str(field) for field in [count, chrom, strand, gene, transcript, gene_type, format_fpkm(fpkms, transcript),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR5"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "exon"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR3")])							
						for feature in ("UTR5", "CDS", "UTR3"):
							for analysis in analyses:
								if analysis == "SPECtre":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_score_posterior")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_signal")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_signal_posterior"))])
								if analysis == "Full":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FULL_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FULL_score_posterior")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FULL_signal"))])
								if analysis == "FLOSS":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FLOSS_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FLOSS_distribution"))])
								if analysis == "ORFscore":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "ORF_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "P_distributions"))])
						output_file.write(line)
						count += 1
	logger.info("print_metrics/main(): Output results to file... [COMPLETE]")

if __name__ == "__main__":
	
	# Parse arguments from command-line:
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, add_help=False)
	parser.add_argument("--help", action="help", help="show this help message and exit")
	parser.add_argument("--full", action="store_true", default="store_false", help="calculate un-windowed spectral coherence")
	parser.add_argument("--floss", action="store_true", default="store_false", help="calculate FLOSS and distribution")
	parser.add_argument("--orfscore", action="store_true", default="store_false", help="calculate ORFScore and distribution")
	parser.add_argument("--sanitize", action="store_true", default="store_false", help="remove overlapping transcripts")
	spectre_args = parser.add_argument_group("parameters for SPECtre analysis:")
	spectre_args.add_argument("--nt", action="store", required=False, nargs="?", default=1, metavar="INT", help="number of threads for multi-processing (default: %(default)s)")
	spectre_args.add_argument("--len", action="store", required=False, nargs="?", default=30, metavar="INT", help="length of sliding window (default: %(default)s)")
	spectre_args.add_argument("--step", action="store", required=False, nargs="?", default=3, metavar="INT", help="distance between sliding windows (default: %(default)s)")
	spectre_args.add_argument("--min", action="store", required=False, nargs="?", default=3, metavar="FLOAT", help="minimum FPKM for active translation (default: %(default)s FPKM)")
	spectre_args.add_argument("--fdr", action="store", required=False, nargs="?", default=0.05, metavar="FLOAT", help="FDR cutoff (default: %(default)s)")
	spectre_args.add_argument("--type", action="store", required=False, nargs="?", default="median", metavar="TYPE", choices=["mean","median","max","nonzero_mean","nonzero_median"], help="metric for SPECtre analysis (choices: mean,[median],max,nonzero_mean,nonzero_median)")
	spectre_args.add_argument("--offsets", action="store", required=False, nargs="?", default="", metavar="LIST", help="comma-delimited user-defined list of read_length:offset_position definitions (eg. 28:12,29:14,30:15,31:15,32:15)")
	file_args = parser.add_argument_group("input and output parameters:")
	file_args.add_argument("--input", action="store", required=True, nargs="?", metavar="BAM", type=str, help="location of BAM alignment file")
	file_args.add_argument("--output", action="store", required=True, nargs="?", default="./spectre_results.txt", metavar="FILE", help="write results to (default: %(default)s)")
	file_args.add_argument("--fpkm", action="store", required=True, nargs="?", metavar="FILE", type=str, help="location of Cufflinks isoforms.fpkm_tracking file")
	file_args.add_argument("--gtf", action="store", required=True, nargs="?", metavar="FILE", type=str, help="location of GTF annotation file")
	file_args.add_argument("--log", action="store", required=False, nargs="?", default="./spectre_results.log", metavar="FILE", help="track progress to (default: %(default)s)")
	args = parser.parse_args()

	# Enable event logging:
	logger = logging.getLogger(__name__)
	logger.setLevel(logging.INFO)
	handler = logging.FileHandler(args.log)
	handler.setLevel(logging.INFO)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	handler.setFormatter(formatter)
	logger.addHandler(handler)

	# Transcript annotation types were pulled from Ensembl under the entry
	# "Biotype" at http://useast.ensembl.org/Help/Glossary. Mitochondrial and
	# those related to rRNA not included:
	coding = ["protein_coding"]
	noncoding = ["miRNA", "misc_RNA", "ncRNA", "scRNA", "snlRNA", "snoRNA", "snRNA"]
	lncrna = ["lincRNA", "lncRNA"]
	pseudogene = ["process_pseudogene", "pseudogene", "transcribed_processed_pseudogene",
					"transcribed_unprocessed_pseudogene", "translated_processed_pseudogene",
					"translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"]

	# Default transcript boundary buffers:
	asite_buffers = {"CDS": (45,15), "UTR": (15,15)}
	psite_buffers = {"CDS": (3,3), "UTR": (3,3)}

	# Offset positions for A/P sites:
	offsets = {"ingolia": {26: 15, 27: 14, 28: 14, 30: 15, 31: 15}, "bazzini": {26: 12, 27: 12, 28: 12, 29: 12, 30: 12, 31: 12}}
	if len(args.offsets) > 0:
		offsets = add_custom_offsets(args.offsets, offsets)

	# Check BAM input for chromosome format:
	chr_prefix = check_chrom_prefix(args.input)

	# Load BAM into an IntervalTree():
	bam = load_alignments_to_tree(args.input)

	# Extract transcripts, transcript intervals and expression using the provided GTF:
	transcript_fpkms = extract_fpkms(args.fpkm)
	transcript_gtf = parse_gtf(args.gtf, transcript_fpkms, int(args.len), asite_buffers.values() + psite_buffers.values(), chr_prefix, args.sanitize)

	# Initialize the types of analyses to be conducted (default: SPECtre):
	analyses = ["SPECtre"]
	if args.full == True:
		analyses.append("Full")
	if args.floss == True:
		analyses.append("FLOSS")
	if args.orfscore == True:
		analyses.append("ORFscore")

	# Calculate the designated transcript-level scores based on the analyses to be conducted:
	transcript_metrics, reference_read_distribution = calculate_transcript_scores(transcript_gtf, transcript_fpkms, float(args.min), asite_buffers, psite_buffers, bam, int(args.len), int(args.step), args.type, analyses, offsets, int(args.nt))
	
	# Perform a second-pass global analysis based on the transcript-level metrics, such as:
	# ROC analyses, posterior probability as a function of empirical FDR, and codon window
	# signal plots (based on windowed spectral coherence).
	experiment_metrics = ExperimentMetrics(transcript_metrics, transcript_fpkms, analyses, float(args.min), float(args.fdr))
	# Print the results table to the output file:
	print_metrics(open(args.output,"w"), transcript_metrics, experiment_metrics, reference_read_distribution, transcript_gtf, transcript_fpkms, analyses, args, args.verbose)













































