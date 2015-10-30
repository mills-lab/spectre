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

from rpy2.robjects import r
from rpy2.robjects.packages import importr

from bx.intervals.intersection import IntervalTree

def hash():
    return collections.defaultdict(hash)

def convert_chromosome(name):
	return re.sub("chr", "", name)

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

class Annotation(object):

	'''
	This class takes as input a user-supplied GTF transcript annotaiton file,
	and extracts the CDS and UTR coordinates of protein-coding genes, and start
	and end coordinates for non-coding transcripts and loads them into a dict().
	In addition, this function will parse the start and end coordinates for each
	transcript into an IntervalTree(). Default transcript annotation types can
	be found in main(), and can be changed to suit user requirements.
	'''

	def __init__(self, gtf_file, fpkms_file, annotations):
		self.gtf_file = gtf_file
		self.fpkms_file = fpkms_file
		self.annotations = annotations

	def anno(self):

		def append_stop_coordinates(gtf):
			# Since Ensembl annotation GTFs have separate entries for stop coordinates, these
			# must be manually appended to the last exon of each protein-coding gene.
			for gene_type in gtf:
				if gene_type == "protein_coding":
					for chr in gtf[gene_type]:
						for strand in gtf[gene_type][chr]:
							for gene in gtf[gene_type][chr][strand]:
								for transcript in gtf[gene_type][chr][strand][gene]:
									for feature in gtf[gene_type][chr][strand][gene][transcript]:
										if feature == "CDS":
											if strand == "+":
												gtf[gene_type][chr][strand][gene][transcript][feature][-1] = (gtf[gene_type][chr][strand][gene][transcript][feature][-1][0],
																												gtf[gene_type][chr][strand][gene][transcript][feature][-1][-1]+3)
											else:
												gtf[gene_type][chr][strand][gene][transcript][feature][0] = (gtf[gene_type][chr][strand][gene][transcript][feature][0][0]-3,
																												gtf[gene_type][chr][strand][gene][transcript][feature][0][-1])
			return gtf

		def partition_utr_coordinates(gtf):
			# Since Ensembl does not differentiate between the 5'UTR and 3'UTR, they must be
			# annotated based on their position relative to the CDS and strand:
			for gene_type in gtf:
				if gene_type == "protein_coding":
					for chr in gtf[gene_type]:
						for strand in gtf[gene_type][chr]:
							for gene in gtf[gene_type][chr][strand]:
								for transcript in gtf[gene_type][chr][strand][gene]:
									# Only process protein-coding transcripts with annotated CDS and UTRs:
									if all (features in gtf[gene_type][chr][strand][gene][transcript] for features in ("CDS", "UTR")):
										utr_coordinates = gtf[gene_type][chr][strand][gene][transcript]["UTR"]
										five_prime, three_prime = list(), list()
										if strand == "+":
											cds_start = gtf[gene_type][chr][strand][gene][transcript]["CDS"][0][0]
											for coordinate in utr_coordinates:
												if coordinate[-1] <= cds_start:
													five_prime.append(coordinate)
												else:
													three_prime.append(coordinate)
										else:
											cds_start = gtf[gene_type][chr][strand][gene][transcript]["CDS"][-1][-1]
											for coordinate in utr_coordinates:
												if coordinate[0] >= cds_start:
													five_prime.append(coordinate)
												else:
													three_prime.append(coordinate)
										# Remove the previous entry for UTR:
										del gtf[gene_type][chr][strand][gene][transcript]["UTR"]
										# Add the partitioned UTR coordinates to the transcript entry:
										gtf[gene_type][chr][strand][gene][transcript]["UTR5"] = five_prime
										gtf[gene_type][chr][strand][gene][transcript]["UTR3"] = three_prime
			return gtf

		# Load transcripts from the GTF file into a dict():
		logger.info("(in anno) Parsing transcripts from GTF: " + self.gtf_file + " into memory... [STARTED]")
		transcripts = dict()
		for line in open(self.gtf_file):
			if not line.startswith("#"):
				seq_name, source, feature, start, end, score, strand, frame, attributes = line.strip().split("\t")
				seq_name = convert_chromosome(seq_name)
				gene_type = re.findall("gene_biotype .*;{1}", attributes)[0].split('"')[1]
				# Parse annotated protein-coding CDS into the GTF dictionary:
				if (gene_type == "protein_coding" and feature in ("CDS", "UTR")) or (gene_type != "protein_coding" and feature == "exon"):
					gene, transcript = re.findall('ENS[A-Z]*G[0-9]*', attributes)[0], re.findall('ENS[A-Z]*T[0-9]*', attributes)[0]
					if gene_type in transcripts:
						if seq_name in transcripts[gene_type]:
							if strand in transcripts[gene_type][seq_name]:
								if gene in transcripts[gene_type][seq_name][strand]:
									if transcript in transcripts[gene_type][seq_name][strand][gene]:
										if feature in transcripts[gene_type][seq_name][strand][gene][transcript]:
											transcripts[gene_type][seq_name][strand][gene][transcript][feature].append((int(start), int(end)))
										else:
											transcripts[gene_type][seq_name][strand][gene][transcript][feature] = [(int(start), int(end))]
									else:
										transcripts[gene_type][seq_name][strand][gene][transcript] = {feature: [(int(start), int(end))]}
								else:
									transcripts[gene_type][seq_name][strand][gene] = {transcript: {feature: [(int(start), int(end))]}}
							else:
								transcripts[gene_type][seq_name][strand] = {gene: {transcript: {feature: [(int(start), int(end))]}}}
						else:
							transcripts[gene_type][seq_name] = {strand: {gene: {transcript: {feature: [(int(start), int(end))]}}}}
					else:
						transcripts[gene_type] = {seq_name: {strand: {gene: {transcript: {feature: [(int(start), int(end))]}}}}}
		transcripts = partition_utr_coordinates(append_stop_coordinates(transcripts))
		logger.info("(in anno) Parsing transcript coordinates from GTF: " + self.gtf_file + " into memory... [COMPLETE]")
		return transcripts

	def interval_tree(self):
		# Load the start and end coordinates for each transcript into an IntervalTree()
		# object. A 5kb buffer is added to each end per Ingolia, et al. (2014) in order
		# to help filter out proximal transcripts. This buffer can be changed or removed
		# according to user requirements, but is not recommended.
		logger.info("(in interval_tree) Parsing transcript intervals from GTF into an IntervalTree()... [STARTED]")
		intervals = dict()
		gtf = self.anno()
		for gene_type in gtf:
			for chr in gtf[gene_type]:
				for strand in gtf[gene_type][chr]:
					for gene in gtf[gene_type][chr][strand]:
						for transcript in gtf[gene_type][chr][strand][gene]:
							for feature in gtf[gene_type][chr][strand][gene][transcript]:
								# For protein-coding genes use the annotated CDS start and end coordinates,
								# for non-coding transcripts use the annotated transcript start and end
								# coordinates.
								if feature in ("CDS", "exon"):
									tree = None
									seq_id = chr + "|" + strand
									start = gtf[gene_type][chr][strand][gene][transcript][feature][0][0] - 5000
									end = gtf[gene_type][chr][strand][gene][transcript][feature][-1][-1] + 5000
									if seq_id in intervals:
										tree = intervals[seq_id]
									else:
										tree = IntervalTree()
										intervals[seq_id] = tree
									tree.add(start, end, (gene, transcript))
		logger.info("(in interval_tree) Parsing transcript intervals from GTF into an IntervalTree()... [COMPLETE]")
		return intervals

	def fpkms(self):
		# This function takes as input a Cufflinks isoforms.fpkm_tracking output file and
		# extracts the transcript-level expression in FPKM. Pleast note that, as written,
		# this function should support a generic tab-delimited file such that the name of
		# the transcript is in the first column, and the column holding the expression
		# values is titled as "FPKM".
		logger.info("(in fpkms) Parsing transcript FPKMs from file: " + self.fpkms_file + " into memory... [STARTED]")
		fpkm = dict()
		position = 0
		for linenum, line in enumerate(open(self.fpkms_file)):
			if linenum == 0:
				for i, j in enumerate(line.strip().split("\t")):
					if j.upper() == "FPKM":
						position = i
			else:
				name, exp = line.strip().split("\t")[0], float(line.strip().split("\t")[position])
				if name not in fpkm:
					fpkm[name] = exp
		logger.info("(in fpkms) Parsing transcript FPKMs from file: " + self.fpkms_file + " into memory.... [COMPLETE]")
		return fpkm

def remove_overlapping_transcripts(gtf, transcript_tracker):

	'''
	This function takes as input a transcript GTF dictionary and a transcript IntervalTree()
	object, and removes from the transcript GTF dictionary those transcripts that overlap an
	annotated transcript in the IntervalTree() object.
	'''

	def count_transcripts(gtf):
		count = 0
		for gene_type in gtf:
			for chr in gtf[gene_type]:
				for strand in gtf[gene_type][chr]:
					for gene in gtf[gene_type][chr][strand]:
						count += len(gtf[gene_type][chr][strand][gene])
		return count

	#logger.info("Number of transcripts prior to sanitizing: ", str(count_transcripts(gtf)))
	cleaned_gtf = hash()
	count = 0
	logger.info("(in remove_overlapping_transcipts) Sanitizing transcript database [STARTED].")
	for gene_type in gtf:
		for chr in gtf[gene_type]:
			for strand in gtf[gene_type][chr]:
				for gene in gtf[gene_type][chr][strand]:
					for transcript in gtf[gene_type][chr][strand][gene]:
						# Use the CDS or transcript sequence as the test region:
						for feature in gtf[gene_type][chr][strand][gene][transcript]:
							if feature in ("CDS", "exon"):
								# Add a 5000 nt buffer to each transcript start and end coordinate:
								seq_id = chr + "|" + strand
								start, end = gtf[gene_type][chr][strand][gene][transcript][feature][0][0] - 5000, gtf[gene_type][chr][strand][gene][transcript][feature][-1][-1] + 5000
								overlaps = list()
								# Test for overlap with other transcripts:
								if seq_id in transcript_tracker:
									overlaps.extend(transcript_tracker[seq_id].find(start, end))
								# Transcript is allowed to overlap with itself, therefore more than two (2) overlaps
								# should be flagged for removal:
								if len(overlaps) > 1:
									if gene_type not in cleaned_gtf:
										if chr not in cleaned_gtf[gene_type]:
											if strand not in cleaned_gtf[gene_type][chr]:
												if gene not in cleaned_gtf[gene_type][chr][strand]:
													if transcript not in cleaned_gtf[gene_type][chr][strand][gene]:
														cleaned_gtf[gene_type][chr][strand][gene][transcript] = gtf[gene_type][chr][strand][gene][transcript]
														count += 1
	#logger.info("Number of transcripts after sanitizing: ", str(count))
	logger.info("(in remove_overlapping_transcipts) Sanitizing transcript database [COMPLETE].")
	return cleaned_gtf

class Alignment(object):

	'''
	This class object takes as input a BAM alignment file, and a set of offsets to adjust
	the 5' aligned position of the read to its A-site or P-site equivalent as protected
	by a ribosome. Read position adjustments are made according to the CIGAR string
	field of the alignment file.
	'''

	# CIGAR string operations:
	# M 0 alignment match (can be a sequence match or mismatch)
	# I 1 insertion to the reference
	# D 2 deletion from the reference
	# N 3 skipped region from the reference
	# S 4 soft clipping (clipped sequences present in SEQ)
	# H 5 hard clipping (clipped sequences NOT present in SEQ)
	# P 6 padding (silent deletion from padded reference)
	# = 7 sequence match
	# X 8 sequence mismatch

	def __init__(self, bam_file, offset_method):
		self.bam_file = bam_file
		self.offset_method = offset_method

	def coverage(self):

		def decode(flag):
			if int(flag) & 16:
				return "-"
			else:
				return "+"

		def calculate_offset(offsets, read_length):
			if read_length in offsets:
				return offsets[read_length]
			else:
				if len(set(offsets.values())) == 1:
					# All Bazzini offset values are the same:
					return offsets.values()[0]
				else:
					# Ingolia offsets are approximately half the read length:
					return read_length / 2

		def extract_coordinates_from_cigar(read_position, cigar):
			# Return a list of expanded coordinates based on the read position and its CIGAR string:
			coordinates = list()
			ops = re.findall("[0-9]*[DHIMNPSX=]{1}", cigar)
			for op in ops:
				increment, modifier = int(op[:-1]), op[-1]
				if not modifier == "N":
					increments = range(increment)
					for x in xrange(len(increments)):
						coordinates.append(read_position + increments[x])
				else:
					position = coordinates[-1] + increment
			return coordinates

		def offset_read_position(pos, strand, length, cigar, offsets):
			if strand == "+":
				return extract_coordinates_from_cigar(pos, cigar)[calculate_offset(offsets, length)-1]
			else:
				return extract_coordinates_from_cigar(pos, cigar)[-calculate_offset(offsets, length)]

		offsets = dict()
		if self.offset_method == "ingolia":
			offsets = {26: 14, 27: 14, 28: 14, 29: 14, 30: 15, 31: 15}
		else:
			offsets = {26: 12, 27: 12, 28: 12, 29: 12, 30: 12, 31: 12}

		logger.info("LOADING " + self.offset_method + " reads: STARTED")	
		logger.info("(in Alignment_coverage) Loading BAM/SAM alignment file: " + self.bam_file + " into coverage [STARTED].")
		# Load reads from BAM alignment file into a coverage IntervalTree() object:
		coverage = dict()
		reads = os.popen("samtools " + "view " + self.bam_file).readlines()
		for read in reads:
			sam = SAM(read)
			# Convert SAM field values into required input formats for coverage functions:
			strand = decode(sam.flag())
			offset_position = offset_read_position(int(sam.pos()), strand, len(sam.seq()), sam.cigar(), offsets)
			# Load the read into the IntervalTree():
			tree = None
			seq_id = sam.rname() + "|" + strand
			if seq_id in coverage:
				tree = coverage[seq_id]
			else:
				tree = IntervalTree()
				coverage[seq_id] = tree
			tree.add(offset_position, offset_position, read)
		logger.info("(in Alignment_coverage) Loading BAM/SAM alignment file: " + self.bam_file + " into coverage [COMPLETE].")
		logger.info("LOADING " + self.offset_method + " reads: COMPLETE")	
		return coverage

class Coverage(object):

	'''
	Read coverage metrics over the range of a transcript.

	Attributes:
		chr: Chromosome ID in Ensembl format.
		strand: The strand of the region to be tested.
		coordinates: The coordinates over which to extract read coverage.
		buffer: Trimming buffer for transcript boundaries.
		coverage_tracker: An IntervalTree() coverage object extract from a BAM/SAM.
	'''

	def __init__(self, chr, strand, coordinates, asite_buffer, psite_buffer, asite_tracker, psite_tracker):
		self.chr = chr
		self.strand = strand
		self.coordinates = coordinates
		self.asite_buffer = asite_buffer
		self.psite_buffer = psite_buffer
		self.asite_tracker = asite_tracker
		self.psite_tracker = psite_tracker

	def normalized_asite_coverage(self):
		'''
		There is a known bug/feature of the find() function in bx.intervals.intersection,
		such that the boundaries of a given region are not included in the search, see
		below code snippet for example:

		from bx.intervals.intersection import Intersector, Interval
		tree = Intersector()
		tree.add_interval(Interval(1,10))
		tree.find(1,1)
		[]
		tree.find(2,2)
		[Interval(1,10)]

		Thus, to find the single nucleotide position of reads that overlap a given position,
		the find() function is to be implemented as follows:

		tree.add_interval(Interval(1,1))
		tree.add_interval(Interval(2,2))
		tree.add_interval(Interval(3,3))
		tree.find(0,2)
		[Interval(1,1)]
		tree.find(1,3)
		[Interval(2,2)]
		tree.find(2,4)
		[Interval(3,3)]
		'''
		if not self.asite_tracker:
			return "NA"
		else:
			region = list()
			for coordinate in self.coordinates:
				# Since Python indexes are 0-based:
				start, end = int(coordinate[0]), int(coordinate[-1])+1
				for position in range(start, end):
					seq_id = self.chr + "|" + self.strand
					if seq_id in self.asite_tracker:
						overlaps = list()
						overlaps.extend(self.asite_tracker[seq_id].find(position-1, position+1))
						region.append(len(overlaps))
					else:
						region.append(0)
			if len(region) == 0:
				return "NA"
			else:
				if max(region) == 0:
					return "NA"
				else:
					normalized_region = list()
					for x in xrange(len(region)):
						normalized_region.append(region[x]/float(max(region)))
					return [normalized_region[::-1], normalized_region][self.strand == "+"]

	def asite_reads(self):
		'''
		Given the chromosome, strand, coordinates, and the FLOSS-based regional buffer,
		return the A-site reads over the designated region.
		'''
		if not(self.asite_tracker):
			return "NA"
		else:
			reads = dict()
			expanded_coordinates = sum([range(start, end+1) for start, end in self.coordinates], [])
			# Trim the region according to the provided buffer:
			if self.asite_buffer[0] == self.asite_buffer[1]:
				if self.asite_buffer[0] == self.asite_buffer[1] == 0:
					pass
				else:
					expanded_coordinates = expanded_coordinates[self.asite_buffer[0]-1:self.asite_buffer[1]]
			else:
				if self.asite_buffer[0] == 0:
					if self.strand == "+":
						expanded_coordinates = expanded_coordinates[0:self.asite_buffer[1]]
					else:
						expanded_coordinates = expanded_coordinates[self.asite_buffer[1]-1:]
				else:
					if self.strand == "+":
						expanded_coordinates = expanded_coordinates[self.asite_buffer[0]-1:]
					else:
						expanded_coordinates = expanded_coordinates[0:self.asite_buffer[0]]
			for position in expanded_coordinates:
				seq_id = self.chr + "|" + self.strand
				if seq_id in self.asite_tracker:
					overlaps = list()
					overlaps.extend(self.asite_tracker[seq_id].find(position-1, position+1))
					for read in overlaps:
						read_length = len(read.split("\t")[9])
						if read_length in reads:
							reads[read_length] += 1
						else:
							reads[read_length] = 1
			return reads

	def psite_frames(self):
		'''
		Given the chromosome, strand, coordinates, and the ORFscore-based regional buffer,
		return the P-site reads over the designated region.
		'''
		if not(self.psite_tracker):
			return "NA"
		else:
			frames = list()
			expanded_coordinates = sum([range(start, end+1) for start, end in self.coordinates], [])
			# Trim the region according to the provided buffer:
			if self.psite_buffer[0] == self.psite_buffer[1]:
				if self.psite_buffer[0] == self.psite_buffer[1] == 0:
					pass
				else:
					expanded_coordinates = expanded_coordinates[self.psite_buffer[0]-1:-self.psite_buffer[1]]
			else:
				if self.psite_buffer[0] == 0:
					if self.strand == "+":
						expanded_coordinates = expanded_coordinates[0:-self.psite_buffer[1]]
					else:
						expanded_coordinates = expanded_coordinates[self.psite_buffer[1]-1:]
				else:
					if self.strand == "+":
						expanded_coordinates = expanded_coordinates[self.psite_buffer[0]-1:]
					else:
						expanded_coordinates = expanded_coordinates[0:-self.psite_buffer[0]]
			# Collect the read depth over the trimmed coordinates:
			region = list()
			for position in expanded_coordinates:
				seq_id = self.chr + "|" + self.strand
				if seq_id in self.psite_tracker:
					overlaps = list()
					overlaps.extend(self.psite_tracker[seq_id].find(position-1, position+1))
					region.append(len(overlaps))
				else:
					region.append(0)
			# Mask putative peaks in the transcript:
			if sum(region) == 0:
				frames = "NA"
			else:
				masked_region = list()
				for x in xrange(len(region)):
					if region[x]/math.fsum(region) >= 0.7:
						masked_region.append(0)
					else:
						masked_region.append(region[x])
				if self.strand == "+":
					frames = [sum(masked_region[i::3]) for i in (0,1,2)]
				else:
					frames = [sum(masked_region[::-1][i::3]) for i in (0,1,2)]
			return frames

class ORF(object):

	def __init__(self, frame_reads):
		self.frame_reads = frame_reads

	def score(self):

		def frame_score(reads, mean_reads):
			return math.pow((reads - mean_reads), 2) / mean_reads
		frames_mean = sum(self.frame_reads) / float(len(self.frame_reads))
		orf_score = math.log(math.fsum([frame_score(self.frame_reads[0], frames_mean),
										frame_score(self.frame_reads[1], frames_mean),
										frame_score(self.frame_reads[2], frames_mean), 1]), 2)
		if (self.frame_reads[0] > self.frame_reads[1]) and (self.frame_reads[0] > self.frame_reads[2]):
			return orf_score
		else:
			return -orf_score

class FLOSS(object):

	def __init__(self, region_reads, reference_distribution):
		self.region_reads = region_reads
		self.reference_distribution = reference_distribution

	def distribution(self):
		'''
		Calculate the read length distribution over a region by summing the number of reads
		of a given length contained in that region. Reads less than or equal to 26 nucleotides
		in length are to be grouped together, and reads longer than or equal to 34 nucleotides
		in length are to be grouped together.
		'''
		dist = {26: 0.0, 27: 0.0, 28: 0.0, 29: 0.0, 30: 0.0, 31: 0.0, 32: 0.0, 33: 0.0, 34: 0.0}
		for read_length in self.region_reads:
			if read_length <= 26:
				dist[26] += self.region_reads[read_length]/float(sum(self.region_reads.values()))
			elif read_length >= 34:
				dist[34] += self.region_reads[read_length]/float(sum(self.region_reads.values()))
			else:
				dist[read_length] += self.region_reads[read_length]/float(sum(self.region_reads.values()))
		return dist

	def score(self):
		floss = 0
		transcript_distribution = self.distribution()
		for read_length in transcript_distribution:
			if read_length in self.reference_distribution:
				floss += abs(transcript_distribution[read_length] - self.reference_distribution[read_length])
		return floss / 2.0

class SPECtre(object):

	def __init__(self, normalized_coverage, window_length, metric):
		self.normalized_coverage = normalized_coverage
		self.window_length = int(window_length)
		self.metric = metric

		'''
		Given the normalized read coverage over a region, the window_length to be tested, a
		minimum fpkm_cutoff, statistical test, and the desired FDR, calculate the spectral
		coherence and Bayesian posterior probability over the region. If the window_length
		variable is zero, this indicates that the full transcript spectral coherence is to be
		calculated.
		'''

	def signal(self):
		if self.window_length == 0:
			# Since no windows are moved along the region, there is only a single data point: the score().
			return self.normalized_coverage
		else:
			# Generate the reference coverage signal based on the length of the normalized region:
			reference_signal = ([4/6.0, 1/6.0, 1/6.0]*int(math.ceil(len(self.normalized_coverage)/3.0)))[0:len(self.normalized_coverage)]
			coherences = list()
			if len(self.normalized_coverage) >= self.window_length:
				for i in range(0, len(self.normalized_coverage))[::3]:
					j = i + self.window_length
					if (math.fsum(self.normalized_coverage[i:j]) == 0) or (len(self.normalized_coverage[i:j]) < self.window_length):
						coherences.append(0.0)
					else:
						r('window.region <- c(%s)' % ",".join(str(n) for n in self.normalized_coverage[i:j]))
						r('window.coding <- c(%s)' % ",".join(str(n) for n in reference_signal[i:j]))
						r('test.spec <- spec.pgram(data.frame(window.region, window.coding), spans=c(3,3), plot=FALSE)')
						coherences.append(r('test.spec$coh[which(abs(test.spec$freq-1/3)==min(abs(test.spec$freq-1/3)))]')[0])
			return coherences

	def score(self):
		if self.window_length == 0:
			# Load the normalized and reference signals into R and calculate the spectral coherence
			# over the full length of the normalized region:
			reference_signal = ([4/6.0, 1/6.0, 1/6.0]*int(math.ceil(len(self.normalized_coverage)/3.0)))[0:len(self.normalized_coverage)]
			r('test.region <- c(%s)' %",".join(str(n) for n in self.normalized_coverage))
			r('test.coding <- c(%s)' %",".join(str(n) for n in reference_signal))
			r('spec.coding <- spec.pgram(data.frame(test.region, test.coding), spans=c(3,3), plot=FALSE)')
			return r('spec.coding$coh[which(abs(spec.coding$freq-1/3)==min(abs(spec.coding$freq-1/3)))]')[0]
		else:
			if self.metric == "mean":
				return math.fsum(self.signal()) / float(len(self.signal()))
			elif self.metric == "median":
				sorted_signal = sorted(self.signal())
				midpoint = (len(sorted_signal)-1) // 2
				if len(sorted_signal) % 2:
					return sorted_signal[midpoint]
				else:
					return math.fsum([sorted_signal[midpoint], sorted_signal[midpoint+1]]) / 2.0
			elif self.metirc == "maximum":
				sorted_signal = sorted(signal)
				return sorted_signal[-1]
			elif self.metric == "nonzero_mean":
				nonzero_score = [n for n in signal if n > 0]
				return math.fsum(nonzero_score) / float(len(norzero_score))
			elif self.metric == "nonzero_median":
				sorted_signal = sorted([n for n in self.signal() if n > 0])
				midpoint = (len(sorted_signal)-1) // 2
				if len(sorted_signal) % 2:
					return sorted_signal[midpoint]
				else:
					return math.fsum([sorted_signal[midpoint], sorted_signal[midpoint+1]]) / 2.0

class Reference(object):

	'''
	Given a set of protein-coding transcripts, calculate the averaged read length distribution
	for each transcript.
	'''

	def __init__(self, gtf, asite_tracker, buffers):
		self.gtf = gtf
		self.asite_tracker = asite_tracker
		self.asite_buffers = asite_buffers

	def distribution(self):
		read = {26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0}
		transcript_count = 0.0
		for chr in self.gtf:
			for strand in self.gtf[chr]:
				for gene in self.gtf[chr][strand]:
					for transcript in self.gtf[chr][strand][gene]:
						for feature in self.gtf[chr][strand][gene][transcript]:
							if feature == "CDS":
								region = Coverage(chr, strand, self.gtf[chr][strand][gene][transcript][feature], asite_buffers["CDS"], (0,0), self.asite_tracker, dict())
								if not region.asite_reads():
									pass
								else:
									read[26] += sum(num for read_length, num in region.asite_reads().items() if read_length <= 26) / float(sum(region.asite_reads().values()))
									read[27] += sum(num for read_length, num in region.asite_reads().items() if read_length == 27) / float(sum(region.asite_reads().values()))
									read[28] += sum(num for read_length, num in region.asite_reads().items() if read_length == 28) / float(sum(region.asite_reads().values()))
									read[29] += sum(num for read_length, num in region.asite_reads().items() if read_length == 29) / float(sum(region.asite_reads().values()))
									read[30] += sum(num for read_length, num in region.asite_reads().items() if read_length == 30) / float(sum(region.asite_reads().values()))
									read[31] += sum(num for read_length, num in region.asite_reads().items() if read_length == 31) / float(sum(region.asite_reads().values()))
									read[32] += sum(num for read_length, num in region.asite_reads().items() if read_length == 32) / float(sum(region.asite_reads().values()))
									read[33] += sum(num for read_length, num in region.asite_reads().items() if read_length == 33) / float(sum(region.asite_reads().values()))
									read[34] += sum(num for read_length, num in region.asite_reads().items() if read_length >= 34) / float(sum(region.asite_reads().values()))
									transcript_count += 1.0
		if transcript_count == 0:
			return read
		else:			
			for length in read:
				read[length] /= transcript_count
			return read

class Checks(object):

	def __init__(self, region, frame_reads, strand, length, buffers, fpkms, transcript):
		self.region = region
		self.frame_reads = frame_reads
		self.strand = strand
		self.length = length
		self.buffers = buffers
		self.fpkms = fpkms
		self.transcript = transcript

	def coverage(self):
		max_left_trim, max_right_trim = max([i[0] for i in self.buffers]), max([j[-1] for j in self.buffers])
		if self.region == "NA":
			return False
		else:
			if self.strand == "+":
				if math.fsum(self.region[max_left_trim-1:-max_right_trim]) > 0.0:
					return True
				else:
					return False
			else:
				if math.fsum(self.region[max_right_trim-1:-max_left_trim]) > 0.0:
					return True
				else:
					return False

	def frames(self):
		if self.frame_reads == "NA":
			return False
		elif sum(self.frame_reads) > 0:
			return True
		else:
			return False

	def fpkm(self):
		# Since transcript identifiers are unique:
		if self.transcript in self.fpkms:
			if self.fpkms[self.transcript] > 0:
				return True
			else:
				return False
		else:
			return False

	def trimming(self):
		# Require that each post-trimmed transcript region be at least the length of the sliding window:
		max_left_trim, max_right_trim = max([i[0] for i in self.buffers]), max([j[-1] for j in self.buffers])
		if len(self.region) - max_left_trim - max_right_trim >= int(self.length):
			return True
		else:
			return False

def build_translated_distributions(transcript_scores, transcript_fpkms, fpkm_cutoff, analysis):
	translated_scores = list()
	untranslated_scores = list()
	for gene_type in transcript_scores:
		if gene_type == "protein_coding":
			for chr in transcript_scores[gene_type]:
				for strand in transcript_scores[gene_type][chr]:
					for gene in transcript_scores[gene_type][chr][strand]:
						for transcript in transcript_scores[gene_type][chr][strand][gene]:
							for feature in transcript_scores[gene_type][chr][strand][gene][transcript]:
								if feature == "CDS":
									if analysis == "SPEC":
										if "SPEC" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]:
											if "score_windowed" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]:
												if transcript in transcript_fpkms:
													if transcript_fpkms[transcript] >= float(fpkm_cutoff):
														translated_scores.append(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_windowed"])
													else:
														untranslated_scores.append(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_windowed"])
									elif analysis == "Full":
										if "SPEC" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]:
											if "score_full" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]:
												if transcript in transcript_fpkms:
													if transcript_fpkms[transcript] >= float(fpkm_cutoff):
														translated_scores.append(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full"])
													else:
														untranslated_scores.append(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full"])
									else:
										if "score" in transcript_scores[gene_type][chr][strand][gene][transcript][feature][analysis]:
											if transcript in transcript_fpkms:
												if transcript_fpkms[transcript] >= float(fpkm_cutoff):
													translated_scores.append(transcript_scores[gene_type][chr][strand][gene][transcript][feature][analysis]["score"])
												else:
													untranslated_scores.append(transcript_scores[gene_type][chr][strand][gene][transcript][feature][analysis]["score"])
	return translated_scores, untranslated_scores

def calculate_transcript_metrics(gtf, fpkms, asite_buffer, asite_tracker, psite_buffer, psite_tracker, cutoff, window_length, reference_distribution, spectre_analysis, methods):

	'''
	For each transcript, this function is to calculate (if so designated) its SPECtre metrics
	(score, and windowed coherence signal), FLOSS metrics (score, and read distribution), and
	ORFscore metrics (score, and distribution of reads over each frame). If directed, the
	SPECtre metrics over a given window length(s) will also be calculated.
	'''

	def calculate_transcript_scores(gtf, fpkms, asite_buffer, asite_tracker, psite_buffer, psite_tracker, window_length, reference_distribution, spectre_analsyis, methods):
		logger.info("(in calculate_transcript_scores) Calculating transcript-level metrics [STARTED].")
		scores = dict()
		for gene_type in gtf:
			for chr in gtf[gene_type]:
				for strand in gtf[gene_type][chr]:
					for gene in gtf[gene_type][chr][strand]:
						for transcript in gtf[gene_type][chr][strand][gene]:
							for feature in gtf[gene_type][chr][strand][gene][transcript]:
								if feature in ("CDS", "exon"):
									region = Coverage(chr, strand, gtf[gene_type][chr][strand][gene][transcript][feature], asite_buffer["CDS"], psite_buffer["CDS"], asite_tracker, psite_tracker)
									check = Checks(region.normalized_asite_coverage(), region.psite_frames(), strand, window_length, (asite_buffer["CDS"], psite_buffer["CDS"]), fpkms, transcript)
								else:
									region = Coverage(chr, strand, gtf[gene_type][chr][strand][gene][transcript][feature], asite_buffer["UTR"], psite_buffer["UTR"], asite_tracker, psite_tracker)
									check = Checks(region.normalized_asite_coverage(), region.psite_frames(), strand, window_length, (asite_buffer["UTR"], psite_buffer["UTR"]), fpkms, transcript)
								if (check.trimming() == True) and (check.coverage() == True) and (check.frames() == True) and (check.fpkm() == True):
									# Calculate SPECtre using given window length:
									window_spectre = SPECtre(region.normalized_asite_coverage(), window_length, spectre_analysis)
									if gene_type in scores:
										if chr in scores[gene_type]:
											if strand in scores[gene_type][chr]:
												if gene in scores[gene_type][chr][strand]:
													if transcript in scores[gene_type][chr][strand][gene]:
														if feature not in scores[gene_type][chr][strand][gene][transcript]:
															scores[gene_type][chr][strand][gene][transcript][feature] = {"SPEC": {"score_windowed": window_spectre.score(), "coherences_windowed": window_spectre.signal()}}
													else:
														scores[gene_type][chr][strand][gene][transcript] = {feature: {"SPEC": {"score_windowed": window_spectre.score(), "coherences_windowed": window_spectre.signal()}}}
												else:
													scores[gene_type][chr][strand][gene] = {transcript: {feature: {"SPEC": {"score_windowed": window_spectre.score(), "coherences_windowed": window_spectre.signal()}}}}
											else:
												scores[gene_type][chr][strand] = {gene: {transcript: {feature: {"SPEC": {"score_windowed": window_spectre.score(), "coherences_windowed": window_spectre.signal()}}}}}
										else:
											scores[gene_type][chr] = {strand: {gene: {transcript: {feature: {"SPEC": {"score_windowed": window_spectre.score(), "coherences_windowed": window_spectre.signal()}}}}}}
									else:
										scores[gene_type] = {chr: {strand: {gene: {transcript: {feature: {"SPEC": {"score_windowed": window_spectre.score(), "coherences_windowed": window_spectre.signal()}}}}}}}
									if "Full" in methods:
										full_spectre = SPECtre(region.normalized_asite_coverage(), 0, spectre_analysis)
										scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full"] = full_spectre.score()
										scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["coherences_full"] = full_spectre.signal()
									if "FLOSS" in methods:
										floss = FLOSS(region.asite_reads(), reference_distribution)
										scores[gene_type][chr][strand][gene][transcript][feature]["FLOSS"] = {"score": floss.score(), "distribution": floss.distribution()}
									if "ORFscore" in methods:
										orf = ORF(region.psite_frames())
										scores[gene_type][chr][strand][gene][transcript][feature]["ORF"] = {"score": orf.score(), "distribution": region.psite_frames()}
		logger.info("(in calculate_transcript_scores) Calculating transcript-level metrics [COMPLETE].")
		return scores

# CHANGE SO THAT SPEC SCORE OF 0.0 RESULTS IN POSTERIOR PROBABILITY OF 0.0 ALSO!!!
	def posterior_probability(score, coding_scores, noncoding_scores):
		if score > max(coding_scores):
			return 1.0
		else:
			prob_score_coding = len([n for n in coding_scores if n >= score])/float(len(coding_scores))
			prob_coding = len(coding_scores)/float(len(coding_scores)+len(noncoding_scores))
			prob_score = len([n for n in coding_scores+noncoding_scores if n >= score])/float(len(coding_scores)+len(noncoding_scores))
		return (prob_score_coding * prob_coding) / prob_score

	def windowed_posterior_probability(scores, coding_scores, noncoding_scores):
		windowed = list()
		for score in scores:
			windowed.append(posterior_probability(score, coding_scores, noncoding_scores))
		return windowed

	logger.info("(in calculate_transcript_metrics) Calculating transcript-level posterior probabilites [STARTED].")
	transcript_scores = calculate_transcript_scores(gtf, fpkms, asite_buffer, asite_tracker, psite_buffer, psite_tracker, window_length, reference_distribution, spectre_analysis, methods)
	translated_scores, untranslated_scores = build_translated_distributions(transcript_scores, fpkms, cutoff, "SPEC")
	for gene_type in transcript_scores:
		for chr in transcript_scores[gene_type]:
			for strand in transcript_scores[gene_type][chr]:
				for gene in transcript_scores[gene_type][chr][strand]:
					for transcript in transcript_scores[gene_type][chr][strand][gene]:
						for feature in transcript_scores[gene_type][chr][strand][gene][transcript]:
							if "SPEC" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]:
								if "score_windowed" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]:
									if "coherences_windowed" in transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]:
										transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_posterior"] = posterior_probability(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_windowed"], translated_scores, untranslated_scores)
										transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["coherences_posterior"] = windowed_posterior_probability(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["coherences_windowed"], translated_scores, untranslated_scores)
										if "Full" in methods:
											translated_scores, untranlsated_scores = build_translated_distributions(transcript_scores, fpkms, cutoff, "Full")
											transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full_posterior"] = posterior_probability(transcript_scores[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full"], translated_scores, untranslated_scores)
	logger.info("(in calculate_transcript_metrics) Calculating transcript-level posterior probabilites [COMPLETE].")
	return transcript_scores

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

	def translation_threshold(self):
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level translation threshold [STARTED].")
		translated, untranslated = build_translated_distributions(self.stats, self.fpkms, self.cutoff, "SPEC")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
		r('scores <- rbind(active, inactive)')
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level translation threshold [COMPLETE].")
		return str(r('quantile(scores$SPEC[scores$biotype=="not_translated"], probs=%s)' %(1-float(self.fdr)))[0])

	def spectre_auc(self):
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level SPECtre AUC [STARTED].")
		if "SPECtre" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = build_translated_distributions(self.stats, self.fpkms, self.cutoff, "SPEC")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level SPECtre AUC [COMPLETE].")

	def floss_auc(self):
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level FLOSS AUC [STARTED].")
		if "FLOSS" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = build_translated_distributions(self.stats, self.fpkms, self.cutoff, "FLOSS")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			return str(r('performance(prediction(-scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level FLOSS AUC [COMPLETE].")

	def orfscore_auc(self):
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level ORFscore AUC [STARTED].")
		if "ORFscore" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = build_translated_distributions(self.stats, self.fpkms, self.cutoff, "ORF")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level ORFscore AUC [COMPLETE].")

	def full_auc(self):
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level Spectral Coherence AUC [STARTED].")
		if "Full" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = build_translated_distributions(self.stats, self.fpkms, self.cutoff, "Full")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
		logger.info("EXPERIMENT_METRICS: Calculating experiment-level Spectral Coherence AUC [COMPLETE].")

def print_metrics(output_file, transcript_stats, experiment_stats, gtf, fpkms, analyses, parameters, verbose_check):

	def format_coordinates(coords):
		return ",".join([str(start) + "-" + str(end) for start, end in coords])

	def return_coordinates(gtf, gene_type, chr, strand, gene, transcript, feature):
		try:
			return format_coordinates(gtf[gene_type][chr][strand][gene][transcript][feature])
		except KeyError:
			return "NA"

	def return_metric(stats, gene_type, chr, strand, gene, transcript, feature, analysis):
		try:
			if analysis == "SPECtre":
				return str(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_windowed"]) + "\t" + str(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_posterior"])
			elif analysis == "Full":
				return str(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full"]) + "\t" + str(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["score_full_posterior"])
			elif analysis == "FLOSS":
				return str(stats[gene_type][chr][strand][gene][transcript][feature]["FLOSS"]["score"])
			elif analysis == "ORFscore":
				return str(stats[gene_type][chr][strand][gene][transcript][feature]["ORF"]["score"])
		except KeyError:
			if analysis in ("SPECtre", "Full"):
				return "NA\tNA"
			else:
				return "NA"

	def format_windowed(signal):
		return ",".join([str(sig) for sig in signal])

	def format_distribution(dist):
		return ",".join([str(length) + ":" + str(proportion) for length, proportion in sorted(dist.items())])

	def return_extras(stats, gene_type, chr, strand, gene, transcript, feature, analysis):
		try:
			if analysis == "SPECtre":
				return format_windowed(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["coherences_windowed"]) + "\t" + format_windowed(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["coherences_posterior"])
			elif analysis == "Full":
				return format_windowed(stats[gene_type][chr][strand][gene][transcript][feature]["SPEC"]["coherences_full"])
			elif analysis == "FLOSS":
				return format_distribution(stats[gene_type][chr][strand][gene][transcript][feature]["FLOSS"]["distribution"])
			elif analysis == "ORFscore":
				return format_distribution(stats[gene_type][chr][strand][gene][transcript][feature]["ORF"]["distribution"])
		except KeyError:
			if analysis == "SPECtre":
				return "NA\tNA"
			else:
				return "NA"

	def format_fpkm(transcript_fpkms, transcript):
		if transcript in transcript_fpkms:
			return transcript_fpkms[transcript]
		else:
			return "NA"

	def write_parameters(parameters, analyses):
		parameters_string = "\n# Analyses = " + str(analyses) + "\n# Window Length = " + str(parameters.len) + "\n# FPKM Cutoff = " + str(parameters.min) + "\n# FDR = " + str(parameters.fdr) + "\n# Metric = " + str(parameters.type)
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

	# Initialize the default header text:
	header = "id\tchr\tstrand\tgene_id\ttranscript_id\tgene_type\tribo_fpkm"
	if verbose_check == True:
		header += "\tcoordinates_5UTR\tcoordinates_CDS\tcoordinates_3UTR"
	for analysis in analyses:
		if analysis == "SPECtre":
			header += "\tSPEC_metric_5UTR\tSPEC_prob_5UTR\tSPEC_metric_CDS\tSPEC_prob_CDS\tSPEC_metric_3UTR\tSPEC_prob_3UTR"
			if verbose_check == True:
				header += "\tSPEC_coherences_5UTR\tSPEC_posteriors_5UTR\tSPEC_coherences_CDS\tSPEC_posteriors_CDS\tSPEC_coherences_3UTR\tSPEC_posteriors_3UTR"
		if analysis == "Full":
			header += "\tFULL_metric_5UTR\tFULL_metric_CDS\tFULL_metric_3UTR"
			if verbose_check == True:
				header += "\tFULL_coherences_5UTR\tFULL_coherences_CDS\tFULL_coherences_3UTR"
		if analysis == "FLOSS":
			header += "\tFLOSS_metric_5UTR\tFLOSS_metric_CDS\tFLOSS_metric_3UTR"
			if verbose_check == True:
				header += "\tFLOSS_dist_5UTR\tFLOSS_dist_CDS\tFLOSS_dist_3UTR"
		if analysis == "ORFscore":
			header += "\tORF_metric_5UTR\tORF_metric_CDS\tORF_metric_3UTR"
			if verbose_check == True:
				header += "\tORF_reads_5UTR\tORF_reads_CDS\tORF_reads_3UTR"

	output_file.write("# Parameters:" + write_parameters(parameters, analyses) + "\n")
	#output_file.write("\n# Experiment Metrics:" + write_experiment_metrics(experiment_stats))

	output_file.write(header)
	count = 1
	for gene_type in transcript_stats:
		for chr in transcript_stats[gene_type]:
			for strand in transcript_stats[gene_type][chr]:
				for gene in transcript_stats[gene_type][chr][strand]:
					for transcript in transcript_stats[gene_type][chr][strand][gene]:
						line = "\t".join(str(field) for field in [count, chr, strand, gene, transcript, gene_type, format_fpkm(fpkms, transcript)])
						if verbose_check == True:
							if gene_type == "protein_coding":
								for feature in ("UTR5", "CDS", "UTR3"):
									line += "\t" + return_coordinates(gtf, gene_type, chr, strand, gene, transcript, feature)
							else:
								line += "\tNA\t" + return_coordinates(gtf, gene_type, chr, strand, gene, transcript, "exon") + "\tNA"
						else:
							for analysis in analyses:
								if analysis in ("SPECtre", "Full"):
									if gene_type == "protein_coding":
										for feature in ("UTR5", "CDS", "UTR3"):
											line += "\t" + return_metric(transcript_stats, gene_type, chr, strand, gene, transcript, feature, analysis)
									else:
										line += "\tNA\tNA\t" + return_metric(transcript_stats, gene_type, chr, strand, gene, transcript, "exon", analysis) + "\tNA\tNA"
									if verbose_check == True:
										if gene_type == "protein_coding":
											for feature in ("UTR5", "CDS", "UTR3"):
												line += return_extras(transcript_stats, gene_type, chr, strand, gene, transcript, feature, analysis)
										else:
											line += return_extras(transcript_stats, gene_type, chr, strand, gene, transcript, "exon", analysis)
								else:
									if gene_type == "protein_coding":
										for feature in ("UTR5", "CDS", "UTR3"):
											line += "\t" + return_metric(transcript_stats, gene_type, chr, strand, gene, transcript, feature, analysis)
									else:
										line += "\tNA\t" + return_metric(transcript_stats, gene_type, chr, strand, gene, transcript, "exon", analysis) + "\tNA"
						output_file.write("\n" + line)
						count += 1

if __name__ == "__main__":
	
	# Parse arguments from command-line:
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, add_help=False)
	parser.add_argument("--help", action="help", help="show this help message and exit")
	parser.add_argument("--full", action="store_true", default="store_false", help="calculate un-windowed spectral coherence")
	parser.add_argument("--floss", action="store_true", default="store_false", help="calculate FLOSS and distribution")
	parser.add_argument("--orfscore", action="store_true", default="store_false", help="calculate ORFScore and distribution")
	parser.add_argument("--sanitize", action="store_true", default="store_false", help="remove overlapping transcripts")
	parser.add_argument("--verbose", action="store_true", default="store_false", help="print optional metrics including coordinates, coherence and posterior values, etc.")
	parser.add_argument("--discovery", action="store_true", default="store_false", help="enable discovery mode, see README for details")
	spectre_args = parser.add_argument_group("parameters for SPECtre analysis:")
	spectre_args.add_argument("--len", action="store", required=False, nargs="?", default=30, metavar="INT", help="length of sliding window (default: %(default)s)")
	spectre_args.add_argument("--min", action="store", required=False, nargs="?", default=5, metavar="FLOAT", help="minimum FPKM for active translation (default: %(default)s read)")
	spectre_args.add_argument("--fdr", action="store", required=False, nargs="?", default=0.05, metavar="FLOAT", help="FDR cutoff (default: %(default)s)")
	spectre_args.add_argument("--type", action="store", required=False, nargs="?", default="median", metavar="TYPE", choices=["mean","median","max","nonzero_mean","nonzero_median"], help="metric for SPECtre analysis (choices: mean,[median],max,nonzero_mean,nonzero_median)")
	file_args = parser.add_argument_group("required input and output parameters:")
	file_args.add_argument("--input", action="store", required=True, nargs="?", metavar="BAM", type=str, help="location of BAM alignment file")
	file_args.add_argument("--output", action="store", required=True, nargs="?", default="./spectre_results.txt", metavar="FILE", help="write results to (default: %(default)s)")
	file_args.add_argument("--fpkm", action="store", required=True, nargs="?", metavar="FILE", type=str, help="location of Cufflinks isoforms.fpkm_tracking file")
	file_args.add_argument("--gtf", action="store", required=True, nargs="?", metavar="FILE", type=str, help="location of GTF annotation file")
	file_args.add_argument("--log", action="store", required=True, nargs="?", default=".spectre_results.log", metavar="FILE", help="track progress to (default: %(default)s)")
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

	# Extract transcripts, transcript intervals and expression using the provided GTF:
	transcripts = Annotation(args.gtf, args.fpkm, coding + noncoding + lncrna + pseudogene)
	transcript_intervals = transcripts.interval_tree()
	transcript_fpkms = transcripts.fpkms()
	transcript_gtf = transcripts.anno()

	# Sanitize transcripts to be tested if filtering of overlapping transcripts has been enabled:
	if (args.sanitize == True):
		transcript_gtf = remove_overlapping_transcripts(transcript_gtf, transcript_intervals)

	# Default transcript boundary buffers:
	asite_buffers = {"CDS": (45,15), "UTR": (15,15)}
	psite_buffers = {"CDS": (3,3), "UTR": (3,3)}

	# Initialize the types of analyses to be conducted (default: SPECtre):
	analyses = ["SPECtre"]
	if args.full == True:
		analyses.append("Full")

	# Convert the aligned reads from the provided BAM/SAM file into an IntervalTree()
	# coverage object (may be necessary to compute coverage based on the default
	# Ingolia (A-site) offsets, and the Bazzini (P-ste) offsets. A-site adjusted
	# coverage is the default for SPECtre analysis.
	asite_reads = Alignment(args.input, "ingolia")
	psite_reads = Alignment(args.input, "bazzini")

	# Calculate the reference read distribution for FLOSS metric calculation:
	reference_reads = Reference(transcript_gtf["protein_coding"], asite_reads.coverage(), asite_buffers)

	if args.floss == True:
		analyses.append("FLOSS")
	if args.orfscore == True:
		analyses.append("ORFscore")

	# Perform a first-pass analysis of each transcript, including: SPECtre score (over each
	# designated window length), SPECtre transcript codon signal (for each window length
	# specified), FLOSS metric and read distribution (if required), and ORFscore metric and reading frame
	# distribution (if required).
	transcript_metrics = calculate_transcript_metrics(transcript_gtf, transcript_fpkms, asite_buffers, asite_reads.coverage(), psite_buffers, psite_reads.coverage(), args.min, args.len, reference_reads.distribution(), args.type, analyses)

	# Perform a second-pass global analysis based on the transcript-level metrics, such as:
	# ROC analyses, posterior probability as a function of empirical FDR, and codon window
	# signal plots (based on windowed spectral coherence).
	
	experiment_metrics = ExperimentMetrics(transcript_metrics, transcript_fpkms, analyses, args.min, args.fdr)
	print "EXPERIMENT METRICS:"
	print "Threshold =", experiment_metrics.translation_threshold()
	print "SPECtre_AUC =", experiment_metrics.spectre_auc()
	print "FULL_AUC =", experiment_metrics.full_auc()
	print "FLOSS_AUC =", experiment_metrics.floss_auc()
	print "ORFscore_AUC =", experiment_metrics.orfscore_auc()

	# Print the results table to the output file:
	print_metrics(open(args.output,"w"), transcript_metrics, experiment_metrics, transcript_gtf, transcript_fpkms, analyses, args, args.verbose)






