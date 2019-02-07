#!/usr/bin/env python

"""Utilities and helper functions for SPECtre analysis.

    This module compiles the various utitlities and functions required for the SPECtre analysis pipeline. Utitlities
    and functions are housed in a central repository to ensure that the supporting modules for annotation, metagene
    generation, scoring, and visualization are efficiently structured and organized. Broadly, the functions and
    utilities contained within this module are structured as follows:

        1) General: miscellaneous functions employed by more than one module
        2) Annotation: specific to the Annotation module
        3) Metagene: specific to the Metagene generation module
        4) Scoring: specific to the Scoring module 

"""

# Import standard libraries:
import itertools
import operator
from multiprocessing import Pool

# Import third-party libraries:
import pandas as pd
import HTSeq as hts
import numpy as np

###########
# GENERAL #
###########
def load_alignment_file(infile=None):
    """ Load a BAM format alignment file into an HTSeq.BAM_Reader() object.

    asdf

    """
    if 
    try:
        if infile is not None:
            alignments = HTSeq.BAM_Reader(infile)
            if not alignments:
                raise ValueError('Could not parse BAM file input')
        else:
            raise NameError('Missing BAM file input')
    except (NameError, ValueError):
        return None
    return alignments

def parse_custom_offsets(offsets_file=None):
    """ Adds custom read position offsets.

    """
    try:
        if all(offsets_file is not None):
        # Custom offsets found:
            read_lengths, read_offsets = list(), list()
            for line in open(offsets_file):
                if all(n.isdigit() for n in line.strip().split('\t')):
                    # Read length is the first column:
                    read_lengths.append(int(line.strip().split('\t')[0]))
                    # Read offset is the second column:
                    read_offsets.append(int(line.strip().split('\t')[-1]))
            offsets = dict(zip(read_lengths, read_offsets))
            if not offsets:
                raise ValueError('Could not parse offsets')
        else:
            raise NameError('Missing offsets file input')
    except (NameError, ValueError):
        return None
    return offsets

def count_mapped_reads(infile=None):
    """ Count the number of aligned reads in a BAM file.

    """
    try:
        if infile is not None:
            stdout = subprocess.Popen(['samtools', 'view', '-c', '-F', '0x904', infile], 
                stdout=subprocess.PIPE)
            if stdout:
                output = stdout.communicate()
                count = int(output[0])
                if not count:
                    raise ValueError('Invalid output error from subprocess output')
            else:
                raise ValueError('Invalid subprocess')
        else:
            raise NameError('Missing BAM file input')
    except (NameError, ValueError):
        return None
    return count

def check_input_chromosomes(bam_file=None, annotation_file=None, annotation_type=None):
    """ Checks the format of chromosomes in input files.

    """
    try:
        if all([bam_file is not None, annotation_file is not None]):
            chroms_bam = ([re.findall("SN:\w{1,5}", line)[0].split(":")[-1] for line in
                os.popen("samtools view -H " + bam_file) if "@SQ" in line] if bam_file else list())
            chroms_anno = (list(set(line.strip() for line in os.popen("cut -f %s %s" % ("1" if annotation_type == "GTF"
                else "2", annotation_file)) if "#" not in line)) if annotation_file else list())
            if len(chroms_bam) >= 1 and len(chroms_anno) >= 1:
                if all(chrom in chroms_anno for chrom in chroms_bam):
                    check = True
                else:
                    raise ValueError('Mismatch in input chromosomes')
            else:
                raise ValueError('Chromosomes could not be parsed from input')
        else:
            raise NameError('Missing BAM or annotation file input')
    except (NameError, ValueError):
        return False
    return check

def convert_cigar_to_reference_coordinates(cigar=None):
    """ Extracts the reference coordinates from an HTSeq CIGAR object.

    """
    try:
        if cigar is not None:
            coordinates = list()
            for op in cigar:
                if not op.type == 'N':
                    coordinates.extend(list(range(op.ref_iv.start, op.ref_iv.end)))
            if not coordinates:
                raise ValueError('CIGAR conversion failure')
        else:
            raise NameError('Missing or invalid CIGAR string input')
    except (NameError, ValueError):
        return None
    return sorted(set(coordinates))

def offset_read_alignment_positions(bam=None, offsets=None):
    """ Adjust the reported position of reads based on the offsets.


    Calculates the offset position of the read based on the read length, if the offset
    is not defined then set the offset position to the midpoint of the read. Based on
    these offsets, the reported reference position of the read is adjusted to the
    requisite A- or P-site position of the ribosome.
    """
    try:
        if bam is not None:
            coverage = hts.GenomicArray(chroms='auto', stranded=True, typecode='i', storage='step')
            for alignment in bam:
                offset = (offsets[len(alignment.read.seq)] if len(alignment.read.seq) in offsets
                    else len(alignment.read.seq) // 2)
                offset_pos = (convert_cigar_to_reference_coordinates(alignment.cigar)[offset-1] if alignment.iv.strand == '+'
                    else convert_cigar_to_reference_coordinates(alignment.cigar)[-offset])
                coverage[HTSeq.GenomicPosition(alignment.iv.chrom, offset_pos, alignment.iv.strand)] += 1
            if len(coverage.chrom_vectors) == 0:
                raise ValueError('Alignment position offset failure')
        else:
            raise NameError('Missing BAM input')
    except (NameError, ValueError):
        return None
    return coverage

##############
# ANNOTATION #
##############


############
# METAGENE #
############



###########
# SCORING #
###########
def initialize_scoring_dataframe(database=None):
    """ Initialize a dataframe for the calculated coverage.

    """
    try:
        if database is not None:
            scoring = database[['transcript_id', 'transcript_type', 'utr5_starts', 'cds_starts', 'utr3_starts']]
            if not scoring:
                raise ValueError('Scoring dataframe initialization failure')
            else:
                # Add columns for 5'UTR, CDS, 3'UTR, and raw and normalized coverage:
                scoring.utr5_raw = None
                scoring.utr3_raw = None
                scoring.cds_raw = None
        else:
            raise NameError('Missing database input')
    except (NameError, ValueError):
        return None
    return scoring

def extract_coverage_over_interval(coverage=None, interval=None):
    """ Extract the read coverage by position over an HTSeq GenomicInterval().
  
    """
    try:
        if all([coverage is not None, interval is not None]): 
            interval_coverage = list()
            for iv, depth in coverage[interval].steps():
                if isinstance(depth, set):
                    interval_coverage.extend([len(depth)]*len([pos for pos in iv.xrange()]))
                else:
                    interval_coverage.extend([depth]*len([pos for pos in iv.xrange()]))
            if not interval_coverage:
                raise ValueError('Coverage could not be parsed from input')
        else:
            raise NameError('Missing coverage or interval input')
    except (NameError, ValueError):
        return None
    return interval_coverage

def extract_coverage_over_region(row=None, region=None):
    """ Extract the coverage over a given transcript region (5'UTR, CDS, 3'UTR).

    """
    try:
        if all([row is not None, region is not None]):
            if region == '5UTR':
                ranges = list(zip([int(n) for n in row.utr5_starts.split(',')], 
                    [int(n) for n in row.utr5_ends.split(',')]))
            elif region == '3UTR':
                ranges = list(zip([int(n) for n in row.utr3_starts.split(',')], 
                    [int(n) for n in row.utr3_ends.split(',')]))
            elif region == 'CDS':
                ranges = list(zip([int(n) for n in row.cds_starts.split(',')], 
                    [int(n) for n in row.cds_ends.split(',')]))
            else:
                raise ValueError('Invalid region input')
        else:
            raise NameError('Missing row or region input')
    except (NameError, ValueError):
        return None
    try:
        if ranges is not None:
            coverage = list()
            for iv_start, iv_end in ranges:
                iv = hts.GenomicInterval(chrom=row.chrom, start=iv_start, end=iv_end, strand=row.strand)
                coverage.extend(extract_coverage_over_interval(coverage=bam, interval=iv))
            if not coverage:
                raise ValueError('Coverage could not be parsed from input')
        else:
            raise NameError('Missing ranges input')
    except (NameError, ValueError):
        return None
    return (','.join([str(d) for d in coverage]) if row.strand == '+' else 
        ','.join([str(d) for d in coverage[::-1]]))

def normalize_coverage(coverage=None, mapped_reads=None):
    """ Calculate the normalized read depth.

    """
    try:
        if coverage is not None:
            # Calculate the max depth, used in case the number of mapped reads
            # is not defined by the user:
            max_depth = max([int(depth) for depth in coverage.split(',')])
            # Calculate the normalized coverage:
            normalized_coverage = ([int(depth)/n_mapped_reads/1e6 for depth in coverage.split(',')] if mapped_reads
                else [int(depth)/max_depth for depth in coverage.split(',')])
            if not normalized_coverage:
                raise ValueError('Normalized coverage could not be calculated from input')
        else:
            raise NameError('Missing or invalid coverage input')
    except (NameError, ValueError):
        return None
    return normalized_coverage

def calculate_normalized_coverage(row=None, region=None, mapped_reads=None):
    """ Calculate the normalized read coverage from the raw coverage.

    """
    try:
        if all([row is not None, region is not None]):
            normalized = (normalize_coverage(coverage=row.utr5_raw, mapped_reads=mapped_reads) if region == '5UTR'
                else normalize_coverage(coverage=row.utr3_raw, mapped_reads=mapped_reads) if region == '3UTR'
                else normalize_coverage(coverage=row.cds_raw, mapped_reads=mapped_reads) if region == 'CDS'
                else None)
            if not normalized:
                raise ValueError('Normalized coverage could not be calculated from input')
        else:
            raise NameError('Missing row or region input')
    except (NameError, ValueError):
        return None
    return ','.join([str(d) for d in normalized])

def calculate_region_length(row=None, region=None):
    try:
        if all([row is not None, region is not None]):
            region_chain = extract_region_chain(row=row, region=region)
            if not region_chain:
                raise ValueError('Region chain length could not be calculated from input')
        else:
            raise NameError('Missing row or region input')
    except (NameError, ValueError):
        return 0
    return len(region_chain)

def extract_read_counts_in_region(row=None, region=None):
    """ For a given transcript/row, extract the number of reads in each region from the coverage database.

    """
    try:
        if all([row is not None, region is not None]):
            if region == '5UTR':
                count = sum([int(n) for n in row.utr5_raw.split(',')])
            elif region == '3UTR':
                count = sum([int(n) for n in row.utr3_raw.split(',')])
            elif region == 'CDS':
                count = sum([int(n) for n in row.cds_raw.split(',')])
            else:
                raise ValueError('Invalid region input')
        else:
            raise NameError('Missing row or region input')
    except (NameError, ValueError):
        return 0
    return count

def calculate_tpm(database=None):
    """ For a given region, calculate the transcripts per million mapped reads.

    """
    try:
        if database is not None:
            # Divide the read count of each region by the length of the region in kilobases:
            database.gene_rpk = database.gene_count / (database.gene_length / 1000)
            database.utr5_rpk = database.utr5_count / (database.utr5_length / 1000)
            database.utr3_rpk = database.utr3_count / (database.utr3_length / 1000)
            database.cds_rpk = database.cds_count / (database.cds_length / 1000)
            # Divide the RPK values over each region by its respective scaler:
            database.gene_tpm = database.gene_rpk / (sum(database.gene_rpk) / 1e6)
            database.utr5_tpm = database.utr5_rpk / (sum(database.utr5_rpk) / 1e6)
            database.utr3_tpm = database.utr3_rpk / (sum(database.utr3_rpk) / 1e6)
            database.cds_tpm = database.cds_rpk / (sum(database.cds_rpk) / 1e6)
            # Clean up:
            database.drop('gene_rpk', axis=1, inplace=True)
            database.drop('utr5_rpk', axis=1, inplace=True)
            database.drop('utr3_rpk', axis=1, inplace=True)
            database.drop('cds_rpk', axis=1, inplace=True)
        else:
            raise ValueError('Missing database input')
    except ValueError:
        pass
    return database

def calculate_coherence(row=None, region=None):
    """ Calculate the Welch's coherence over a transcript or region.

    """
    try:
        if all([row is not None, region is not None]):
            region_chain = extract_region_chain(row=row, region=region)
            if not region_chain:
                raise NameError('Region chain could not be generated')
        else:
            raise ValueError('Missing row or region input')
    except (NameError, ValueError):
        return None
    try:
        # Instantiate the idealized signal, and truncate the signal to the length of the test region:
        ideal_signal = ([4/6,1/6,1/6]*(math.ceil(len(region_chain)/3)))[:len(region_chain)]
        if ideal_signal is not None:
            if region == 'gene':
                test_signal = [float(s) for s in row.gene_norm.split(',')]
            elif region == '5UTR':
                test_signal = [float(s) for s in row.utr5_norm.split(',')]
            elif region == '3UTR':
                test_signal = [float(s) for s in row.utr3_norm.split(',')]
            elif region == 'CDS':
                test_signal = [float(s) for s in row.cds_norm.split(',')]
            else:
                raise NameError('Invalid region input')
        else:
            raise ValueError('Ideal signal could not be generated')
    except (NameError, ValueError):
        return None
    try:
        f, Cxy = signal.coherence(test_signal, ideal_signal, nperseg=30, noverlap=27)
        if all([f is not None, Cxy is not None]):
            i = np.where(f == 1/3)
            MSC = Cxy[i]
        else:
            raise ValueError('Could not calculate coherence from test/ideal inputs')
    except ValueError:
        return None
    return MSC

def calculate_error_rate(df=None):
    """ Calculate the empirical error based on the labeled translational status.

    """
    try:
        if df is not None:
            n_inactive = len(list(df.spec[df['status'] == 'inactive']))
            n_active = len(list(df.spec[df['status'] == 'active']))
            if all([n_inactive is None, n_active is None]):
                raise NameError('Number of inactive/active genes could not be parsed')
        else:
            raise ValueError('Missing dataframe input')
    except (NameError, ValueError):
        return None
    try:
        rate = n_inactive / (n_inactive + n_active)
        if not rate:
            raise ValueError('Cannot calculate error rate')
    except ValueError:
        return None
    return rate

def posterior_probability(row=None, model=None, score=None):
    """ Calculate the posterior probability of translation for a region.

    """
    try:
        if all([row is not None, model is not None, score is not None]):
            # Extract the scores for inactive and active regions:
            inactive_scores = list(model.cds_coh[model['status'] == 'inactive'])
            active_scores = list(model.cds_coh[model['status'] == 'active'])
            if all([inactive_scores is not None, active_scores is not None]):
                # Additional tests for extreme scores:
                if score == 0:
                    post = 0.0
                elif score >= max(active_scores):
                    post = 1.0
                else:
                    # The Bayesian posterior is calculated as:
                    # P(active|score) = (P(score|active) * P(active)) / P(score)
                    p_score_active = len([s for s in active_scores if s >= score]) / len(active_scores)
                    p_active = len(active_scores) / (len(active_scores) + len(inactive_scores))
                    p_score = len([s for s in active_scores + inactive_scores if s >= score]) / (len(active_scores) + len(inactive_scores))
                    # Calculate the Bayesian posterior for the score:
                    post = (p_score_active * p_active) / p_score
            else:
                raise NameError('Inactive/active scores could not be parsed')
        else:
            raise ValueError('Missing row, model or score input')
    except (NameError, ValueError):
        return None
    return post
