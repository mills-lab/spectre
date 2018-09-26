#!/usr/bin/env python

"""Utilities and helper functions for SPECtre analysis.

"""

# Import standard libraries:
import itertools
import operator

# Import third-party libraries:
import pandas as pd
import HTSeq as hts
import numpy as np

def initialize_annotation_dataframe(columns=None, dtypes=None, index=None):
    """Instantiates an empty dataframe.

    Creates a Pandas DataFrame object with specified columns and data types. Please refer to:
    https://bit.ly/2nloFkG
    """
    try:
        if all([columns, dtypes]):
            df = pd.DataFrame(index=index)
            for column, datatype in zip(columns, dtypes):
                df[column] = pd.Series(dtype=datatype)
        else:
            raise NameError("Invalid inputs for database initialization")
    except NameError:
        return None
    return df

def convert_coordinates(coords=None):
    """Convert coordinate string into list of integers.

    Takes a comma-delimited string of positional coordinates, and converts them tino a list()
    of integers.
    """
    try:
        if coords == None:
            raise ValueError('Invalid coordinate input')
        else:
            return [int(pos) for pos in coords.split(',') if pos.isdigit()]
    except ValueError:
        return None

def collapse_chain_to_ranges(chain=None):
    """Convert a list of integer coordinates to a de-limited of start, end pairs.

    From a list of positional coordinates, convert to a de-limited list of tuples() defined as
    consecutive (start, end) coordinates.
    """
    try:
        if chain == None:
            raise ValueError('Invalid coordinate chain input')
        else:
            ranges = list()
            for key, group in itertools.groupby(enumerate(chain), lambda i: i[0]-i[1]):
                group = list(map(operator.itemgetter(1), group))
                ranges.append((group[0], group[-1])) if len(group) > 1 else ranges.append(group[0])
            return ranges
    except ValueError():
        return None

def expand_exons_to_chain(exons=None):
    """Convert a list of (start, end) coordinates, into a chain.

    From a de-limited list of paired (start, end) tuples(), expand the (start, end) pairs into a
    list() of consecutive integer coordinates.
    """
    try:
        if exons == None:
            raise ValueError('Missing exons or invalid input')
        else:
            return list(itertools.chain(*[list(range(start, end+1)) for start, end in exons]))
    except:
        return None

def parse_custom_offsets(offsets_file=None, default=None):
    """Adds custom read position offsets.

    """
    if all(offsets_file is not None, default is not None):
        try:
            # Custom offsets found:
            read_lengths, read_offsets = list(), list()
            for line in open(offsets_file):
                if all(n.isdigit() for n in line.strip().split('\t')):
                    # Read length is the first column:
                    read_lengths.append(int(line.strip().split('\t')[0]))
                    # Read offset is the second column:
                    read_offsets.append(int(line.strip().split('\t')[-1]))
            default = dict(zip(read_lengths, read_offsets))
        except:
            pass
    return default

def check_input_chromosomes(bam_file=None, annotation_file=None, annotation_type=None):
    """Checks the format of chromosomes in input files.

    """
    chroms_bam = ([re.findall("SN:\w{1,5}", line)[0].split(":")[-1] for line in
        os.popen("samtools view -H " + bam_file) if "@SQ" in line] if bam_file else list())
    chroms_anno = (list(set(line.strip() for line in os.popen("cut -f %s %s" % ("1" if annotation_type == "GTF"
        else "2", annotation_file)) if "#" not in line)) if annotation_file else list())
    # Check that the chromosomes from the BAM and annotation inputs were properly extracted:
    assert len(bam_chroms) > 1, "Error in parsing chromosomes from your BAM input: %s" % (bam_file)
    assert len(anno_chroms) > 1, "Error in parsing chromosomes from your annotation file: %s" % (annotation_file)
    # Check that the the chromosomes from the BAM, annotation and targeted regions are consistent:
    assert all(chrom in chroms_anno for chrom in chroms_bam), "Mismatch in input chromosomes."

def convert_cigar_to_reference_coordinates(cigar=None):
    """Extracts the reference coordinates from an HTSeq CIGAR object.

    """
    try:
        if not cigar:
            raise ValueError('Missing CIGAR string input')
        else:
            coordinates = list()
            for op in cigar:
                if not op.type == 'N':
                    coordinates.extend(list(range(op.ref_iv.start, op.ref_iv.end)))
    except ValueError:
        return None
    return sorted(set(coordinates))

def offset_read_alignment_positions(bam=None, offsets=None):
    """Adjust the reported position of reads based on the offsets.


    Calculates the offset position of the read based on the read length, if the offset
    is not defined then set the offset position to the midpoint of the read. Based on
    these offsets, the reported reference position of the read is adjusted to the
    requisite A- or P-site position of the ribosome.
    """
    try:
        if bam is not None:
            coverage = HTSeq.GenomicArray(chroms='auto', stranded=True, typecode='i', storage='step')
            for alignment in bam:
                offset = (offsets[len(alignment.read.seq)] if len(alignment.read.seq) in offsets
                    else len(alignment.read.seq) // 2)
                offset_pos = (convert_cigar_to_reference_coordinates(alignment.cigar)[offset-1] if alignment.iv.strand == '+'
                    else convert_cigar_to_reference_coordinates(alignment.cigar)[-offset])
                coverage[HTSeq.GenomicPosition(alignment.iv.chrom, offset_pos, alignment.iv.strand)] += 1
        else:
            raise ValueError('Missing BAM input')
    except ValueError:
        return None
    return coverage


