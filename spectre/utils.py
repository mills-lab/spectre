#!/usr/bin/env python

"""Utilities and helper functions for SPECtre analysis.

This module compiles the various utitlities and functions required for the SPECtre analysis
pipeline. Utitlities and functions are housed in a central repository to ensure that the
supporting modules for annotation, metagene generation, scoring, and visualization are
efficiently structured and organized. Broadly, the functions and utilities contained within
this module are structured as follows:

1) General - miscellaneous functions employed by more than one module
2) Annotation - specific to the Annotation module
3) Metagene - specific to the Metagene generation module
4) Scoring - specific to the Scoring module 
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

def extract_region_chain(row=None, region=None):
    """Calculate the length of a region using the start and end coordinates.

    """
    try:
        if all([row is not None, region is not None]):
            if region == 'gene':
                region_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.exon_starts),
                    convert_coordinates(row.exon_ends))))
            elif region == '5UTR':
                region_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.utr5_starts),
                    convert_coordinates(row.utr5_ends))))
            elif region == '3UTR':
                region_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.utr5_starts),
                    convert_coordinates(row.utr5_ends))))
            elif region == 'CDS':
                region_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.utr5_starts),
                    convert_coordinates(row.utr5_ends))))
            else:
                raise TypeError('Invalid region input')
        else:
            raise ValueError('Missing row or region input')
    except (TypeError, ValueError):
        return None
    return region_chain

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

##############
# ANNOTATION #
##############
def parse_regions_from_exons(region=None, row=None):
    """Parse coding and non-coding regions from exon coordinates.
    
    This function takes in an individual transcript in the form of an annotation DataFrame
    row, and parses the exon starts and ends into a continiguous coordinate chain, then
    partitions the coordinates into the 5'UTR, CDS, and 3'UTR based on the annotated CDS
    start and end coordinates.
    """
    coordinates = None
    if not all([region is not None, row is not None]):
        pass
    else:
        try:
            # Attempt to build the coordinate chain from the exon starts and ends:
            exon_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.exon_starts), 
                convert_coordinates(row.exon_ends))))
            if not exon_chain:
                raise ValueError('Exon conversion to coordinates chain failed')
        except ValueError:
            pass
        if region == 'utr5_starts':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(start) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos < int(row.cds_start)])]))
        if region == 'utr5_ends':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(end) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos < int(row.cds_start)])]))
        if region == 'cds_starts':
            coordinates = (','.join([str(start) for start, end in collapse_chain_to_ranges([pos for pos 
                in exon_chain])])) if row.cds_start == row.cds_end else (','.join([str(start) for start, end 
                    in collapse_chain_to_ranges([pos for pos in exon_chain if (pos >= int(row.cds_start) and pos <= 
                        int(row.cds_end))])]))
        if region == 'cds_ends':
            coordinates = (','.join([str(end) for start, end in collapse_chain_to_ranges([pos for pos 
                in exon_chain])])) if row.cds_start == row.cds_end else (','.join([str(end) for start, end 
                    in collapse_chain_to_ranges([pos for pos in exon_chain if (pos >= int(row.cds_start) and pos <= 
                        int(row.cds_end))])]))
        if region == 'utr3_starts':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(start) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos > int(row.cds_end)])]))
        if region == 'utr3_ends':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(end) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos > int(row.cds_end)])]))
    return None if not coordinates else coordinates

def add_ensembl_record(record=None, database=None, cols=None):
    """Adds an Ensembl-formatted record into the transcript annotation database.

    This function parses an Ensembl GTF annotation file and extracts transcript coordinate
    information into a DataFrame object. Transcript records are first scanned for a valid
    'transcript_id' attribute, then parsed according to the type of record. Start and stop
    codons are annotated to each transcript, and exon coordinates are appended to a transient
    'exon_starts' and 'exon_ends' column, from which CDS coordinates are inferred at a later
    point in the transcript annotation pipeline.
    """
    
    def get_transcript_type(record):
        return 'biotype' if 'transcript_biotype' in record.attr else 'type' if 'transcript_type' in record.attr else None
    
    def modify_transcript(rec=None, db=None):
        try:
            if not all([rec is not None, db is not None]):
                raise ValueError('Missing record or database input')
            else:
                if rec.type == 'start_codon':
                    # Since CDS start and end coordinates are reported according to strand:
                    if rec.iv.strand == '+':
                        db.cds_start[db.transcript_id == rec.attr['transcript_id']] = rec.iv.start
                    else:
                        db.cds_end[db.transcript_id == rec.attr['transcript_id']] = rec.iv.end
                elif rec.type == 'stop_codon':
                    if rec.iv.strand == '+':
                        db.cds_end[db.transcript_id == rec.attr['transcript_id']] = rec.iv.end
                    else:
                        db.cds_start[db.transcript_id == rec.attr['transcript_id']] = rec.iv.start
                elif rec.type == 'exon':
                    # Exon coordinates to be re-sorted later:
                    if not db.exon_starts[db.transcript_id == rec.attr['transcript_id']].values[0]:
                        db.exon_starts[db.transcript_id == rec.attr['transcript_id']] = (
                            str(rec.iv.start))
                        db.exon_ends[db.transcript_id == rec.attr['transcript_id']] = (
                            str(rec.iv.end))
                    else:
                        db.exon_starts[db.transcript_id == rec.attr['transcript_id']] += (
                            ',' + str(rec.iv.start))
                        db.exon_ends[db.transcript_id == rec.attr['transcript_id']] += (
                            ',' + str(rec.iv.end))
                else:
                    pass
        except ValueError:
            pass
        return db
    try:
        if not all([record is not None, database is not None]):
            raise ValueError('Invalid record or database input')
        else:
            if 'transcript_id' in record.attr:
                # Only records with a transcript_id attribute are to be parsed:
                if not (database.transcript_id.any() == record.attr['transcript_id']) and record.type == 'transcript':
                    database = database.append(dict(zip(cols, [
                        record.attr['gene_name'],
                        record.attr['gene_id'],
                        record.attr['transcript_id'],
                        record.attr['gene_' + get_transcript_type(record)],
                        record.attr['transcript_' + get_transcript_type(record)],
                        record.source,
                        record.iv.chrom,
                        record.iv.strand,
                        record.iv.start,
                        record.iv.end,
                        record.iv.start,
                        record.iv.end,
                        str(),
                        str(),
                        str(),
                        str(),
                        str(),
                        str(),
                        str(),
                        str()])), ignore_index=True)
                else:
                    database = modify_transcript(rec=record, db=database)
    except ValueError:
        pass
    return database

def add_ucsc_record(record=None, database=None, cols=None):
    """Adds a UCSC record into the transcript annotation database.

    This function parses an UCSC knownGene transcript annotation file into a DataFrame object
    for scoring by SPECtre. Since UCSC knownGene records are contained within a single line,
    much of the annotation parsing is comprised of transferring them into the DataFrame and
    parsing out 5'UTR, CDS, and 3'UTR coordinates.
    """

    def get_transcript_type(rec):
        return 'non_coding' if rec['cds_start'] == rec['cds_end'] else 'protein_coding'
    try:
        # UCSC transcripts are contained within a single line:
        if all([record, database is not None]):
            if not database.transcript_id.any() == record['name']:
                # Add the gene to the database:
                database = database.append(dict(zip(cols, [
                    record['protein_id'],
                    None,
                    record['name'],
                    get_transcript_type(record),
                    get_transcript_type(record),
                    'UCSC',
                    record['chrom'],
                    record['strand'],
                    record['tx_start'],
                    record['tx_end'],
                    record['cds_start'],
                    record['cds_end'],
                    record['exon_starts'][:-1], # Strips trailing comma
                    record['exon_ends'][:-1], # Strips trailing comma
                    str(),
                    str(),
                    str(),
                    str(),
                    str(),
                    str()])), ignore_index=True)
            else:
                raise ValueError('Transcript already exists in database')
        else:
            raise NameError('Invalid record or database input')
    except (NameError, ValueError):
        pass
    return database

def reorder_coordinates(region=None, row=None):
    try:
        if all([region is not None, row is not None]):
            if region == 'exon_starts':
                ordered = ','.join([str(pos) for pos in sorted([int(pos) for pos in row['exon_starts'].split(',')])])
            if region == 'exon_ends':
                ordered = ','.join([str(pos) for pos in sorted([int(pos) for pos in row['exon_ends'].split(',')])])
            if not ordered:
                raise ValueError('Coordinates could not be re-ordered based on input')
        else:
            raise ValueError('Invalid coordinates input')
    except ValueError:
        return None
    return ordered

def map_gene_name(mappings=None, row=None):
    # Maps UCSC identifiers to a designated gene identifier and symbol:
    mapped = None
    if all([mappings is not None, row is not None]):
        if row.transcript_id in mappings.ucsc.values:
            mapped = mappings.symbol[mappings.ucsc == row.transcript_id].values[0]
    return mapped
    
def map_gene_id(mappings=None, row=None):
    # Maps UCSC identifiers to a designated gene identifier and symbol:
    mapped = None
    if all([mappings is not None, row is not None]):
        if row.transcript_id in mappings.ucsc.values:
            mapped = mappings.ensembl[mappings.ucsc == row.transcript_id].values[0]
    return mapped

############
# METAGENE #
############



###########
# SCORING #
###########
def initialize_coverage_dataframe(database=None):
    """Initialize a dataframe for the calculated coverage.

    """
    try:
        if database is not None:
            coverage = database[['transcript_id', 'transcript_type', 'utr5_starts', 'cds_starts', 'utr3_starts']]
            # Add columns for 5'UTR, CDS, 3'UTR, and raw and normalized coverage:
            coverage.utr5_raw = None
            coverage.utr3_raw = None
            coverage.cds_raw = None
        else:
            raise ValueError('Missing database input')
    except ValueError:
        return None
    return coverage

def extract_coverage_over_interval(coverage=None, interval=None):
    """Extract the read coverage by position over an HTSeq GenomicInterval().
  
    """
    try:
        if all([coverage is not None, interval is not None]): 
            interval_coverage = list()
            for iv, depth in coverage[interval].steps():
                if isinstance(depth, set):
                    interval_coverage.extend([len(depth)]*len([pos for pos in iv.xrange()]))
                else:
                    interval_coverage.extend([depth]*len([pos for pos in iv.xrange()]))
        else:
            raise TypeError('Missing coverage or interval input')
    except TypeError:
        return None
    return interval_coverage

def extract_coverage_over_region(row=None, region=None):
    """Extract the coverage over a given transcript region (5'UTR, CDS, 3'UTR).

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
                raise TypeError('Invalid region input')
        else:
            raise ValueError('Missing row or region input')
    except (TypeError, ValueError):
        return None
    try:
        if ranges is not None:
            coverage = list()
            for iv_start, iv_end in ranges:
                iv = hts.GenomicInterval(chrom=row.chrom, start=iv_start, end=iv_end, strand=row.strand)
                coverage.extend(extract_coverage_over_interval(coverage=bam, interval=iv))
        else:
            raise ValueError('Missing ranges input')
    except ValueError:
        return None
    return (','.join([str(d) for d in coverage]) if row.strand == '+' else 
        ','.join([str(d) for d in coverage[::-1]]))

def normalize_coverage(coverage=None, n_mapped_reads=None):
    """Calculate the normalized read depth.

    """
    try:
        if coverage is None:
            raise ValueError('Missing coverage input')
        else:
            if n_mapped_reads is not None:
                normalized_coverage = [int(depth)/n_mapped_reads/1e6 for depth in coverage.split(',')]
            else:
                max_depth = max([int(depth) for depth in coverage])
                normalized_coverage = [int(depth)/max_depth for depth in coverage.split(',')]
    except ValueError:
        return None
    return normalized_coverage

def calculate_normalized_coverage(row=None, region=None, mapped_reads=None):
    """Calculate the normalized read coverage from the raw coverage.

    """
    try:
        if all([row is not None, region is not None]):
            if region == '5UTR':
                normalized = normalize_coverage(coverage=row.utr5_raw, n_mapped_reads=mapped_reads)
            elif region == '3UTR':
                normalized = normalize_coverage(coverage=row.utr3_raw, n_mapped_reads=mapped_reads)
            elif region == 'CDS':
                normalized = normalize_coverage(coverage=row.cds_raw, n_mapped_reads=mapped_reads)
            else:
                raise TypeError('Invalid region input')
        else:
            raise ValueError('Missing row or region input')
    except (TypeError, ValueError):
        return None
    return ','.join([str(nd) for nd in normalized])

def calculate_region_length(row=None, region=None):
    try:
        if all([row is not None, region is not None]):
            region_chain = extract_region_chain(row=row, region=region)
        else:
            raise ValueError('Missing row or region input')
    except ValueError:
        return 0
    return len(region_chain)

def extract_read_counts_in_region(row=None, region=None):
    """For a given transcript/row, extract the number of reads in each region from the coverage database.

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
                raise TypeError('Invalid region input')
        else:
            raise ValueError('Missing row or region input')
    except (TypeError, ValueError):
        return 0
    return count

def calculate_tpm(database=None):
    """For a given region, calculate the transcripts per million mapped reads.

    """
    try:
        if database is not None:
            # Divide the read count of each region by the length of the region in kilobases:
            database.gene_rpk = database.gene_count / (database.gene_length / 1000)
            database.utr5_rpk = database.utr5_count / (database.utr5_length / 1000)
            database.utr3_rpk = database.utr3_count / (database.utr3_length / 1000)
            database.cds_rpk = database.cds_count / (database.cds_length / 1000)
            # Calculate the RPK values over each transcript:
            gene_scaler = sum(database.gene_rpk) / 1e6
            utr5_scaler = sum(database.utr5_rpk) / 1e6
            utr3_scaler = sum(database.utrc_rpk) / 1e6
            cds_scaler = sum(database.cds_rpk) / 1e6
            # Divide the RPK values over each region by its respective scaler:
            database.gene_tpm = database.gene_rpk / gene_scaler
            database.utr5_tpm = database.utr5_rpk / utr5_scaler
            database.utr3_tpm = database.utr3_rpk / utr3_scaler
            database.cds_tpm = database.cds_rpk / cds_scaler
            # Clean up:
            database.drop('gene_rpk', axis=1, inplace=True)
            database.drop('utr5_rpk', axis=1, inplace=True)
            database.drop('utr3_rpk', axis=1, inplace=True)
            database.drop('cds_rpk', axis=1, inplace=True)
            del([gene_scaler, utr5_scaler, cds_scaler, utr3_scaler])
        else:
            raise ValueError('Missing database input')
    except ValueError:
        pass
    return database

def calculate_coherence(row=None, region=None):
    """Calculate the Welch's coherence over a transcript or region.

    """
    try:
        if all([row is not None, region is not None]):
            region_chain = extract_region_chain(row=row, region=region)
        else:
            raise ValueError('Missing row or region input')
    except ValueError:
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
                raise TypeError('Invalid region input')
        else:
            raise ValueError('Ideal signal could not be generated')
    except (TypeError, ValueError):
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
    """Calculate the empirical error based on the labeled translational status.

    """
    try:
        if df is not None:
            n_inactive = len(list(df.spec[df['status'] == 'inactive']))
            n_active = len(list(df.spec[df['status'] == 'active']))
        else:
            raise ValueError('Missing dataframe input')
    except ValueError:
        return None
    try:
        rate = n_inactive / (n_inactive + n_active)
        if not rate:
            raise ValueError('Cannot calculate error rate')
    except ValueError:
        return None
    return rate

def posterior_probability(row=None, model=None, score=None):
    """Calculate the posterior probability of translation for a region.

    """
    try:
        if all([row is not None, model is not None, score is not None]):
            # Extract the scores for inactive and active regions:
            inactive_scores = list(model.cds_coh[model['status'] == 'inactive'])
            active_scores = list(model.cds_coh[model['status'] == 'active'])
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
            raise ValueError('Missing row, model or score input')
    except ValueError:
        return None
    return post
