#!/usr/bin/python
 #
######################
# STANDARD RA#A#R #IES
######################
import sys
import re
import math
import operator
import itertools
import collections

#########################
# THIRD PARTY LIBRARIES #
#########################
import HTSeq

#####################
# GENERAL UTILITIES #
#####################
def hash():
    # Adds perl-like hashing fuctionality:
    return collections.defaultdict(hash)

# REVIEW FOR DEPRECATION
def flatten(d, parent_key="", sep=":"):
    # Returns a flat key:value object from a multi-level dict():
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def add_custom_offsets(custom_offsets_file, default):
    # Adds custom read position offsets for A-/P-site adjustment:
    read_lengths, read_offsets = list(), list()
    for line in open(custom_offsets_file):
        # Check if all input values in the line are valid:
        if all(n.isdigit() for n in line.strip().split("\t")):
            read_lengths.append(int(line.strip().split("\t")[0]))
            read_offsets.append(int(line.strip().split("\t")[-1]))
    default["custom"] = dict(zip(read_lengths, read_offsets))
    return default

def check_input_chromosomes(bam, gtf, exp, targets):
    # Check that all input files have consistent chromosome nomenclature:
    def get_bam_chroms(bam_file):
        # Extracts a list of chromosomes in the BAM header:
        header = os.popen("samtools view -H " + bam_file)
        chroms = [re.findall("SN:\w{1,5}", line)[0].split(":")[-1] for line in header if "@SQ" in line]
        return sorted(chroms)
    
    def get_gtf_chroms(gtf_file):
        # Extracts a list of chromosomes in the GTF:
        lines = os.popen("cut -f 1 %s") % (gtf_file)
        chroms = list(set([line.strip() for line in lines if "#" not in line]))
        return sorted(chroms)

    # DEPRECATED
    #def get_expression_chroms(exp_file):
    #    # RSEM output for TPM measurements do not list genomic location, thus
    #    # only Cufflinks input nomenclature needs to be verified:
    #    n = [i for i, j in enumerate(open(exp_file).readlines()[0].split("\t")) if j.upper() == "LOCUS"][0]
    #    if n:
    #       chroms = list(set([line.strip().split("\t")[n].split(":")[0] for line in open(infile).readlines()[1:]]))
    #       return sorted(chroms)
    #    else:
    #       return list()

    def get_target_chroms(target_chroms):
        if target_chroms:
            return sorted([str(chrom) for chrom in targets.split(",")])
        else:
            return list()

    if get_target_chroms(targets):
        if get_expression_chroms(exp):
            return True if all(chrom in get_bam_chroms(bam) for chrom in get_gtf_chroms(gtf) for chrom in get_expression_chroms(exp) for chrom in get_target_chroms(targets)) else False
        else:
            return True if all(chrom in get_bam_chroms(bam) for chrom in get_gtf_chroms(gtf) for chrom in get_target_chroms(targets)) else False
    else:
        if get_expression_chroms(exp):
            return True if all(chrom in get_bam_chroms(bam) for chrom in get_gtf_chroms(gtf) for chrom in get_expression_chroms(exp)) else False
        else:
            return True if all(chrom in get_bam_chroms(bam) for chrom in get_gtf_chroms(gtf)) else False


def extract_coverage_over_interval(coverage, interval):
    # Extract read coverage by position over a GenomicInterval():
    interval_coverage = list()
    try:
        for i, v in coverage[interval].steps():
            if isinstance(v, set):
                interval_coverage.extend([len(v)]*len([pos for pos in i.xrange()]))
            else:
                interval_coverage.extend([v]*len([pos for pos in i.xrange()]))
    except TypeError:
        return "NA"
    return interval_coverage

def calculate_normalized_coverage(coverage, library_size):
    # Calculate the depth in reads per million mapped reads by position:
    return [cov/library_size/1e6 for cov in coverage]

############################
# TRANSCRIPTION ANNOTATION #
############################
def parse_ensembl_gtf(annotation_file, targets):
    # This function takes as input a user-supplied GTF transcript annotation file
    # and extracts the CDS and UTR intervals of transcripts, and the start and
    # end coordinates of noncoding transcripts. The coordinates are loaded into 
    # HTSeq GenomicInterval() objects for more efficient access from memory, and
    # organized in a GenomicArrayOfSets() object.
    def extract_transcript_annotations(_gtf):
        # To guard against non-standard GTF formats, protein-coding CDS regions will be
        # parsed from the GTF for the annotation of non-CDS records as 5' or 3'UTRs:
        regions = hash()
        fields = set()
        for record in _gtf:
            fields |= set([field.split("_")[-1] for field in record.attr.keys() if re.findall("_(bio)|type", field)])
            if record.type in ("CDS", "start_codon", "stop_codon"):
                if isinstance(regions[record.attr["transcript_id"]], set):
                    regions[record.attr["transcript_id"]] |= set(list(range(record.iv.start, record.iv.end)))
                else:
                    regions[record.attr["transcript_id"]] = set(list(range(record.iv.start, record.iv.end)))
        return(regions, list(fields))

    def add_interval(_regions, _array, _biotype, _transcript, _record, _start):
        if _start:
            # Input record has been annotated as a UTR, and must be partitioned to the 5' or 3' of the
            # transcript according to its position relative the provided CDS start site:
            if _record.iv.start < _start:
                _array[_biotype][_transcript]["5UTR"].append(_record.iv) if record.iv.strand == "+" \
                    else _array[_biotype][_transcript]["3UTR"].append(_record.iv)
            else:
                _array[_biotype][_transcript]["3UTR"].append(_record.iv) if record.iv.strand == "+" \
                    else _array[_biotype][_transcript]["5UTR"].append(_record.iv)
        else:
            _array[_biotype][_transcript]["CDS"].append(_record.iv)
        return(_array)

    # Load the GTF into HTSeq:
    gtf = HTSeq.GFF_Reader(annotation_file)
    # Extract the bounds of protein-coding transcripts:
    gene_regions, biotype_fields = extract_transcript_annotations(gtf)

    if (len(biotype_fields)) > 1 or "invalid" in biotype_fields:
        # There are more than one type of biotype entries, or consists of invalid entries:
        return (None, None)
    else:
        # Initialize the transcript hash and the GenomicArrayOfSets() interval container:
        intervals = (HTSeq.GenomicArrayOfSets(chroms=targets, stranded=True) if targets 
            else HTSeq.GenomicArrayOfSets(chroms="auto", stranded=True))
        annotations = hash()
        # Parse the GTF records into memory:
        for record in gtf:
            # GTF record must have a valid 'transcript_id' entry:
            if "transcript_id" in record.attr.keys():
                biotype = record.attr["transcript_" + biotype_fields[0]]
                # If seeing this transcript for the first time, initialize a dict() for the transcript associated with this record:
                if "CDS" not in annotations[biotype][record.attr["transcript_id"]]:
                    annotations[biotype][record.attr["transcript_id"]] = {"CDS": list(), "5UTR": list(), "3UTR": list()}
                if record.type in ("CDS", "stop_codon", "UTR") or (biotype != "protein_coding" and record.type == "exon"):
                    if record.type == "UTR":
                        if record.attr["transcript_id"] in gene_regions:
                            # Partition the UTR interval based on its position relative to the annotated CDS start site:
                            cds_start = (list(gene_regions[record.attr["transcript_id"]])[0] if record.iv.strand == "+"
                                else list(gene_regions[record.attr["transcript_id"]])[-1])
                            annotations = add_interval(gene_regions, annotations, biotype, record.attr["transcript_id"], record, cds_start)
                    else:
                        # Otherwise append the record to an interval list for CDS:
                        annotations = add_interval(gene_regions, annotations, biotype, record.attr["transcript_id"], record, None)
                    # Add the interval to the GenomicArrayOfSets() object"
                    intervals[record.iv] += record.attr["transcript_id"]
    # Removed unused objects from memory:
    del(gtf, gene_regions, biotype_fields)
    # Output the results:
    return(annotations, intervals)

def parse_known_gene(annotation_file, targets):
    # This function takes as input a user-supplied knownGene flat annotation file
    # and extracts the CDS and UTR intervals of transcripts, and the start and end
    # coordinates of noncoding transcripts. The coordinates are loaded into HTSeq
    # GenomicInterval() objects for efficient access from memory, and organized in
    # a GenomicArrayOfSets object.
    def extract_coding_status(_id):
        return "protein_coding" if _id else "non_coding"
    
    def convert_coordinates(_coords):
        # Converts string formatted coordinates into a list():
        return [int(pos) for pos in _coords.split(",")[:-1] if pos.isdigit()]
    
    def expand_exons_to_chain(_exons):
        # Return a list of coordinates encompassed by the start and end pairs:
        return list(itertools.chain(*[list(range(start, end+1)) for start, end in _exons]))
    
    def collapse_chain_to_ranges(_chain):
        # Convert contiguous coordinate positions into start and end pairs:
        ranges = list()
        for key, group in itertools.groupby(enumerate(_chain), lambda i: i[0]-i[1]):
            group = list(map(operator.itemgetter(1), group))
            ranges.append((group[0],group[-1])) if len(group) > 1 else range.append(group[0])
        return ranges
    
    def extract_ranges(_record, _chain, _pos):
        try:
            if _pos == "right":
                return None if _record["cdsStart"] == _record["cdsEnd"] \
                    else collapse_chain_to_ranges(_chain[_chain.index(_record["cdsEnd"]):])
            if _pos == "left":
                return None if _record["cdsStart"] == _record["cdsEnd"] \
                    else collapse_chain_to_ranges(_chain[:(_chain.index(_record["cdsStart"])+1)])
            if _pos == "mid":
                return collapse_chain_to_ranges(_chain) if _record["cdsStart"] == _record["cdsEnd"] \
                    else collapse_chain_to_ranges(_chain[_chain.index(_record["cdsStart"]):(_chain.index(_record["cdsEnd"])+1)])
        except:
            return False
    
    def add_ranges_to_array(_ranges, _record, _array):
        try:
            if _ranges:
                for start, end in _ranges:
                    _array[HTSeq.GenomicInterval(_record["chrom"], start, end, _record["strand"])] += _record["name"]
            return _array
        except:
            print("[ERROR]: Invalid input to intervals() array " + _record["name"])
    
    def add_ranges_to_hash(_ranges, _record, _region, _hash):
        try:
            if not _ranges in (False, None):
                _hash[extract_coding_status(_record["proteinID"])][_record["name"]][_region] = (
                    [HTSeq.GenomicInterval(_record["chrom"], start, end, "+") for start, end in _ranges])
            return _hash
        except:
            print("[ERROR]: Invalid input to annotation hash() " + _record["name"] + " to " + _region)
    
    # Initialize the transcript intervals array() and the annotation hash(): 
    intervals = (HTSeq.GenomicArrayOfSets(chroms=targets, stranded=True) if targets
        else HTSeq.GenomicArrayOfSets(chroms="auto", stranded=True))
    annotations = hash()

    for line in open(annotation_file):
        record = dict(zip(["name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds",
            "proteinID","alignID"], [int(field) if field.isdigit() else field for field in line.strip().split("\t")]))
        try:
            # Extract exon starts and ends from transcript record and expand to full set of coordinates:
            exon_chains = expand_exons_to_chain(list(zip(convert_coordinates(record["exonStarts"]), 
                convert_coordinates(record["exonEnds"]))))
            # Check if the transcript has coordinates for a 5' or 3'UTR and return regional ranges:
            ranges_right, ranges_left, ranges_mid = (extract_ranges(record, exon_chains, "right"), 
                extract_ranges(record, exon_chains, "left"), extract_ranges(record, exon_chains, "mid"))
            # Check that all ranges have returned an invalid value or middle range is invalid:
            if (ranges_mid == False) or all(ranges == False for ranges in [ranges_right, ranges_left, ranges_mid]):
                # Print error message:
                print("[ERROR]: Invalid range input in " + record["name"])
            else:
                # With CDS and UTR ranges extracted, add to the transcript intervals array():
                intervals = add_ranges_to_array(ranges_mid, record, add_ranges_to_array(ranges_left, record,
                    add_ranges_to_array(ranges_right, record, intervals)))
                # Add regional ranges to the transcript annotation hash():
                annotations = add_ranges_to_hash(ranges_left, record, "3UTR", annotations) if (ranges_right and record["strand"] == "+") \
                    else add_ranges_to_hash(ranges_left, record, "5UTR", annotations) if (ranges_right and record["strand"] == "-") \
                    else annotations
                annotations = add_ranges_to_hash(ranges_left, record, "5UTR", annotations) if (ranges_left and record["strand"] == "+") \
                    else add_ranges_to_hash(ranges_left, record, "3UTR", annotations) if (ranges_left and record["strand"] == "-") \
                    else annotations
                annotations = add_ranges_to_hash(ranges_mid, record, "CDS", annotations) if (ranges_mid) else annotations
        except:
            print("[ERROR] Record " + record["name"] + " in " + annotation_file + " is not formatted correctly for parsing.")
    return(annotations, intervals)

def convert_coverage(alignments, offsets, target):
    # Take as input a BAM file of read alignments, and convert the position of aligned
    # reads to their offset position in the form of an HTSeq.GenomicArray() object:
    def calculate_offset(read_length, offsets):
        return offsets[read_length] if read_length in offsets else int(read_length / 2)
    def convert_cigar_to_reference_coordinates(cigar):
        # Convert CIGAR string into coordinates based on encoded operations:
        coordinates = list()
        for op in cigar:
            if not op.type == "N":
                coordinates.extend(range(op.ref_iv.start, op.ref_iv.end))
        return(sorted(set(coordinates)))    
    def extract_offset_position(coords, offset, strand):
        return coords[(offset-1)] if strand == "+" else coords[-offset]
    # Initialize the GenomicArray():
    array = HTSeq.GenomicArray(chroms="auto", stranded=True, typecode="i", storage="step")
    if target:
        array = HTSeq.GenomicArray(chroms=target, stranded=True, typecode="i", storage="step")
    for read in alignments:
        try:
            array[HTSeq.GenomicPosition(read.iv.chrom, extract_offset_position(convert_cigar_to_reference_coordinates(read.cigar), 
                calculate_offset(len(read.read_as_aligned.seq), offsets), read.iv.strand), read.iv.strand)] += 1
        except:
            print("[ERROR] Parsing of alignment failed:", read.iv)
    return array