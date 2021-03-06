#!/usr/bin/env python

""" Annotation file parser.

This code is part of the spectre distribution and governed by its license. Please see the LICENSE file that should 
have been included as part of this package.

Transcript annotation formats accepted for use as part of the spectre analytical pipeline are limited to the following:

1) Ensembl GFF/GTF (https://ensembl.org/info/website/upload/gff.html):

    Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" 
    columns should be denoted with a '.'

    1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
    2. source - name of the program that generated this feature, or the data source (database or project name)
    3. feature - feature type name, e.g. Gene, Variation, Similarity
    4. start - start position of the feature, with sequence numbering starting at 1.
    5. end - end position of the feature, with sequence numbering starting at 1.
    6. score - a floating point value.
    7. strand - defined as + (forward) or - (reverse).
    8. frame - one of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, 
               '1' that the second base is the first base of a codon, and so on.
    9. attribute - a semicolon-separated list of tag-value pairs, providing additional information about each feature.

    NOTE: The 'gene_type' or 'gene_biotype' attribute must be provided for spectre to build its distribution of coding
    and non-coding regions, failure to provide this attribute will result in program termination.

2) UCSC knownGene (refer to UCSC Table Browser schema for additional documentation):

    Fields must be tab-delimited.

    1. name - name of gene
    2. chrom - reference sequence chromosome or scaffold
    3. strand - + or - for strand
    4. txStart - transcript start position (or end position for minus strand item)
    5. txEnd - transcript end position (or start position for minus strand item)
    6. cdsStart - coding region start (or end position if minus strand item)
    7. cdsEnd - coding region end (or start position if minus strand item)
    8. exonCount - number of exons
    9. exonStarts - exon start positions (or end positions for minus strand item)
    10. exonEnds - exon end positions (or start positions for minus strand item)
    11. proteinID - UniProt display ID, UniProt accession, or RefSeq protein ID
    12. alignID - unique identifier (for example, GENCODE transcript ID)

    NOTE: The 'gene_type' will be inferred from a combination of protein identity, and the presence of annotated
    5' and 3'UTRs (to be inferred by non-matching txStart/txEnd with cdsStart/cdsEnd).

Gene and transcript-level information will be parsed into a Pandas DataFrame object, which can then pickled for
downstream use by other portions of the spectre package. The information parsed and tabulated will be structured
as follows:

    1. gene_name - name of gene (eg. TP53, ACTB, etc.)
    2. gene_id - gene identifier (eg. ENSG ID, RefSeq ID)
    3. transcript_id - transcript identifier (eg. ENST ID, UCSC ID)
    4. gene_type - gene type annotation
    5. transcript_type - transcript type annotation
    6. source - data source or name of program that generated this item
    7. chrom - reference sequence chromosome or scaffold
    8. strand - one of '+' for forward, or '-' for reverse
    9. start - equivalent to UCSC txStart
    10. end - equivalent to UCSC txEnd
    11. cds_start - equivalent to UCSC cdsStart
    12. cds_end - equivalent to UCSC cdsEnd
    13. exon_starts - equivalent to UCSC exonStarts
    14. exon_ends - equivalent to UCSC exonEnds
    15. utr5_starts - 5'UTR start positions (or end positions for minus strand item)
    16. utr5_ends - 5'UTR end positions (or start positions for minus strand item)
    17. cds_starts - coding region start positions (or end positions for minus strand item)
    18. cds_ends - coding region end position (or start positions for minus strand item)
    19. utr3_starts - 3'UTR start positions (or end positions for minus strand item)
    20. utr3_ends - 3'UTR end positions (or start positions for minus strand item)

This provides spectre a common data structure for UCSC and Ensembl/GENCODE alignments and workflows.

The UCSC Table Browser may be used to map UCSC gene identifiers to HGNC gene symbols and Ensembl gene locus id-
entifiers, using the following settings for hg38 (as of 9/21/2018). Provision of these mappings is optional,
but saves the end-user post-analysis annotation of UCSC gene identifiers:

    1) Navigate to genome.ucsc.edu/cgi-bin/hgTables
    2) Select the following options in the drop-down menus:

        clade: Mammal
        genome: Human
        assembly: Dec. 2013 (GRCh38/hg38)
        group: Genes and Gene Predictions
        track: GENCODE v24
        table: knownGene
        region: genome
        identifiers (names/accessions): N/A
        filter: N/A
        intersection: N/A
        correlation: N/A
        output format: selected fields from primary and related tables
        output file: <filename to output results>
        filter type returned: plain text

    3) Click on 'get_output'
    4) On the next page, select the following Linked Tables:

        hg38/knownAttrs (Fields in Gencode attrs table that aren't in kgXref)

    5) At the bottom of the page, click on 'allow selection from checked tables'
    6) Then select the following options on the following page:

        From hg38.knownGene:
        name/Name of gene

        From hg38.kgXref fields:
        geneSymbol/Gene Symbol

        From hg38.knownAttrs fields:
        geneId/ENSG* locus identifier

    7) Click on 'get_output', to save your results.

    A sample mapping file of UCSC identifiers to gene symbol and ENSG locus identifiers may be found in the data/
    folder.

"""

# Import SPECtre utilities:
from utils import *

# Turn off warnings for chained assignments:
pd.options.mode.chained_assignment = None

##################################################################
# Global variables for transcript annotation database generation #
##################################################################
# Transcript annotation fields:
annotation_fields = ['gene_name','gene_id','transcript_id','gene_type','transcript_type','source','chrom','strand',
    'start','end','cds_start','cds_end','exon_starts','exon_ends','utr5_starts','utr5_ends','cds_starts','cds_ends',
    'utr3_starts','utr3_ends']

# Transcript annotation data types:
annotation_datatypes = [np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),
    np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('uint8'),np.dtype('U1'),np.dtype('uint32'),np.dtype('uint32'),
    np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),
    np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32')]

# Transcript annotation fields for knownGene-formatted files:
known_gene_fields = ['name','chrom','strand','tx_start','tx_end','cds_start','cds_end','exon_count','exon_starts',
    'exon_ends','protein_id','align_id']

def parse_regions_from_exons(region=None, row=None):
    """ Parse coding and non-coding regions from exon coordinates.
    
    This function takes in an individual transcript in the form of an annotation DataFrame row, and parses the exon
    starts and ends into a continiguous coordinate chain, then partitions the coordinates into the 5'UTR, CDS, and
    3'UTR based on the annotated CDS start and end coordinates.

    """
    assert region in ('utr5_starts','utr5_ends','cds_starts','cds_ends','utr3_starts','utr3_ends'), (
        'Region name is required')
    assert row is not None, 'No annotation row found'
    # Expand exon coorindates to chain:
    exon_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.exon_starts), 
        convert_coordinates(row.exon_ends))))
    assert exon_chain is not None, 'Expansion of exon ranges to chain failure'
    try:
        if region == 'utr5_starts':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(start) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos < int(row.cds_start)])]))
        if region == 'utr5_ends':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(end) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos < int(row.cds_start)])]))
        if region == 'cds_starts':
            coordinates = (','.join([str(start) for start, end in collapse_chain_to_ranges([pos for pos in 
                exon_chain])])) if row.cds_start == row.cds_end else (','.join([str(start) for start, end 
                in collapse_chain_to_ranges([pos for pos in exon_chain if (pos >= int(row.cds_start) 
                and pos <= int(row.cds_end))])]))
        if region == 'cds_ends':
            coordinates = (','.join([str(end) for start, end in collapse_chain_to_ranges([pos for pos in 
                exon_chain])])) if row.cds_start == row.cds_end else (','.join([str(end) for start, end 
                in collapse_chain_to_ranges([pos for pos in exon_chain if (pos >= int(row.cds_start) 
                and pos <= int(row.cds_end))])]))
        if region == 'utr3_starts':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(start) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos > int(row.cds_end)])]))
        if region == 'utr3_ends':
            coordinates = None if row.cds_start == row.cds_end else (','.join([str(end) for start, end in 
                collapse_chain_to_ranges([pos for pos in exon_chain if pos > int(row.cds_end)])]))
        if not region:
            raise ValueError('Parsing of regions from exon ranges failure')
    except ValueError:
        return None
    return coordinates

def add_ensembl_record(record=None, database=None, cols=None):
    """ Adds an Ensembl-formatted record into the transcript annotation database.

    This function parses an Ensembl GTF annotation file and extracts transcript coordinate information into a DataFrame
    object. Transcript records are first scanned for a valid 'transcript_id' attribute, then parsed according to the
    type of record. Start and stop codons are annotated to each transcript, and exon coordinates are appended to a
    transient 'exon_starts' and 'exon_ends' column, from which CDS coordinates are inferred at a later point in
    the transcript annotation pipeline.
    """
    
    def get_transcript_type(record=None):
        assert record is not None, 'Missing Ensembl annotation record'
        try:
            annotation = ('biotype' if 'transcript_biotype' in record.attr else 'type' 
                if 'transcript_type' in record.attr else None)
            if not annotation:
                raise ValueError('Failure to parse bio|type from attributes')
        except ValueError:
            return None
        return annotation
    
    def modify_transcript(rec=None, db=None):
        assert rec is not None, 'Transcript record not found'
        assert db is not None, 'Input database not found'
        try:
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
                    db.exon_starts[db.transcript_id == rec.attr['transcript_id']] = (str(rec.iv.start))
                    db.exon_ends[db.transcript_id == rec.attr['transcript_id']] = (str(rec.iv.end))
                else:
                    db.exon_starts[db.transcript_id == rec.attr['transcript_id']] += (','+str(rec.iv.start))
                    db.exon_ends[db.transcript_id == rec.attr['transcript_id']] += (','+str(rec.iv.end))
            else:
                # Does not need further modification:
                pass
        except:
            # Enforce more rigorous 
            pass
        return db

    try:
        if all([record is not None, database is not None]):
            if 'transcript_id' in record.attr:
                # Only records with a transcript_id attribute are to be parsed:
                if not (database.transcript_id.any() == record.attr['transcript_id']) and record.type == 'transcript':
                    database = database.append(dict(zip(cols, [record.attr['gene_name'],record.attr['gene_id'],
                        record.attr['transcript_id'],record.attr['gene_' + get_transcript_type(record)],
                        record.attr['transcript_' + get_transcript_type(record)],record.source,record.iv.chrom,
                        record.iv.strand,record.iv.start,record.iv.end,record.iv.start,record.iv.end,
                        str(),str(),str(),str(),str(),str(),str(),str()])), ignore_index=True)
                else:
                    database = modify_transcript(rec=record, db=database)
            else:
                raise NameError('Invalid record or database input')
    except NameError:
        pass
    return database

def add_ucsc_record(record=None, database=None, cols=None):
    """ Adds a UCSC record into the transcript annotation database.

    This function parses an UCSC knownGene transcript annotation file into a DataFrame object
    for scoring by SPECtre. Since UCSC knownGene records are contained within a single line,
    much of the annotation parsing is comprised of transferring them into the DataFrame and
    parsing out 5'UTR, CDS, and 3'UTR coordinates.
    """

    def get_transcript_type(record=None):
        try:
            if record is not None:
                annotation = ('non_coding' if rec['cds_start'] == rec['cds_end'] 
                    else 'protein_coding')
                if not annotation:
                    raise ValueError('Failure to parse bio|type from record')
            else:
                raise NameError('Missing or invalid record input')
        except (NameError, ValueError):
            return None
        return annotation

    try:
        # UCSC transcripts are contained within a single line:
        if all([record, database is not None]):
            if not database.transcript_id.any() == record['name']:
                # Add the gene to the database:
                database = database.append(dict(zip(cols, [record['protein_id'],None,record['name'],
                    get_transcript_type(record),get_transcript_type(record),'UCSC',record['chrom'],record['strand'],
                    record['tx_start'],record['tx_end'],record['cds_start'],record['cds_end'],
                    record['exon_starts'][:-1], # Strips trailing comma
                    record['exon_ends'][:-1], # Strips trailing comma
                    str(),str(),str(),str(),str(),str()])), ignore_index=True)
            else:
                # Since UCSC transcripts are contained in one line, redundant
                # transcripts should not exist:
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
            raise NameError('Invalid coordinates input')
    except (NameError, ValueError):
        return None
    return ordered

def map_gene_name(mappings=None, row=None):
    # Maps UCSC identifiers to a designated gene identifier and symbol:
    try:
        if all([mappings is not None, row is not None]):
            if row.transcript_id in mappings.ucsc.values:
                mapped = mappings.symbol[mappings.ucsc == row.transcript_id].values[0]
                if not mapped:
                    raise ValueError('Mapping of gene name from UCSC identifier failed')
            else:
                raise ValueError('Transcript identifier missing from mappings')
        else:
            raise NameError('Invalid mappings or missing row input')
    except (NameError, ValueError):
        return None
    return mapped
    
def map_gene_id(mappings=None, row=None):
    # Maps UCSC identifiers to a designated gene identifier and symbol:
    try:
        if all([mappings is not None, row is not None]):
            if row.transcript_id in mappings.ucsc.values:
                mapped = mappings.ensembl[mappings.ucsc == row.transcript_id].values[0]
                if not mapped:
                    raise ValueError('Mapping of gene locus from UCSC identifier failed')
            else:
                raise ValueError('Transcript identifier missing from mappings')
        else:
            raise NameError('Invalid mappings or missing row input')
    except (NameError, ValueError):
        return None
    return mapped
def initialize_annotation_dataframe(columns=None, dtypes=None, index=None):
    """ Instantiates an empty dataframe.

    Creates a Pandas DataFrame object with specified columns and data types. Please refer to:
    https://bit.ly/2nloFkG

    Args:
        columns (:py:object:list): List of expected data fields for the annotation database
        dtypes (:py:object:list): List of expected data types for the annotation databse
        index (:py:object:str): Column to index the Pandas DataFrame

    Returns:
        df (:pandas:object:DataFrame): Pre-formatted empty Pandas DataFrame to store annotation records

    """
    assert columns is not None, 'Missing data headers to Pandas DataFrame'
    assert dtypes is not None, 'Missing data types to Pandas DataFrame'
    try:
        df = pd.DataFrame(index=index)
        for col, datatype in zip(columns, dtypes):
            df[col] = pd.Series(dtype=datatype)
        if df is None:
            raise ValueError('Failure to initiate annotation dataframe')
    except ValueError:
        return None
    return df

def convert_coordinates(coords=None):
    """ Convert positional coordinate string into a list of integers.

    Takes a comma-delimited string of positional coordinates, and converts them into a list() of integers.

    Args:
        coords (:py:object:str): Comma-delimited string of characters to be coerced to int()

    Returns:
        converted (:py:object:list): List of integer coordinates

    """
    assert coords is not None, 'Missing coordinates for conversion'
    try:
        converted = [int(pos) for pos in coords.split(',') if pos.isdigit()]
        if not converted:
            raise ValueError('Non-numerical input for coordinate conversion')
    except ValueError:
        return None
    return converted

def collapse_chain_to_ranges(chain=None):
    """ Convert a list of integer coordinates to a de-limited of start, end pairs.

    From a list of positional coordinates, convert to a de-limited list of tuples() defined as consecutive
    (start, end) coordinates.


    """
    assert chain is not None, 'Missing chain for collapse to exons'
    try:
        ranges = list()
        for key, group in itertools.groupby(enumerate(chain), lambda i: i[0]-i[1]):
            group = list(map(operator.itemgetter(1), group))
            ranges.append((group[0], group[-1])) if len(group) > 1 else ranges.append(group[0])
        if not ranges:
            raise ValueError('Chain collapse error')
    except ValueError:
        return None
    return ranges

def expand_exons_to_chain(exons=None):
    """ Convert a list of (start, end) coordinates, into a chain.

    From a de-limited list of paired (start, end) tuples(), expand the (start, end) pairs into a list() of consecutive
    integer coordinates.

    """
    assert exons is not None, 'Missing exons for chain expansion'
    try:
        chain = list(itertools.chain(*[list(range(start, end+1)) for start, end in exons]))
        if not chain:
            raise ValueError('Exon expansion error')
    except ValueError:
        return None
    return chain

def extract_region_chain(row=None, region=None):
    """ Calculate the length of a region using the start and end coordinates.

    """
    assert row is not None, 'Missing annotation record input'
    assert region in ('gene','5UTR','3UTR','CDS'), 'Invalid annotaton region input'
    try:
        if region == 'gene':
            region_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.exon_starts),
                convert_coordinates(row.exon_ends))))
        if region in ('5UTR','CDS','3UTR'):
            region_chain = expand_exons_to_chain(list(zip(convert_coordinates(row.utr5_starts),
                convert_coordinates(row.utr5_ends))))
        if not region_chain:
            raise ValueError('Chain extraction failure')
    except ValueError:
        return None
    return region_chain

def load_ensembl_annotations(infile=None):
    """Loads Ensembl-formatted transcript annotations to a Pandas DataFrame.

    This function executes the following: 1) instantiates a DataFrame object to store gene- and transcript-level anno-
    tation information, 2) parses individual transcript records and loads them into the database, 3) re-orders the
    exon start and end coordinates into a strand agnostic format approximating UCSC knownGene records, 4) defines
    the coordinates that make up the 5'UTR, CDS and 3'UTR of a transcript, then 5) returns only those records
    defined by the Ensembl consortium.

    Args:
        infile (:py:object:str): Path to an Ensembl-formated GTF gene annotation file.

    Returns:
       db (:pandas:object:DataFrame): Pandas DataFrame of parsed GTF annotation records.

    """
    assert infile is not None, 'Ensembl GTF not provided'
    # Initialize a DataFrame with gene- and transcript-level information:
    db = initialize_annotation_dataframe(columns=annotation_fields, dtypes=annotation_datatypes)
    try:
        gtf = hts.GFF_Reader(infile)
        for entry in gtf:
            db = add_ensembl_record(record=entry, database=db, cols=annotation_fields)
        if not db:
            raise ValueError('Ensembl record parsing failed')
    except ValueError:
        return None
    # Re-order the exon starts and ends:
    db.exon_starts = db.apply(lambda x: reorder_coordinates(region='exon_starts', row=x), axis=1)
    db.exon_ends = db.apply(lambda x: reorder_coordinates(region='exon_ends', row=x), axis=1)
    # Parse out the 5'UTR, CDS and 3'UTR coordinates:
    db.utr5_starts = db.apply(lambda x: parse_regions_from_exons(region='utr5_starts', row=x), axis=1)
    db.utr3_starts = db.apply(lambda x: parse_regions_from_exons(region='utr3_starts', row=x), axis=1)
    db.cds_starts = db.apply(lambda x: parse_regions_from_exons(region='cds_starts', row=x), axis=1)
    db.utr5_ends = db.apply(lambda x: parse_regions_from_exons(region='utr5_ends', row=x), axis=1)
    db.utr3_ends = db.apply(lambda x: parse_regions_from_exons(region='utr3_ends', row=x), axis=1)
    db.cds_ends = db.apply(lambda x: parse_regions_from_exons(region='cds_ends', row=x), axis=1)
    # Exon start and end coordinates are no longer needed:
    # db.drop('exon_starts', axis=1, inplace=True)
    # db.drop('exon_ends', axis=1, inplace=True)
    # Remove the GTF from memory:
    if gtf:
        del(gtf)
    # For now, only Ensembl annotated transcripts are supported:
    try:
        db = db[db.source=='ensembl']
    except ValueError:
        return None
    return db

def load_ucsc_annotations(infile=None, mapfile=None):
    """Loads UCSC-formatted knownGene transcript annotations to a database.

    This function executes the following: 1) instantiates a DataFrame object to store gene- and transcript-level anno-
    tation information, 2) parses individual transcript records and loads them into the database, 3) infers the coding
    status of each transcript, 4) defines the coordinates that constitute the 5'UTR, CDS, and 3'UTR of each transcript,
    then 5) maps UCSC transcript identifiers to gene symbol and Ensembl locus identifiers.

    To get gene-level groupings, like isoforms of the same gene, Ensembl locus identifiers and known gene symbols may
    be optionally mapped to the 'gene_id' and 'gene_name' columns, respectively. This is based on annotated UCSC id-
    entifiers mapped to Ensembl identifiers and known gene symbols using the UCSC Table Browser as described in the
    modules top-level documentation. If no map file is designated, then the default 'gene_name' (proteinID in the
    knownGene formatted table structure) will be mapped to the 'gene_id' column. See the data/ folder for sample
    hg38-Ensembl mappings.

    """
    assert infile is not None, 'UCSC knownGene file not provided'
    assert mapfile is not None, 'hg38-Ensembl mapping file not provided'
    # Initialize a DataFrame with gene- and transcript-level information:
    db = initialize_annotation_dataframe(columns=annotation_fields, dtypes=annotation_datatypes)
    try:
        with open(infile) as f:
            for line in f:
                entry = dict(zip(known_gene_fields, line.strip().split('\t')))
                db = add_ucsc_record(record=entry, database=db, cols=annotation_fields)
        if db is not None:
            raise ValueError('Failure to parse knownGene records')
    except ValueError:
        return None
    # Parse out the 5'UTR, CDS, and 3'UTR coordinates:
    db.utr5_starts = db.apply(lambda x: parse_regions_from_exons(region='utr5_starts', row=x), axis=1)
    db.utr3_starts = db.apply(lambda x: parse_regions_from_exons(region='utr3_starts', row=x), axis=1)
    db.cds_starts = db.apply(lambda x: parse_regions_from_exons(region='cds_starts', row=x), axis=1)
    db.utr5_ends = db.apply(lambda x: parse_regions_from_exons(region='utr5_ends', row=x), axis=1)
    db.utr3_ends = db.apply(lambda x: parse_regions_from_exons(region='utr3_ends', row=x), axis=1)
    db.cds_ends = db.apply(lambda x: parse_regions_from_exons(region='cds_ends', row=x), axis=1)
    # Map the UCSC transcript identifiers to known gene symbols and Ensembl identifiers:
    db.gene_id = db.gene_name
    try:
        ucsc_mappings = pd.read_csv(filepath_or_buffer=mapfile, sep='\t', header=None, skiprows=1, names=['ucsc',
            'symbol','ensembl'], dtype={'ucsc': np.dtype('unicode_'), 'symbol': np.dtype('unicode_'), 
            'ensembl': np.dtype('unicode_')})
        if not ucsc_mappings:
            raise FileNotFoundError('Missing or invalid map file input')
    except FileNotFoundError:
        return None
    db.gene_name = db.apply(lambda x: map_gene_name(mappings=ucsc_mappings, row=x), axis=1)
    db.gene_id = db.apply(lambda x: map_gene_id(mappings=ucsc_mappings, row=x), axis=1)
    # Exon start and end coordinates are no longer needed:
    # db.drop('exon_starts', axis=1, inplace=True)
    # db.drop('exon_ends', axis=1, inplace=True)
    # Remove UCSC-Ensembl mappings from memory:
    if ucsc_mappings is not None:
        del(ucsc_mappings)
    return db

# Used for testing purposes:
#map_file = '/mnt/c/Users/stonyc/Documents/repos/spectre/spectre/data/hg38toEnsembl.txt'
#ensembl_gtf_file = '/home/stonyc/references/ensembl/v78/annotation/test.gtf'
#known_genes_file = '/home/stonyc/references/hg38/annotation/knownTest.txt'
