#!/usr/bin/env python

"""Annotation file parser.

This code is part of the spectre distribution and governed by its
license. Please see the LICENSE file that should have been included
as part of this package.

Transcript annotation formats accepted for use as part of the spectre analytical pipeline
are limited to the following:

1) Ensembl GFF/GTF (https://ensembl.org/info/website/upload/gff.html):

	Fields must be tab-separated. Also, all but the final field in each feature line must
 	contain a value; "empty" columns should be denoted with a '.'

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

	NOTE: The 'gene_type' will be inferred from a combination of protein identity, and the presence of annotated 5'
	and 3'UTRs (to be inferred by non-matching txStart/txEnd with cdsStart/cdsEnd).

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
	8. strand - one of '+' for forward, or '-' for minus
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

The UCSC Table Browser may be used to map UCSC gene identifiers to HGNC gene symbols and Ensembl gene locus identifiers,
using the following settings for hg38 (as of 9/21/2018). Provision of these mappings is optional, but saves the end-user
post-analysis annotation of UCSC gene identifiers:

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
filte type returned: plain text

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

A sample mapping file of UCSC identifiers to gene symbol and ENSG locus identifiers may be found in the data/ folder.

"""

# Import SPECtre utilities:
from utils import *

# Turn off warnings for chained assignments:
pd.options.mode.chained_assignment = None

##################################################################
# Global variables for transcript annotation database generation #
##################################################################
# Transcript annotation fields:
annotation_fields = ['gene_name','gene_id','transcript_id','gene_type','transcript_type','source','chrom',
	'strand','start','end','cds_start','cds_end','exon_starts','exon_ends','utr5_starts','utr5_ends',
	'cds_starts','cds_ends','utr3_starts','utr3_ends']

# Transcript annotation data types:
annotation_datatypes = [np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),
	np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('uint8'),np.dtype('U1'),np.dtype('uint32'),
	np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),
	np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),
	np.dtype('uint32')]

# Transcript annotation fields for knownGene-formatted files:
known_gene_fields = ['name','chrom','strand','tx_start','tx_end','cds_start','cds_end','exon_count',
	'exon_starts','exon_ends','protein_id','align_id']

def load_ensembl_annotations(infile=None):
	"""Loads Ensembl-formatted transcript annotations to a database.

	This function executes the following: 1) instantiates a DataFrame object to store gene- and
	transcript-level annotation information, 2) parses individual transcript records and loads
	them into the database, 3) re-orders the exon start and end coordinates into a strand
	agnostic format approximating UCSC knownGene records, 4) defines the coordinates that make
	up the 5'UTR, CDS and 3'UTR of a transcript, then 5) returns only those records defined
	by the Ensembl consortium.
	"""
	
	# Initialize a DataFrame with gene- and transcript-level information:
	db = initialize_annotation_dataframe(columns=annotation_fields, dtypes=annotation_datatypes)
	try:
		if not infile:
			raise ValueError('GTF input file not specified or missing')
		else:
			gtf = hts.GFF_Reader(infile)
			for entry in gtf:
				db = add_ensembl_record(record=entry, database=db, cols=annotation_fields)
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
	db.drop('exon_starts', axis=1, inplace=True)
	db.drop('exon_ends', axis=1, inplace=True)
	# Remove the GTF from memory:
	if gtf:
		del(gtf)
	# For now, only Ensembl annotated transcripts are supported:
	db = db[db.source=='ensembl']
	return db

def load_ucsc_annotations(infile=None, mapfile=None):
	"""Loads UCSC-formatted knownGene transcript annotations to a database.

	This function executes the following: 1) instantiates a DataFrame object to store gene- and
	transcript-level annotation information, 2) parses individual transcript records and loads
	them into the database, 3) infers the coding status of each transcript, 4) defines the
	coordinates that constitute the 5'UTR, CDS, and 3'UTR of each transcript, then 5) maps
	UCSC transcript identifiers to gene symbol and Ensembl locus identifiers.

	To get gene-level groupings, like isoforms of the same gene, Ensembl locus identifiers and
	known gene symbols may be optionally mapped to the 'gene_id' and 'gene_name' columns,
	respectively. This is based on annotated UCSC identifiers mapped to Ensembl identifiers
	and known gene symbols using the UCSC Table Browser as described in the modules top-level
	documentation. If no map file is designated, then the default 'gene_name' (proteinID in the
	knownGene formatted table structure) will be mapped to the 'gene_id' column. See the data/
	folder for sample hg38-Ensembl mappings.
	"""
	
	# Initialize a DataFrame with gene- and transcript-level information"
	try:
		db = initialize_annotation_dataframe(columns=annotation_fields, dtypes=annotation_datatypes)
		if not infile:
			raise ValueError('knownGene input file not specified or missing')
		else:
			with open(infile) as f:
				for line in f:
					entry = dict(zip(known_gene_fields, line.strip().split('\t')))
					db = add_ucsc_record(record=entry, database=db, cols=annotation_fields)
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
	if mapfile:
		try:
			ucsc_mappings = pd.read_csv(filepath_or_buffer=mapfile, sep='\t', header=None, skiprows=1, names=['ucsc','symbol','ensembl'],
				dtype={'ucsc': np.dtype('unicode_'), 'symbol': np.dtype('unicode_'), 'ensembl': np.dtype('unicode_')})
			if ucsc_mappings is None:
				raise FileNotFoundError('Missing or invalid map file input')
		except FileNotFoundError:
			pass
		db.gene_name = db.apply(lambda x: map_gene_name(mappings=ucsc_mappings, row=x), axis=1)
		db.gene_id = db.apply(lambda x: map_gene_id(mappings=ucsc_mappings, row=x), axis=1)
	# Exon start and end coordinates are no longer needed:
	db.drop('exon_starts', axis=1, inplace=True)
	db.drop('exon_ends', axis=1, inplace=True)
	# Remove UCSC-Ensembl mappings from memory:
	if ucsc_mappings is not None:
		del(ucsc_mappings)
	return db

# Used for testing purposes:
#map_file = '/mnt/c/Users/stonyc/Documents/repos/spectre/spectre/data/hg38toEnsembl.txt'
#ensembl_gtf_file = '/home/stonyc/references/ensembl/v78/annotation/test.gtf'
#known_genes_file = '/home/stonyc/references/hg38/annotation/knownTest.txt'
