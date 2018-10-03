#!/usr/bin/env python

"""SPECtre

Usage: 
    spectre.py [options] <annotations> <bam>

Arguments:
	annotations	Ensembl GTF or UCSC knownGene file
	bam	Indexed BAM alignment file
	
Options:
    -h, --help	Show this screen
    -o Path to output results [default: ./spectre_results]
    -m Path to the UCSC-Ensembl map file [default: None]
    -
    -t Specify annotation format as GTF [default: GTF]
    -n NUM 	Use <n> threads of analysis [default: 1]
"""


# Import standard modules:
import sys
import re
import math
import pickle
import operator
import itertools
import collections

from functools import partial
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

# Import third-party modules:
from docopt import docopt
import HTSeq as hts
import pandas as pd
import numpy as np

# Import SPECtre modules:
from utils import *
from anno import *
#from metagene import *
#from scoring import *

def main():
	# Read in arguments/parameters:
	arguments = docopt(__doc__)
	annotation_file = arguments['<annotation_file>']
	bam_file = arguments



	print(arguments)
	annofile = arguments['<annotation_file>']
	bamfile = arguments['<bam_file>']
	outpath = arguments['-o']

	annotype = arguments['-t']
	print(outpath)
	print(annofile)
	print(bamfile)

	# Validate the chromosome formats in the supplied input files:
	validated = check_input_chromosomes(bamfile=bam_file, annofile=annotation_file, annotype=annotation_type)
	if not validated:
		# Validation failure is a fatal error:
		print('Input validation failed, please check that the choromosomes in your BAM file match your reference annotations.')
		sys.exit()
	else:	
		# Count the number of mapped reads:
		n_mapped_reads = count_mapped_reads(infile=bamfile)

		# Parse the read alignment offsets:
		read_offset_positions = parse_custom_offsets(offsets=offsets_file) if offsets_file else None

		# Load BAM file into memory as an HTSeq object:
		alignments = HTSeq.BAM_Reader(bam_file)

		# Adjust the alignments based on the read offset positions:
		offset_alignments = offset_read_alignment_positions(bam=alignments, offsets=read_offset_positions)

		# Build the transcript annotation database:
		anno_db = (load_ensembl_annotations(infile=annotation_file) if annotation_type == 'GTF'
			else load_ucsc_annotations(infile=annotation_file, mapfile=map_file) if (map_file is not None
			and annotation_type == 'UCSC') else load_ucsc_annotations(infile=annotation_file, mapfile=None)
			if (map_file is None and annotation_type == 'UCSC')) else None

		# Extract the coverage over each transcript/region, then normalize coverage to number of mapped reads:
		coverage = calculate_coverage_over_transcript_regions(db=anno_db, bam=bam_file, nreads=n_mapped_reads, threads=nt)

		# Calculate the SPECtre score over each transcript/region:
		coverage = calculate_coherence_over_regions(db=coverage)

		# Calculate the transcripts per million mapped reads (TPM) over each transcript/region:
		coverage = calculate_transcripts_per_million(db=coverage)

		# Build the active/inactive score distributions:
		model = build_translational_probability_model(db=coverage, tpm_minimum=cutoff)
		threshold = binary_search_translational_threshold(df=model)

		# Calculate the posterior probability of translation for each transcript/region:
		coverage = calculate_posterior_probability_by_region(db=coverage, model=model, cutoff=cutoff)

		# Save and/or pickle all databases and results to standard output:

if __name__ == "__main__":
	main()
