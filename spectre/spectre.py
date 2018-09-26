#!/usr/bin/env python

"""SPECtre

Usage: 
    spectre.py [options] <annotation_file> <bam_file>

Arguments:
	annotation_file	Ensembl GTF or UCSC knownGene file
	bam_file	Indexed BAM alignment file
	
Options:
    -h, --help	Show this screen
    -o PATH 	Output path [default: .]
    -n NUM 	Use <n> threads of analysis [default: 1]
    -t 	Specify annotation format as GTF [default: GTF]
    -m 	Path to the UCSC-Ensembl map file [default: .]
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
	print(arguments)
	annofile = arguments['<annotation_file>']
	bamfile = arguments['<bam_file>']
	outpath = arguments['-o']

	annotype = arguments['-t']
	print(outpath)
	print(annofile)
	print(bamfile)

	# Validate the chromosome formats in the input files:
	validated = check_input_chromosomes(bamfile=bamfile, annotationfile=annofile, annotationtype=annotype)
	if not validated:
		# Write out an error message:
		#sys.exit()
		print('Not validated')
	else:
		# Parse the read alignment offsets:
		read_offset_positions = parse_custom_offsets(offsets_file=offset_file, default=None)
		# Load BAM file into memory as an HTSeq object:
		alignments = HTSeq.BAM_Reader(bam_file)
		# Adjust the coverage to the A-/P-site position:
		offset_alignments = offset_read_alignment_positions(bam=alignments, offsets=read_offset_positions, targets=target_regions)
		#if not offset_alignments:
			# Write out an error message:
			#sys.exit()
		#else:
			# Build the transcript annotation database from the input file:
			#if annotype == 'UCSC':
			#	anno_db = (load_ucsc_annotations(infile=annofile, mapfile=mapfile) if mapfile 
			#		else load_ucsc_annotations(infile=annofile, mapfile=None))
			#elif annotype == 'GTF':
			#	anno_db = load_ensembl_annotations(infile=annofile)
			#else:
			#	sys.exit()

			# Extract the coverage over each transcript in the annotation database:

			# Calculate the normalized coverage over each transcript/region:

			# Calculate the transcripts per million mapped reads (TPM) over each transcript and region:

			# Calculate the SPECtre score over each transcript/region:

			# Build the coding/non-coding score distributions:

			# Calculate the posterior probability of translation for each transcript/region:

			# Write out summary statistics, and print results to output files:







if __name__ == "__main__":
	main()
