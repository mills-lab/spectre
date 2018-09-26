#!/usr/bin/env python

"""Transcript scoring using coherence.

This modules takes an input DataFrame of gene- and transcript-level annotation
information, including start and end coordinates for 5'UTRs, CDS, and 3'UTRs.
The coverage over these coordinates is extracted from a user-defiend BAM file
and the SPECtre score is calculated.
"""

# Import standard libraries:
import math
from scipy import signal

# Import SPECtre utilities:
from utils import *

def calculate_coverage_over_transcript_regions(db=None, bam=None, nreads=None, threads=1):
	"""Extract read coverage over the regions in a transcript.

	"""
	try:
		coverage_db = initialize_coverage_dataframe(database=db)
		if coverage_db is not None:
			coverage_db.utr5_raw = coverage_db.apply(lambda x: extract_coverage_over_region(x, region='5UTR'))
			coverage_db.utr3_raw = coverage_db.apply(lambda x: extract_coverage_over_region(x, region='3UTR'))
			coverage_db.cds_raw = coverage_db.apply(lambda x: extract_coverage_over_region(x, region='CDS'))
		else:
			raise ValueError('Missing transcript database input')
	except ValueError:
		return None
	# Calculate the normalized read coverage over each region in the transcript:
	coverage_db.utr5_norm = coverage_db.apply(lambda x: calculate_normalized_coverage(x, region='5UTR'))
	coverage_db.utr3_norm = coverage_db.apply(lambda x: calculate_normalized_coverage(x, region='3UTR'))
	coverage_db.cds_norm = coverage_db.apply(lambda x: calculate_normalized_coverage(x, region='CDS'))
	# Return the transcript coverage dataframe:
	return coverage_db

def calculate_transcripts_per_million(db=None):
	"""Calculate the transcripts per million mapped reads (TPM) for each transcript/region.

	"""
	try:
		if db is not None:
			# Calculate the length of each region:
			db.utr5_length = db.apply(lambda x: calculate_region_length(x, region='5UTR'))
			db.utr3_length = db.apply(lambda x: calculate_region_length(x, region='3UTR'))
			db.cds_length = db.apply(lambda x: calculate_region_length(x, region='CDS'))
			# Calculate the length of the transcripts from its constituent regions:
			db.gene_length = sum([db.utr5_length, db.cds_length, db.utr3_length])
			# Extract the read counts over each region:
			db.utr5_count = db.apply(lambda x: extract_read_counts_in_region(x, region='5UTR', coverage_db=coverage))
			db.utr3_count = db.apply(lambda x: extract_read_counts_in_region(x, region='3UTR', coverage_db=coverage))
			db.cds_count = db.apply(lambda x: extract_read_counts_in_region(x, region='CDS', coverage_db=coverage))
			# Calculate the read counts over the entire transcript:
			db.gene_count = sum([db.utr5_count, db.cds_count, db.utr3_count])
			# Calculate the TPM over each region:
			db = calculate_tpm(database=db)
        else:
        	raise ValueError('Missing database input')
    except ValueError:
    	return None
	return db

def calculate_coherence_over_region(db=None):
	"""Calculate the Welch's coherence over each transcript/region.

	"""
	try:
		if db is not None:
			db.gene_coh = db.apply(lambda x: calculate_coherence(row=x, region='gene'), axis=1)
			db.utr5_coh = db.apply(lambda x: calculate_coherence(row=x, region='5UTR'), axis=1)
			db.utr3_coh = db.apply(lambda x: calculate_coherence(row=x, region='3UTR'), axis=1)
			db.cds_coh = db.apply(lambda x: calculate_coherence(row=x, region='CDS'), axis=1)
		else:
			raise ValueError('Missing database input')
	except ValueError:
		return None
	return db

#def build_translational_probability_model():


#def calculate_posterior_probability():















