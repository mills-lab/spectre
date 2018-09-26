#!/usr/bin/env python

"""Transcript scoring using coherence.

This modules takes an input DataFrame of gene- and transcript-level annotation
information, including start and end coordinates for 5'UTRs, CDS, and 3'UTRs.
The coverage over these coordinates is extracted from a user-defiend BAM file
and the SPECtre score is calculated.
"""

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

def calculate_transcripts_per_million(coverage=None, db=None):
	"""Calculate the transcripts per million mapped reads (TPM) for each transcript/region.

	"""
	try:
		tpm_db = initialize_tpm_dataframe(database=db)
		if all([coverage is not None, tpm_db is not None]):
			# Calculate the length of each region:
			tpm_db.utr5_length = tpm_db.apply(lambda x: calculate_region_length(x, region='5UTR'))
			tpm_db.utr3_length = tpm_db.apply(lambda x: calculate_region_length(x, region='3UTR'))
			tpm_db.cds_length = tpm_db.apply(lambda x: calculate_region_length(x, region='CDS'))
			# Calculate the length of the transcripts from its constituent regions:
			tpm_db.gene_length = sum([tpm_db.utr5_length, tpm_db.cds_length, tpm_db.utr3_length])
			# Extract the read counts over each region:
			tpm_db.utr5_count = tpm_db.apply(lambda x: extract_read_counts_in_region(x, region='5UTR', coverage_db=coverage))
			tpm_db.utr3_count = tpm_db.apply(lambda x: extract_read_counts_in_region(x, region='3UTR', coverage_db=coverage))
			tpm_db.cds_count = tpm_db.apply(lambda x: extract_read_counts_in_region(x, region='CDS', coverage_db=coverage))
			# Calculate the read counts over the entire transcript:
			tpm_db.gene_count = sum([tpm_db.utr5_count, tpm_db.cds_count, tpm_db.utr3_count])
			# Calculate the TPM over each region:
			tpm_db = calculate_tpm(database=tpm_db)
		else:
			raise ValueError('Missing coverage or TPM database input')
	except ValueError:
		return None
	return tpm_db


#def calculate_spectral_coherence():


#def build_translational_model():


#def score_transcripts():















