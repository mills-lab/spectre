#!/usr/bin/env python

"""Transcript scoring using coherence.

This modules takes an input DataFrame of gene- and transcript-level annotation
information, including start and end coordinates for 5'UTRs, CDS, and 3'UTRs.
The coverage over these coordinates is extracted from a user-defiend BAM file
and the SPECtre score is calculated.
"""

def calculate_coverage_over_transcript_regions(db=None, bam=None, nreads=None, threads=1):
	"""Extract read coverage over the regions in a transcript.

	"""
	def initialize_coverage_dataframe(database=None):
		"""Initialize a dataframe for the calculated coverage.

		"""
		try:
			if database is not None:
				coverage = database.transcript_id
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
		return (','.join([str(d) for d in coverage]) if row.strand == '+' else 
			','.join([str(d) for d in coverage[::-1]]))

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

def calculate_transcripts_per_million(database=None):
	"""Calculate the transcripts per million mapped reads (TPM) for each transcript/region.

	"""
	def initialize_tpm_dataframe(database=None):
		"""Initialize a dataframe for the calculate TPMs.

		"""
		try:
			if database is not None:
				tpms = database[['transcript_id', 'transcript_type', 'utr5_starts', 'cds_starts', 'utr3_starts']]
				# Add length and read count columns for 5'UTR, CDS, 3'UTR:
				tpms.utr5_length = None
				tpms.utr3_length = None
				tpms.cds_length = None
				tpms.utr5_count = None
				tpms.utr3_count = None
				tpms.cds_count = None
			else:
				raise ValueError('Missing database input')
		except ValueError:
			return None
		return tpms



#def calculate_spectral_coherence():




















