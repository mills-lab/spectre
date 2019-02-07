#!/usr/bin/env python

"""Aggregate transcript window parser

This code is part of the spectre distribution and governed by its license. Please see the LICENSE file that should have
been included as part of this package.

This module takes as input a Pandas DataFrame of parsed transcript annotation structures and builds an aggregate tran-
script window database of shared regions across a gene family (ie. transcript isoforms). The purpose of these
aggregated regions is to ensure that the SPECtre scoring distributions are built off of unique sections of
the transcriptome.

For example, where U, E, and I indicate untranslated regions, exons, and introns respectively, the shared aggregate
regions (S) of a family of transcript isoforms may be considered as:

                AUG
POSN    0123456789012345678901234567890123456789012345678901234567890123456789
ISO1    UUUUUUUUEEEEEEEEEEEIIIIIIIEEEEEEEEIIIIIIIIIIEEEEEEEEEEEEEIIIIEEEEEUUUU
ISO2    ....UUUUEEEEEEEEEEEIIIIIIIEEEEEEEEIIIIIIIIIIEEEEEIIIIIIIIIIIIEEEEEUUUU
ISO3    ..UUUUUUEEEEEEEEEEEIIIIIIIEEEEEEEEIIIIEEEEEEEEEEEEEEEEIIIIIIIEEEEEUUUU
META    ....SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS....SSSSS....
SPEC    ........SSSSSSSSSSS.......SSSSSSSS....SSSSSSSSSSSSSSSSSSS....SSSSS....

Where ISO1, ISO2 and ISO3 are annotated isoforms of the same gene, and META defines the shared aggregate regions across
all three isoforms, and SPEC indicate the coding portion of those shared isoform regions. Typically, SPECtre will only
score the shared coding regions of the aggregated metagene to build its scoring distributions.

"""

def check_annotations(gene_types=None):
	if gene_types is None:
		return False
	else:
		return True if (len(list(set(gene_types))) == 1 and 
			list(set(gene_types))[0] == 'protein_coding') else False

def get_minimal_region()




metagene_fields = ['chrom','start','end','name','status','strand','meta_starts','meta_ends']
metagene_dtypes = [np.dtype('unicode_'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('unicode_'),
	np.dtype('uint8'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_')]

mdf = initialize_annotation_dataframe(columns=metagene_fields, dtypes=metagene_dtypes, index=None)


starts = list(set(d.cds_start))
for start in starts:
	member_ids = [name for name in list(d.transcript_id[d.cds_start == start])]
	member_types = list(d.gene_type[d.cds_start == start])
	if check_annotations(gene_types=member_types) == True:
		if len(member_ids) > 1:
			print(member_ids)
			chains = list()
			sets = list()
			for member in member_ids:
				print('cds_starts:')
				print(d.cds_starts[d.transcript_id == member].values[0])
				print('cds_ends:')
				print(d.cds_ends[d.transcript_id == member].values[0])
				s = convert_coordinates(d.cds_starts[d.transcript_id == member].values[0])
				e = convert_coordinates(d.cds_ends[d.transcript_id == member].values[0])
				chain = expand_exons_to_chain(exons=list(zip(s,e)))
				sets.append(set(chain))
				chains.append(chain)
			shared = sorted(list(set.union(*sets)))
			blocks = difflib.SequenceMatcher(None, *chains)
			shared_index = blocks.get_matching_blocks()[0].size
			print('union:')
			print(shared)
			print('shared:')
			print(shared[:shared_index])




