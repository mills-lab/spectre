#SPECtre

##Description
This software is designed to identify regions of active translation from ribosome profiling sequence data. This analytical pipeline scores the translational status of each annotated region (5'UTR, CDS, exon, 3'UTR) as a function of its spectral coherence over user-defined N nucleotide sliding windows to an idealized reference coding signal. It is used by calling *SPECtre.py* with the listed parameters.

##Required third-party resources:
```
R:			https://www.r-project.org/
Rpy:		http://rpy.sourceforge.net
ROCR:		https://rocr.bioinf.mpi-sb.mpg.de/
Python:		https://www.python.org/
samtools:	http://samtools.sourceforge.net/
```

##Quick Start
Download and Install:
```
git clone git@github.com:mills-lab/spectre.git
cd spectre
chmod +x SPECtre.py
```

Index Alignments:
```
samtools index <in.bam>
```

Run SPECtre with default parameters:
```
python SPECtre.py \
		--input <in.bam> \
		--output <spectre_results.txt> \
		--log <spectre_results.log> \
		--fpkm <isoforms.fpkm_tracking> \
		--gtf <ensembl.gtf>
```

##Supporting Files
Sample BAM alignment file, Cufflinks output, and Ensembl-formatted GTF are available for testing purposes in the folder *Example*.

##Output
*SPECtre* outputs transcript-level and experiment-level translational metrics in tab-delimited text format. Example output is shown in the folder *Example*.

##Usage
```
python SPECtre.py [parameters]
```

###Parameters:

####Required File Arguments:
```
	--input, alignment file in BAM format
	--output, file to output results
	--log, file to track progress
	--fpkm, location of isoforms.fpkm_tracking file from Cufflinks output
	--gtf, location of annotation file in GTF format (only Ensembl supported currently)
```

####User-defined Analytical Arguments:
```
	--len <INTEGER>, length in nucleotides of sliding window for SPECtre analysis (default: 30 nt)
	--min <FLOAT>, minimum FPKM or reads required for classification as active translation (default: 1 read)
	--fdr <FLOAT>, FDR cutoff to use for calculation of posterior probabilities (default: 0.05)
	--type <STRING>, summary statistic to use for SPECtre score (default: median)
```

####Optional Arguments:
```
	--full, enables calculation of un-windowed spectral coherence over full length of transcript
	--floss, enables calculation of FLOSS metric (Ingolia, 2014) for each transcript
	--orfscore, enables calculation of ORFscore (Bazzini, 2014) for each transcript
	--discovery, enables "Discovery Mode" see README for details (not yet supported)
	--sanitize, removes overlapping and proximally-located transcripts from analyses
	--verbose, print additional optional metrics to the output file
```
