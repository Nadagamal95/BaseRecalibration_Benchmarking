# Verification of the significance of base quality score re-calibration for different datasets of human genome reads

**I. Introduction**

Phred quality score (Q score), is the most common metric used to assess the accuracy of a sequencing platform. It indicates the probability that a given base is called incorrectly by the sequencer, i.e base quality scores express how confident the machine was that it called the correct base each time.	
Variant calling algorithms rely heavily on the quality score assigned to the individual base calls in each sequence read. However, quality scores produced by the machines are subject to various sources of systematic (non-random) technical error, leading to over- or under-estimated base quality scores in the data. Some of these errors are due to the physics or the chemistry of how the sequencing reaction works, and some are probably due to manufacturing flaws in the equipment. The technique itself causes intrinsic errors such as color or laser cross-talk, cross-talk between adjacent clusters, phasing, and dimming. The problem is that systematic errors can easily be mistaken for heterozygous sites in individuals, or for SNPs in population analyses. Systematic errors are particularly problematic in low coverage experiments, or in estimates of allele-specific expression from RNA-Seq data.
Base quality score recalibration (BQSR) is a data pre-processing step to variant calling that detects systematic errors made by the sequencer when it estimates the quality score of each base call. It is a method of adjusting quality scores to be more accurate through considering every base in the data file as a whole not solely. It is a process in which machine learning models these errors empirically and adjusts the quality scores accordingly. This allows more accurate base qualities overall, which in turn improves the accuracy of variant calls.
The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants, then it adjusts the base quality scores in the data based on the model.

**II. Aim**

In our project, we planned to compare variant calling on different types of datasets and examine the influence of BQSR on the quality of variant calling. There are no clear guidelines until now to state whether BQSR is an essential step in variant calling pipeline or not. BQSR is recommended, however, on many occasions it is computationally expensive and takes time. We would like to invest our knowledge to determine a suggestive cut off for variant calling results on which on can decide if they need to apply BQSR on their sample reads or not. We have searched literature for advice on BQSR, unfortunately, we didn't find any.
		
**III. Material and methods**

*Materials*

i- Software:
•	Burrows-Wheeler Alignment Tool (BWA)
•	Sequence Alignment/Map tools (SAMtools)
•	Picard Tools
•	Genome Analysis Toolkit (GATK)
•	Tabix Tool
•	Real Time Genomics (RTG)
•	R and R libraries, ggplot2 

ii- Hardware: 
A computer with as much memory and computing power as possible. 
iii- Datasets
The applied datasets in this project were downloaded from: 
•	https://www.ncbi.nlm.nih.gov/sra/SRX4941314[accn] 
•	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP053196
•	ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz 

Each dataset was aligned to a reference genome by using BWA mem (https://github.com/molgenis/NGS_DNA); resulting in the generation of sam files which later were converted to sorted bam files using samtools.

*Methods*

**V. Discussion**

*Challenges*

1. VCF file of known variants have been very challenging. We believe that VCF files are not well-curated even if they are uploaded on reliable data resources like Ensemble for example. We struggled to solve many errors found in whole-genome VCF file (latest release) and also individual VCF files for chromosomes 1, 15, and 21.
2. To get reliable mapping quality and pipeline, a high computational power is needed (RAM and disk space).
3. Not all errors of VCF files are troubleshooted online.
4. For reliable BQSR results, sample sizes and reads should be high enough. With inadqutae disk space and slow Internet connection. It takes several hours to download a reasonable size dataset.

*Recommendations*

There is an





**VI. References**

1. https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
2. http://zenfractal.com/2014/01/25/bqsr/
3. Frazer Meacham el al. Identification and correction of systematic error in high-throughput sequence data. BMC Bioinformatics (2011), 12:451
4. Franziska Pfeiffer et al. Systematic evaluation of error rates and causes in short samples in next-generation sequencing. Scientific Reports (2018), volume 8, Article number: 10950
5. Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics (2009) 25, 2078-9. [PMID: 19505943].
6. McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. GENOME RESEARCH (2010), 20:1297-303.
7.  R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/
8. Teder, Hindrek et al. TAC-seq: targeted DNA and RNA sequencing for precise biomarker molecule counting. NPJ genomic medicine vol. 3 34. 18 Dec. 2018, doi:10.1038/s41525-018-0072-5.























 
  








 










