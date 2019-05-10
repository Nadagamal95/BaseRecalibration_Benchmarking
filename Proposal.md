**Verification of the significance of base quality score recalibration for different datasets of human genome reads** 

I. Introduction
``` java
    Phred quality score (Q score), is the most common metric used to assess the accuracy of a sequencing platform. It indicates the probability that a given base is called incorrectly by the sequencer, i.e base quality scores express how confident the machine was that it called the correct base each time.

Variant calling algorithms rely heavily on the quality score assigned to the individual base calls in each sequence read. However, quality scores produced by the machines are subject to various sources of systematic (non-random) technical error, leading to over- or under-estimated base quality scores in the data. Some of these errors are due to the physics or the chemistry of how the sequencing reaction works, and some are probably due to manufacturing flaws in the equipment. The technique itself causes intrinsic errors such as colour or laser cross-talk, cross-talk between adjacent clusters, phasing, and dimming. The problem is that systematic errors can easily be mistaken for heterozygous sites in individuals, or for SNPs in population analyses. Systematic errors are particularly problematic in low coverage experiments, or in estimates of allele-specific expression from RNA-Seq data.

Base quality score recalibration (BQSR) is a data pre-processing step to variant calling that detects systematic errors made by the sequencer when it estimates the quality score of each base call. It is a method of adjusting quality scores to be more accurate through considering every base in the data file as a whole not solely. It is a process in which machine learning models these errors empirically and adjusts the quality scores accordingly. This allows more accurate base qualities overall, which in turn improves the accuracy of variant calls.
The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants, then it adjusts the base quality scores in the data based on the model.
```

II. Aim
``` java
    In our project, we planned to compare variant calling on different types of datasets and examine the influence of BQSR on the quality of variant calling. There are no clear guidelines until now to state whether BQSR is an essential step in variant calling pipeline or not. BQSR is recommended, however, on many occasions it is computationally expensive and takes time. We would like to invest our knowledge to determine a suggestive cut off for variant calling results on which on can decide if they need to apply BQSR on their sample reads or not. We have searched literature for advice on BQSR, unfortunately, we didn't find any.
```

III. Material and methods

# Materials:
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

# Methods:
The following steps were applied on each of the three datasets in parallel 
Merge replicates:
Merging replicates was done by using the Picard tools were used to merge the replicates.
The PicardTools are built with java, according to that when running a jar file (e.g., java -jar picard.jar <PicardTool>) a memory limit to java can be added, for example requiring java to use no more than 2GB memory: -Xmx2g. This can help ensure our program does not use more memory than we request (https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html).

 Mapping QC
Samtools were used for the quality control of the mapped reads.
Mark duplicate
Same DNA fragments may be sequenced many times during the sequencing process. The resulted duplicate reads may cause the propagation of errors across all the subsequent duplicate reads. The duplicate reads can violate the assumptions of variant calling (https://github.com/yuhuanq/GATK-Pipeline), so we marked the duplicate reads and deleted it using the Picard tools. 
Indexing 
i.  Indexing dedupped bam samples file
The BuildBamIndex tool generates a BAM index ".bai" file for the input BAM, allowing a fast look-up of data in a BAM file (https://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex). In this project, BuildBamIndex tool was used to  index the dedupped bam samples file. 
ii. Indexing reference for GATK
GATK 4 tools require that the main FASTA file be accompanied by a dictionary file (.dict) and an index file (.fai), as it allows efficient random access to the reference bases. The GATK look for these index files based on their name, so it's important that they have the same basename as the FASTA file (https://software.broadinstitute.org/gatk/documentation/article?id=11013). We generated these files by applying the following command lines: 
Download known variants
The known variants were downloaded for the used reference from the following link: 
ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz 
Bases Recalibration  (BQSR)    
The goal of the Bases Recalibration (BQSR) procedure is to correct the systematic bias that might affect the assignment of base quality scores by the sequencer (https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr). The procedure of BQSR consists two main passes:
the first pass consists of calculating error empirically and finding patterns in how error varies with basecall features over all bases. This step is performed by using the BaseRecalibrator tool and the relevant observations are written to a recalibration table (https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php) . 
The second pass is performed by using the ApplyBQSR tool and consists of applying numerical corrections to each individual basecall based on the patterns identified in the first step (recorded in the recalibration report) and write out the recalibrated data to a new BAM file (https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php). 

######################################################################## The variant calling steps were performed on both un-recalibrated and                  #
 #   recalibrated samples.                                                                                                  #
######################################################################

# Variant Calling
HaplotypeCaller was used to call variants using BAM files of tis sample. The output files were GVCF files which has raw, unfiltered SNP and indel calls for all sites, variant or invariant (https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php). 
# VCF statistics:
i- index the VCF file:
The generated vcf files were indexed using Tabix. It indexes position sorted files in TAB-delimited formats and creates an index file (.gz.tbi ) (Heng Li, 2011, https://doi.org/10.1093/bioinformatics/btq671). 
ii- Calculate some statistics:
The performance of some quick statistics was done by using the Real Time Genomic (RTG) tools rtg vcfstats command. These tools include useful utilities for dealing with vcf files and sequence data. The most  interesting is the vcfeval command that performed comparison of vcf files (https://github.com/RealTimeGenomics/rtg-tools). 
######################################################################     #In the pipeline we intended to split SNPs and Indels, assess different filters and  #                          #plot its figures using R-script, in addition to SNP and Indel Variant filtrations.   #
 # However, we stopped at the statistics step as the stat.txt file contains zero reading 
# for all parameters except the Passed Filters.                                                            #
#######################################################################
Split SNPs and Indels:
A vcf file containing variants needs to be subsetted in order to facilitate certain analyses (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain requirements). SelectVariants GATK tool can be used to split the variants into SNPs and Indels (https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php).
The splitting procedure was meant to be done via these command lines:
i- For SNPs:
gatk --java-options "-Xmx2G" SelectVariants \
-R Homo_sapiens.GRCh38.dna.chromosome.21.fa\
-V SRR8115017_bqsr_ann.vcf\
--select-type-to-include SNP \
-O SRR8115017_bqsr_ann _SNP.vcf
 
ii- For Indels:
gatk --java-options "-Xmx2G" SelectVariants \
-R Homo_sapiens.GRCh38.dna.chromosome.21.fa \
-V SRR8115017_bqsr_ann.vcf \
--select-type-to-include INDEL \
-O SRR8115017_bqsr_ann _INDEL.vcf
# Assess different filters in both known and novel
For the filtration step, nine filters were used (https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering):
•	AN: Total number of alleles in called genotypes 
•	DP: The unfiltered depth of coverage across all samples
•	QD: QualByDepth
Variant quality score divided by the depth of the alternate allele. Recommendation: SNPs: 2.0, INDELS: 2.0
•	FS: FisherStrand
Phred-scaled p-values using Fisher's Exact Test for strand bias. Higher values are more likely false positive calls. Recommendation: SNPs: 60.0, INDELS: 200.0
•	MQ: RMSMappingQuality
Root Mean Square of the mapping quality of reads across samples. Recommendation: SNPs: 40.0
•	MQRankSum: MappingQualityRankSumTest
U-based z-approximation from Mann-Whitney Rank Sum Test for mapping qualities, comparing reads with reference bases and those with alternate alleles. Recommendation: SNPs: -12.5
•	ReadPosRankSum: ReadPosRankSumTest
U-based z-approximation from Mann-Whitney Rank Sum Test for distance from end of reads for those reads with an alternate allele. As bases near the ends of reads are more likely to contain errors, if all reads with the allele are near the end of the reads this may be indicative of an error. Recommendation: SNPs: -8.0, INDELS: -20.0
•	SOR: StrandOddsRatio
High values indicate strand bias in the data Recommendation: SNPs: 3.0, INDELS: 10.
 For performing the assessment of different filters, the following for loop will be used that will produce 18 files for SNPs and 18 files for INDELs:
for var in "SNP" "INDEL";do
 input=" SRR8115017_bqsr_ann _"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done; done
Plotting figures 
All plots are generated using the ggplot2 library in R. On the x-axis are the annotation values, and on the y-axis are the density values. The area under the density plot gives the probability of observing the annotation values (https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations). 
SNP and Indel Variant filteration
SNPs or Indels matching the recommended parameters will be considered bad and filtered out, i.e. marked with a filter name (which will be specified in the filtering command) in the output VCF file. SNPs or Indels that do not match any of these parameters will be considered good and marked PASS in the output VCF file (https://software.broadinstitute.org/gatk/documentation/article?id=2806).
i- SNP Variant filteration
gatk --java-options "-Xmx2G" VariantFiltration \
-R Homo_sapiens.GRCh38.dna.chromosome.21.fa \
-V SRR8115017_bqsr_ann _SNP.vcf \
--filter-name "snpQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "snpMQ" \
--filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filter-name "snpMQRankSum" \
--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filter-name "snpFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
--filter-name "snpSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filter-name "snpReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filter-name "snpDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_SNP_clean.vcf
ii- INDEL Variant filtration
gatk --java-options "-Xmx2G" VariantFiltration \
-R Homo_sapiens.GRCh38.dna.chromosome.21.fa \
-V SRR8115017_bqsr_ann _INDEL.vcf \
--filter-name "indelQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "indelFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 200.0" \
--filter-name "indelSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 10.0" \
--filter-name "indelReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
--filter-name "indelInbreedingCoeff" \
--filter-expression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
--filter-name "indelDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_INDEL_clean.vcf


V. Discussion
``` java



*Challenges*
1. VCF file of known variants have been very challenging. We believe that VCF files are not well-curated even if they are uploaded on reliable data resources like Ensemble for example. We struggled to solve many errors found in whole-genome VCF file (latest release) and also individual VCF files for chromosomes 1, 15, and 21.
2. To get reliable mapping quality and pipeline, a high computational power is needed (RAM and disk space).
3. Not all errors of VCF files are troubleshooted online.
4. For reliable BQSR results, sample sizes and reads should be high enough. With inadqutae disk space and slow Internet connection. It takes several hours to download a reasonable size dataset.

*Recommendations*

There is an



```

VI. References
``` java
1. https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
2. http://zenfractal.com/2014/01/25/bqsr/
3. Frazer Meacham el al. Identification and correction of systematic error in high-throughput sequence data. BMC Bioinformatics (2011), 12:451
4. Franziska Pfeiffer et al. Systematic evaluation of error rates and causes in short samples in next-generation sequencing. Scientific Reports (2018), volume 8, Article number: 10950
5. Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics (2009) 25, 2078-9. [PMID: 19505943].
6. McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. GENOME RESEARCH (2010), 20:1297-303.
7.  R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/
8. Teder, Hindrek et al. TAC-seq: targeted DNA and RNA sequencing for precise biomarker molecule counting. NPJ genomic medicine vol. 3 34. 18 Dec. 2018, doi:10.1038/s41525-018-0072-5.
```






















 
  








 










