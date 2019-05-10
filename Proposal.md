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

III. Materials
``` java
i- Software:
- Bowtie 2: Fast and sensitive read alignment
- Sequence Alignment/Map tools (SAMtools)
- Picard Tools
- Genome Analysis Toolkit (GATK)
- Tabix Tool
- Real Time Genomics (RTG)
- R and R libraries, ggplot2

ii- Hardware: 
Combination of individual computers with variable RAM and disk memory

iii- Datasets
We first considered applying base quality score recalibration on a model organism other than human. We chose E.coli. We found "PathSeq" in GATK, however, it is still a beta version and we will not be able to rely on it. Haplotype caller can be used to identify the ploidy of an unknown sample, however, this is not our target here. In order for us to apply base quality score recalibration, we have to create our own set of known variants for E.coli for which we don't have quite the time. In addition, GATK best practice for germline short variant discovery mentioned that variant recalibration's algorithm requires high-quality sets of known variants to use as training and truth resources, which for many organisms are not yet available. It also requires quite a lot of data in order to learn the profiles of good vs. bad variants, so it can be difficult or even impossible to use on small datasets that involve only one or a few samples, on targeted sequencing data, on RNAseq, and on non-model organisms. GATK best practices, recommendation and testing are based on human variations. Therefore, we choose to check the significance of applying BQSR vs not applying it on human reads for DNA sequences mapped to human genome.

We decided to choose diverse datasets and GATK pipelines. We chose germline short variant discovery (SNPs + Indels) and somatic short variant discovery (SNVs + Indels) pipelines. Datasets were cell line (47,XX, +21) genomic DNA, whole-exome germline DNA, and whole-genome somatic DNA. In materials and methods, the pipeline of the first dataset is expained in detail. 
The applied data set (47,XX,+21) genomic DNA was downloaded from https://www.ncbi.nlm.nih.gov/sra/SRX4941314[accn] . It was aligned to a reference of whole genome one time and to a reference of chromosome 21 another time using bowtie (http://bowtie-bio.sourceforge.net/index.shtml); resulting in the generation of sam files which later were converted to sorted bam files using samtools.
```
IV. Methods:
Merge replicates:
The Picard tools were used to merge the replicates by using the following command:
java  -Xmx2g -jar $picard_path/picard.jar MergeSamFiles I=BD143_TGACCA_L005.sorted.bam I=BD143_TGACCA_L006.sorted.bam OUTPUT=BD143_TGACCA_merged.sorted.bam

The PicardTools are built with java, according to that when running a jar file (e.g., java -jar picard.jar <PicardTool>) a memory limit to java can be added, for example requiring java to use no more than 2GB memory: -Xmx2g. This can help ensure our program does not use more memory than we request (https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html).

 Mapping QC
Using the samt tools, quality control for the mapped reads was applied by running the following command:
for bamFile in *.sorted.bam;do
  output=${bamFile%.sorted.bam}
  samtools depth $bamFile | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > $output.cov
  samtools flagstat $bamFile > $output.stat
done
Mark duplicate
Same DNA fragments may be sequenced many times during the sequencing process. The resulted duplicate reads may cause the propagation of errors across all the subsequent duplicate reads. The duplicate reads can violate the assumptions of variant calling (https://github.com/yuhuanq/GATK-Pipeline), so we marked the duplicate reads and deleted it using the Picard tools. The following for loop was used to accomplish this step: 
for sample in *.sorted.bam;do
  name=${sample%.sorted.bam}
  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done
Indexing 
i.  Indexing dedupped bam samples file
The BuildBamIndex tool generates a BAM index ".bai" file for the input BAM, allowing a fast look-up of data in a BAM file (https://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex). In this project we indexed the dedupped bam samples file as follow: 
for sample in *.dedup.bam;do
  #name=${sample%.dedup.bam}
  java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$sample
done
ii. Indexing reference for GATK
GATK 4 tools require that the main FASTA file be accompanied by a dictionary file (.dict) and an index file (.fai), as it allows efficient random access to the reference bases. The GATK look for these index files based on their name, so it's important that they have the same basename as the FASTA file (https://software.broadinstitute.org/gatk/documentation/article?id=11013). We generated these files by applying the following command lines: 
ln -s ~/workdir/sample_data/dog_chr5.fa .
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=dog_chr5.fa O=dog_chr5.dict
samtools faidx dog_chr5.fa
Download known variants
The known variants were downloaded for the used reference from the following link: 
wget ftp://ftp.ensembl.org/pub/release-89/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz
And the variants on chromosome 22 were selected using the following command lines: 
gunzip canis_familiaris.vcf.gz
grep "^#" canis_familiaris.vcf > canis_fam_chr5.vcf
grep "^5" canis_familiaris.vcf | sed 's/^5/chr5/' >> canis_fam_chr5.vcf
gatk IndexFeatureFile -F canis_fam_chr5.vcf
Bases Recalibration  (BQSR)    
The goal of the Bases Recalibration (BQSR) procedure is to correct the systematic bias that might affect the assignment of base quality scores by the sequencer (https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr). The procedure of BQSR consists two main passes:
the first pass consists of calculating error empirically and finding patterns in how error varies with basecall features over all bases. This step is performed by using the BaseRecalibrator tool and the relevant observations are written to a recalibration table (https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php) . 
The second pass is performed by using the ApplyBQSR tool and consists of applying numerical corrections to each individual basecall based on the patterns identified in the first step (recorded in the recalibration report) and write out the recalibrated data to a new BAM file (https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php). 
The following for loop was applied to accomplish the BQSR procedure:
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R dog_chr5.fa -I $sample --known-sites canis_fam_chr5.vcf \
-O $name.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R dog_chr5.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done

	The variant calling steps were performed on both un-recalibrated and recalibrated samples.





Variant Calling
HaplotypeCaller was used to call variants using BAM files of tis sample. The output files were GVCF files which has raw, unfiltered SNP and indel calls for all sites, variant or invariant(https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php). The command  lines used to perform this step are: 
~/path/to/gatk HaplotypeCaller \
--java-options "-Xmx4g -XX:ParallelGCThreads=1" \ 
-R reference.fasta \
-I sample_1.dedup.sorted.bam \
-O sample_1.raw.g.vcf \
--emit-ref-confidence GVCF
VCF statistics:
i- index the VCF file:
The generated vcf files were indexed using Tabix. It indexes position sorted files in TAB-delimited formats and create an index file (.gz.tbi ) (Heng Li, 2011, https://doi.org/10.1093/bioinformatics/btq671). 
The indexing was performed through these command lines:
conda install -c bioconda tabix
bgzip -c raw_variants_ann.vcf > raw_variants_ann.vcf.gz
tabix -p vcf raw_variants_ann.vcf.gz
ii- Calculate some statistics:
The performance of some quick statistics were done by using the Real Time Genomic (RTG) tools rtg vcfstats command. These tools includes useful utilities for dealing with vcf files and sequence data. The most  interesting is the vcfeval command that performed comparison of vcf files (https://github.com/RealTimeGenomics/rtg-tools). 
conda install -c bioconda rtg-tools
rtg vcfstats raw_variants_ann.vcf.gz > stats.txt
In the pipeline we intended to split SNPs and Indels,assess different filters and plot its figures using R-script, in addition to SNP and Indel Variant filtrations. However we stopped at the statistics step as the stat.txt file contains zero reading for all parameters except the Passed Filters.
Split SNPs and Indels:
A vcf file containing variants needs to be subsetted in order to facilitate certain analyses (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain requirements). SelectVariants GATK tool can be used to split the variants into SNPs and Indels (https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php).
The splitting procedure was done via these command lines:
i- For SNPs:
gatk --java-options "-Xmx2G" SelectVariants \
-R dog_chr5.fa \
-V raw_variants_ann.vcf \
--select-type-to-include SNP \
-O raw_variants_ann_SNP.vcf
 
ii- For Indels:
gatk --java-options "-Xmx2G" SelectVariants \
-R dog_chr5.fa \
-V raw_variants_ann.vcf \
--select-type-to-include INDEL \
-O raw_variants_ann_INDEL.vcf

Assess different filters in both known and novel
For the filtration step, nine filters were used (https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering):
•	AN: Total number of alleles in called genotypes 
•	DP: the unfiltered depth of coverage across all samples
•	
•	QD : QualByDepth
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
 For performing the assessment of different filters, the following for loop will be used:
for var in "SNP" "INDEL";do
 input="raw_variants_ann_"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done; done
Plotting figures 
All plots are generated using the ggplot2 library in R. On the x-axis are the annotation values, and on the y-axis are the density values. The area under the density plot gives the probability of observing the annotation values (https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations). 
wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
sudo Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done
SNP and Indel Variant filteration
SNPs or Indels matching the recommended parameters will be considered bad and filtered out, i.e. marked with a filter name (which will be specified in the filtering command) in the output VCF file. SNPs or Indels that do not match any of these parameters will be considered good and marked PASS in the output VCF file (https://software.broadinstitute.org/gatk/documentation/article?id=2806).
i- SNP Variant filteration
cd ~/workdir/GATK_tutorial
gatk --java-options "-Xmx2G" VariantFiltration \
-R dog_chr5.fa \
-V raw_variants_ann_SNP.vcf \
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
ii- INDEL Variant filteration
gatk --java-options "-Xmx2G" VariantFiltration \
-R dog_chr5.fa \
-V raw_variants_ann_SNP.vcf \
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






*Challenges:
1. VCF file of known variants have been very challenging. We believe that VCF files are not well-curated even if they are uploaded on reliable data resources like Ensemble for example. We struggled to solve many errors found in whole-genome VCF file (latest release) and also individual VCF files for chromosomes 1, 15, and 21.
2. To get reliable mapping quality and pipeline, a high computational power is needed (RAM and disk space).
3. Not all errors of VCF files are troubleshooted online.
4. For reliable BQSR results, sample sizes and reads should be high enough. With inadqutae disk space and slow Internet connection. It takes several hours to download a reasonable size dataset.

*Recommendations:
There is an





*References:
1. https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
2. http://zenfractal.com/2014/01/25/bqsr/
3. Frazer Meacham el al. Identification and correction of systematic error in high-throughput sequence data. BMC Bioinformatics (2011), 12:451
4. Franziska Pfeiffer et al. Systematic evaluation of error rates and causes in short samples in next-generation sequencing. Scientific Reports (2018), volume 8, Article number: 10950
5. Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics (2009) 25, 2078-9. [PMID: 19505943].
6. McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. GENOME RESEARCH (2010), 20:1297-303.
7.  R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/
8. Teder, Hindrek et al. TAC-seq: targeted DNA and RNA sequencing for precise biomarker molecule counting. NPJ genomic medicine vol. 3 34. 18 Dec. 2018, doi:10.1038/s41525-018-0072-5.























 
  








 









