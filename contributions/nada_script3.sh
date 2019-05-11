#somatic best practice
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146
#https://github.com/gatk-workflows/gatk4-somatic-snvs-indels
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146#
source activate ngs1

#download whole sra data

mkdir -p ~/ngs2_project/sample_data && cd ~/ngs2_project/sample_data

#lung cancer sample 
#https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8617597
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR861/SRR8617597/SRR8617597.sra

#lung cancer sample2
#https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8617625
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR861/SRR8617625/SRR8617625.sra

####################################################

#Download the SRAtoolkit


wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz

#unzip the SRAtoolkit

tar -xzf sratoolkit.2.9.6-ubuntu64.tar.gz

#####################################################

#Download human reference chromosome 

mkdir -p ~/ngs2_project/sample_data && cd ~/ngs2_project/sample_data

#chr11
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr11.fa.gz
gunzip chr11.fa.gz

#chr5

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr5.fa.gz
gunzip chr5.fa.gz

##########################################################

#Extract the First 5M Reads From basic file (PAIRED ends)and convert them from SRA file to FASTQ file

source activate ngs1

cd ~/ngs2_project/sample_data
fastq-dump --split-files -X 1000000 SRR8617597.sra 
fastq-dump --split-files -X 1000000 SRR8617625.sra 

#counting number of sequences in the fastq file
grep '@' SRR8617597.fastq | wc -l

##########################################################

#BWA Alignment

#install bwa

source activate ngs1
conda install -c bioconda bwa

##########################################################

#index your genome

#for chr5

#mkdir -p ~/ngs2_project/bwa_align/bwaIndex && cd ~/ngs2_project/bwa_align/bwaIndex
#ln -s ~/ngs2_project/sample_data/chr11.fa .
#bwa index -a bwtsw chr11.fa

#for chr5

mkdir -p ~/ngs2_project/bwa_align/bwaIndex2 && cd ~/ngs2_project/bwa_align/bwaIndex2
#ln -s ~/ngs2_project/sample_data/chr5.fa .
#bwa index -a bwtsw chr5.fa
##########################################################

#sequence alignment

cd ~/ngs2_project/bwa_align

#for chromosome 11

R1="$HOME/ngs2_project/sample_data/SRR8617597_1.fastq"
R2="$HOME/ngs2_project/sample_data/SRR8617597_2.fastq"
/usr/bin/time -v bwa mem bwaIndex/chr11.fa $R1 $R2 > SRR8617597.sam

#for chromosome 5

R1="$HOME/ngs2_project/sample_data/SRR8617597_1.fastq"
R2="$HOME/ngs2_project/sample_data/SRR8617597_2.fastq"
/usr/bin/time -v bwa mem bwaIndex/chr5.fa $R1 $R2 > SRR8617597_2.sam

##########################################################

#install samtools

source activate ngs1
conda install samtools

#Index your alignment file
#http://www.htslib.org/doc/samtools.html
#for chr11

# Convert the SAM file into a BAM file that can be sorted and indexed:

samtools view -hbo SRR8617597.bam SRR8617597.sam

# Sort the BAM file by position in genome:

samtools sort SRR8617597.bam -o SRR8617597.sorted.bam

# Index the BAM file so that we can randomly access it quickly:

samtools index SRR8617597.sorted.bam

############################################################

#Visualize mapping using tview in samtools

#samtools tview -p SRR8617597.sorted.bam bwaIndex/chr11.fa

#error
#samtools tview: cannot read index for "bwaIndex/chr11.fa"

############################################################

#Call variants

Install BCFTools
conda install bcftools

bcftools mpileup -Ou -f bwaIndex/chr11.fa SRR8617597.sorted.bam |\
bcftools call -Ov -mv > SRR8617597.vcf

############################################################

#mapping QC

samtools depth SRR8617597.sorted.bam | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > output.cov
samtools flagstat SRR8617597.sorted.bam > output.stat

############################################################

#Mark duplicate
#https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/picard_sam_markduplicates_MarkDuplicates.php

#for sample in *.sorted.bam;do
#  name=${sample%.sorted.bam}
#  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
#done

  java  -Xmx2g -jar ~/miniconda3/pkgs/picard-2.19.2-0/share/picard-2.19.2-0/picard.jar MarkDuplicates INPUT=SRR8617597.sorted.bam OUTPUT=SRR8617597.dedup.bam METRICS_FILE=SRR8617597.metrics.txt;

############################################################

#Install GATK
#conda install -c bioconda gatk4

#indexing

  java -Xmx2g -jar ~/miniconda3/pkgs/picard-2.19.2-0/share/picard-2.19.2-0/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR8617597.dedup.bam
  
##############################################################

# Reference

ln -s ~/ngs2_project/sample_data/chr11.fa .
java -Xmx2g -jar ~/miniconda3/pkgs/picard-2.19.2-0/share/picard-2.19.2-0/picard.jar CreateSequenceDictionary R=chr11.fa O=chr11.dict
samtools faidx chr11.fa

##############################################################

#Download known varinats

# Download known polymorphic sites

mkdir -p ~/ngs2_project/VCF && cd ~/ngs2_project/VCF
wget 'ftp://ftp.ensembl.org/pub/release-89/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz' -O Homo_sapiens.vcf.gz
gunzip Homo_sapiens.vcf.gz

#downloading the vcf dbsnp for chr11
#https://www.ibm.com/downloads/cas/LY1OY9XJ

https://www.ncbi.nlm.nih.gov/genome/guide/human/
ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz

#it's size is tooo large two be downloaded and it will take two much time my laptop want be able to download it.

# Select variants on chr5 and correct chr name

grep "^#" Homo_sapiens.vcf > Homo_adjusted.vcf
perl -pe 's/^([^#])/chr\1/' Homo_sapiens.vcf >> Homo_adjusted.vcf

#https://www.biostars.org/p/201603/

grep "^chr11" Homo_adjusted.vcf >> Homo_adjusted_chr11.vcf
#or
grep -w "^chr11" Homo_adjusted.vcf >> Homo_adjusted_chr11.vcf

#another way 

#grep "^#" Homo_sapiens.vcf > Homo_adjusted.vcf
grep "^11" Homo_adjusted_2.vcf | sed 's/^11/chr11/' >> Homo_adjusted_chr11.vcf
gatk IndexFeatureFile -F Homo_adjusted_chr11.vcf

sudo apt install gawk

gawk 'BEGIN{FS="\t"; OFS="\t"}{if (NF>1 && $5=="") {$5="."; print $0} else print $0}' Homo_adjusted.vcf > Homo_adjusted_2.vcf

gatk IndexFeatureFile -F Homo_adjusted_chr11.vcf

A USER ERROR has occurred: Cannot read file:///home/ngs-01/ngs2_project/VCF/Homo_adjusted_chr11_2.vcf because no suitable codecs found


# I tried alot to get a solution for this error but I failed 

#https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/org_broadinstitute_hellbender_tools_IndexFeatureFile.php
#https://github.com/broadinstitute/gatk/issues/4184
#https://www.ibm.com/downloads/cas/LY1OY9XJ

##################################################################

#Recalibrate Bases BQSR

gatk --java-options "-Xmx2G" BaseRecalibrator \
-R dog_chr5.fa -I $sample --known-sites canis_fam_chr5.vcf \
-O $name.report


