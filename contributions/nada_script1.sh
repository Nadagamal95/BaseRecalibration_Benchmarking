source activate ngs1

#download whole sra data

mkdir -p ~/ngs2_project/sample_data && cd ~/ngs2_project/sample_data

#neuroblastoma sample #but it is too large, so we can not work on it 
#https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8670773
#wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR867/SRR8670773/SRR8670773.sra

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

#add SRAtoolkit to the PATH:

PATH=$PATH: ~/ngs1_project/sratoolkit.2.9.6-ubuntu64/bin

####################################################

#Download human reference

mkdir -p ~/ngs2_project/sample_data && cd ~/ngs2_project/sample_data

#wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

#wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

#uncompress the fa.gz files:

#gunzip human_g1k_v37.fasta.gz

#it does not work on 

###################################

wget -c ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz

#uncompress the fa.gz files:

gunzip Homo_sapiens.GRCh38.dna.alt.fa.gz

#does not work out
#########################################################

# Download the primary assembly genome fasta File

#wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
#gunzip GRCh38.primary_assembly.genome.fa.gz

#doe not work
##########################################################

#download the human genome GRCH38p12 fasta file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.p12.genome.fa.gz

gunzip GRCh38.p12.genome.fa.gz

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

#mkdir -p ~/ngs2_project/bwa_align/bwaIndex3 && cd ~/ngs2_project/bwa_align/bwaIndex3
#ln -s ~/ngs2_project/sample_data/Homo_sapiens.GRCh38.dna.alt.fa .
#bwa index -a bwtsw Homo_sapiens.GRCh38.dna.alt.fa

#error: 
#[bwt_gen] Finished constructing BWT in 688 iterations.
#[bwa_index] 5825.29 seconds elapse.
#[bwa_index] Update BWT... [bwt_bwtupdate_core] Failed to allocate 3099750752 bytes at bwtindex.c line 158: Cannot allocate memory

##########################################################

#mkdir -p ~/ngs2_project/bwa_align/bwaIndex && cd ~/ngs2_project/bwa_align/bwaIndex
#ln -s ~/ngs2_project/sample_data/GRCh38.primary_assembly.genome.fa .
#bwa index -a bwtsw GRCh38.primary_assembly.genome.fa
#error:
#[bwt_gen] Finished constructing BWT in 688 iterations.
#[bwa_index] 5825.29 seconds elapse.
#[bwa_index] Update BWT... [bwt_bwtupdate_core] Failed to allocate 3099750752 bytes at bwtindex.c line 158: Cannot allocate memory

##########################################################

mkdir -p ~/ngs2_project/bwa_align/bwaIndex && cd ~/ngs2_project/bwa_align/bwaIndex
ln -s ~/ngs2_project/sample_data/GRCh38.p12.genome.fa .
bwa index -a bwtsw GRCh38.p12.genome.fa

[bwt_gen] Finished constructing BWT in 718 iterations.
[bwa_index] 6062.21 seconds elapse.
[bwa_index] Update BWT... [bwt_bwtupdate_core] Failed to allocate 3252208928 bytes at bwtindex.c line 158: Cannot allocate memory

##########################################################

#sequence alignment

cd ~/ngs2_project/bwa_align

R1="$HOME/ngs2_project/sample_data/SRR8617597_1.fastq"
R2="$HOME/ngs2_project/sample_data/SRR8617597_2.fastq"
/usr/bin/time -v bwa mem bwaIndex/GRCh38.primary_assembly.genome.fa $R1 $R2 > GRCh38.primary_assembly.genome.sam

##########################################################

#install samtools

source activate ngs1
conda install samtools

##########################################################

#Index your alignment file

# Convert the SAM file into a BAM file that can be sorted and indexed:
#samtools view -hbo BD143_TGACCA_L005.bam BD143_TGACCA_L005.sam

# Sort the BAM file by position in genome:
#samtools sort BD143_TGACCA_L005.bam -o BD143_TGACCA_L005.sorted.bam

# Index the BAM file so that we can randomly access it quickly:
#samtools index BD143_TGACCA_L005.sorted.bam

##########################################################

#generate & sort BAM file
#cd ~/ngs2_project/bwa_align

#for samfile in *.sam;do
#  sample=${samfile%.sam}
#  samtools view -hbo $sample.bam $samfile
#  samtools sort $sample.bam -o $sample.sorted.bam
#  samtools index $sample.sorted.bam
#done

##########################################################

#Visualize mapping using tview in samtools

#samtools tview -p chr5:62155107 BD143_TGACCA_L005.sorted.bam bwaIndex/dog_chr5.fa

##########################################################

#Call variants

#Install BCFTools

#conda install bcftools
#bcftools mpileup -Ou -f bwaIndex/dog_chr5.fa BD143_TGACCA_L005.sorted.bam |\
#bcftools call -Ov -mv > BD143_TGACCA_L005.vcf

##########################################################

