##BWA alignment
#dowloading data
mkdir -p ~/ngs2_pro/fqData && cd ~/ngs2_pro/fqData
fastq-dump SRR1839137 #data size(~10G)
fastq-dump SRR1838529 #data size(~7G)

#Download reference file
mkdir -p ~/ngs2_pro/sample_data && cd ~/ngs2_pro/sample_data
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

mkdir -p ~/ngs2_pro/sample_data1 && cd ~/ngs2_pro/sample_data1
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz

#download vcf files
wget ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr3.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-89/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz

#index chr1
mkdir -p ~/ngs2_pro/bwa_align/bwaIndex && cd ~/ngs2_pro/bwa_align/bwaIndex
ln -s ~/ngs2_pro/sample_data/Homo_sapiens.GRCh38.dna.chromosome.1.fa .
bwa index -a bwtsw Homo_sapiens.GRCh38.dna.chromosome.1.fa
#index chr3
#mkdir -p ~/ngs2_pro/bwa_align1/bwaIndex1 && cd ~/ngs2_pro/bwa_align1/bwaIndex1
#ln -s ~/ngs2_pro/sample_data1/ Homo_sapiens.GRCh38.dna.chromosome.3.fa .
#bwa index -a bwtsw Homo_sapiens.GRCh38.dna.chromosome.3.fa
   
#index whole
#mkdir -p ~/ngs2_pro/bwa_align/bwaIndex_whole && cd ~/ngs2_pro/bwa_align/bwaIndex_whole
#ln -s ~/ngs2_pro/sample_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa .
#bwa index -a bwtsw Homo_sapiens.GRCh38.dna.primary_assembly.fa

#count no of bases in each fastq file # prefer more than 100M bases per read group
cat SRR1839137.fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c #contains 3,414.9376 Mbases 
cat SRR1838529.fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c #contains 2,314.9568 Mbases

#Add [Read group information] and align all reads
mkdir -p ~/ngs/gatk && cd ~/ngs/gatk
for R1 in ~/ngs/fqData/SRR1839137.fastq;do
    SM=$(basename $R1 | cut -d"_" -f1)                                          
    LB=$(basename $R1 | cut -d"_" -f1,2)                                        
    PL="Illumina"                                                           
    RGID=$(cat $R1 | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       
    PU=$RGID.$LB                                                            
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"
    index="$HOME/ngs/bwa_align/bwaIndex/Homo_sapiens.GRCh38.dna.chromosome.1.fa"
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" $index $R1 > $(basename $R1 _SRR1839137.fastq).sam
done
#generate & sort BAM file
samtools view -hbo SRR1838529.bam SRR1838529.sam                # Convert the SAM file into a BAM file that can be sorted and indexed
samtools sort SRR1838529.bam -o SRR1838529.sorted.bam           # Sort the BAM file by position in genome
#samtools index  SRR1838529.sorted.bam                          # Index the BAM file so that we can randomly access it quickly

#mapping QC
for bamFile in *.sorted.bam;do
  output=${bamFile%.sorted.bam}
  samtools depth $bamFile | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > $output.cov
  samtools flagstat $bamFile > $output.stat #mapping quality=38.82%
done

# Install Picard tools
conda install -c bioconda picard 
picard_path=$CONDA_PREFIX/share/picard-2.19.0-0
#Mark duplicate
for sample in *.sorted.bam;do
  name=${sample%.sorted.bam}
  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done
##indexing
#sample
java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR1838529.dedup.bam

# Reference
ln -s ~/ngs/sample_data/Homo_sapiens.GRCh38.dna.chromosome.1.fa .
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.chromosome.1.fa O=Homo_sapiens.GRCh38.dna.chromosome.1.dict
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.1.fa

# Download known polymorphic sites
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr1_GRCh38.genotypes.20170504.vcf.gz -O homo_sapiens-chr1.vcf.gz
wget 'ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr1.vcf.gz' -O homo_sapiens-chr1.vcf.gz

# Select variants on chr5 and correct chr name
gunzip homo_sapiens-chr1.vcf.gz
grep "^#" homo_sapiens-chr1.vcf > homo_chr1.vcf
grep "^1" homo_sapiens-chr1.vcf | sed 's/^1/chr1/' >> homo_chr1.vcf 

#wc -l  homo_chr1.vcf #51576057
sed '138813d' homo_chr1.vcf > modi_homo_chr1.vcf
sed '497718d' modi_homo_chr1.vcf > homo_chr1.vcf
sed '1211016d' homo_chr1.vcf > homo1_chr1.vcf
sed '1471889d' homo1_chr1.vcf > homo2_chr1.vcf
sed '1670525d' homo2_chr1.vcf > homo3_chr1.vcf
sed '1706851d' homo3_chr1.vcf > homo4_chr1.vcf
sed '1850075d' homo5_chr1.vcf > homo6_chr1.vcf
sed '2028471d' homo6_chr1.vcf > homo7_chr1.vcf
sed '2054425d' homo7_chr1.vcf > homo8_chr1.vcf
sed '2329518d' homo8_chr1.vcf > homo9_chr1.vcf
sed '2426030d' homo9_chr1.vcf > homo10_chr1.vcf

gatk IndexFeatureFile -F homo10_chr1.vcf

#java -jar $picard_path/picard.jar LiftoverVcf -I homo_sapiens-chr1.vcf -O lifted_over.vcf -CHAIN b37tohg38.chain -REJECT rejected_variants.vcf -R Homo_sapiens.GRCh38.dna.primary_assembly.fa
