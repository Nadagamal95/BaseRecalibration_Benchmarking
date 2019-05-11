#install gatk4

conda install -c bioconda gatk4 

##Download known varinats

# Download known polymorphic sites

mkdir -p ~/ngs2_project/VCF && cd ~/ngs2_project/VCF
wget 'ftp://ftp.ensembl.org/pub/release-89/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz' -O Homo_sapiens.vcf.gz
gunzip Homo_sapiens.vcf.gz

#gatk IndexFeatureFile -F Homo_sapiens.vcf

3error in indexing the vcf files 
#A USER ERROR has occurred: Error while trying to create index for /home/ngs-01/ngs2_project/VCF/Homo_sapiens.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 14625537: unparsable vcf record with allele W

#I'm trying to solve the error 
#https://www.biostars.org/p/245848/
#https://samtools.github.io/hts-specs/VCFv4.2.pdf

#this command to visualize the place of my error to solve it

sed -n '14625537p' Homo_sapiens.vcf  

10	52771028	rs17287498	C	W,A,G,T	.	.	dbSNP_149;TSA=SNV;E_Freq;E_1000G;E_Cited;Reference_seq=C;Variant_seq=W,A,G,T;variation_id=8772965;allele_string=C,W,A,G,T;AA=A


#trying to substitute the "W" letter by "N" 

#sed 's/W/N/' Homo_sapiens.vcf >> Homo_sapiens2.vcf 
#does not work out

#sed 's/W/N/g' Homo_adjusted.vcf > Homo_corrected.vcf
#does not work out

#sed '/W/c N ' Homo_adjusted.vcf > Homo_corrected.vcf
#but it does not work

###############################

# correcting all the chromosomes in genome have been corrected from numbers to (chr) + (chr.no.)

grep "^#" Homo_sapiens.vcf > Homo_adjusted.vcf
perl -pe 's/^([^#])/chr\1/' Homo_sapiens.vcf >> Homo_adjusted.vcf

# for testing and cheking that all the chromosomes in genome have been corrected from numbers to (chr) + (chr.no.)
grep -w "^chr5" Homo_adjusted.vcf

#gatk IndexFeatureFile -F Human_corrected.vcf

#again the same error 

#A USER ERROR has occurred: Error while trying to create index for /home/ngs-01/ngs2_project/VCF/Homo_adjusted.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 14625537: unparsable vcf record with allele W

#perl -pe 's/^([W])/N\1/' Homo_adjusted.vcf >> Homo_adjusted2.vcf
# it does not work also

###############################
#grep "^#" Homo_sapiens.vcf > Human_corrected.vcf
grep "W" Homo_adjusted.vcf | sed 's/W/N/' >> Human_corrected.vcf
####################################

#grepping the lines that has the char "W" to remove them
grep -nr "W" Homo_adjusted.vcf

#removing the position that has this char which is what makes the error 
sed -i '14625549d' Homo_adjusted.vcf

#indexing the vcf files after removing the error
gatk IndexFeatureFile -F Homo_adjusted.vcf

#it works but another error appears 

A USER ERROR has occurred: Error while trying to create index for /home/ngs-01/ngs2_project/VCF/Homo_adjusted.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 19332647: empty alleles are not permitted in VCF records

#awk -F '\t' '($0 ~ /^#/ || $5!=".")' Homo_adjusted.vcf > Homo_adjusted2.vcf
#it does not work, same error 

sudo apt install gawk
gawk 'BEGIN{FS="\t"; OFS="\t"}{if (NF>1 && $5=="") {$5="."; print $0} else print $0}' Homo_adjusted.vcf > Homo_adjusted2.vcf

#gatk IndexFeatureFile -F Homo_adjusted2.vcf 
#it works 

# another error in the next chr 
The provided VCF file is malformed at approximately line number 28156388: unparsable vcf record with allele H

grep -nr "H" Homo_adjusted2.vcf

sed -i '28156388d' Homo_adjusted2.vcf
