Contribution
``` java
1- I choose a model organism "_E.coli_" but not all team members approved this choice so we discard it.
- i provide resourses such as;
https://gatkforums.broadinstitute.org/gatk/discussions/tagged/bacteria/p1
https://github.com/ekg/alignment-and-variant-calling-tutorial/blob/master/README.md
https://gatkforums.broadinstitute.org/gatk/discussion/9790/haplotypecaller-ploidy-for-pooled-samples-of-bacterial-genomes
https://software.broadinstitute.org/gatk/documentation/article?id=10913

2- after we agreed on GATK germline pipline and whole genome for alignment step,so I choose another dataset(whole exome from germ line DNA) :
paper:
http://www.bloodjournal.org/content/127/21/2598.long?sso-checked=true#abstract-2
data set:
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP053196 #we agreed on first SRR1839137
workflow:
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

- i choose two samples and i aligned them against chr1 & chr3 "GRCh38",(due to computational limitation)
sample1:SRR1839137 #data size(~10G) #affected female
sample2:SRR1838529 #data size(~7G)  #unaffected female

3- i read bench marking paper to set inclusion criteria
paper:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006494
```



_ _ _ _ _ _ _ _ _ _ _ _ _ _ 

## Troubleshooting
```java
1- error1: trying to index chr1 & 2 with BWA 
two many levels of symbolic links

- i make each index alone and remove to external hard till i need it again
```
```java
2- error2:
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.chromosome.1.fa O=Homo_sapiens.GRCh38.dna.chromosome.1.dict

- i modified path:
picard_path=$CONDA_PREFIX/share/picard-2.19.0-0
```
```java
3- error3:
Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 138813: Duplicate allele added to VariantContext: CTTTCTTTCTTTCT

- sed '138813' homo_chr1.vcf > homo1_chr1.vcf

https://gatkforums.broadinstitute.org/gatk/discussion/11547/error-duplicate-allele-added-to-variantcontext
https://www.biostars.org/p/228915/
https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_vcf_LiftoverVcf.php
```
``` java
4-error4:
A USER ERROR has occurred: Error while trying to create index for /home/emmyneutron/ngs/indexing/homo6_chr1.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 7920778: empty alleles are not permitted in VCF records

- sed '7920778d' homo1_chr1.vcf > modi_homo_chr1.vcf
https://www.biostars.org/p/243810/
```














