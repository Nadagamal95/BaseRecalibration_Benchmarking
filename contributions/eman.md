#contribution
``` java
1- I choose a model organism "_E.coli_" but not all team members approved this choice so we discard it.
- i provide resourses such as;
https://gatkforums.broadinstitute.org/gatk/discussions/tagged/bacteria/p1
https://github.com/ekg/alignment-and-variant-calling-tutorial/blob/master/README.md
https://gatkforums.broadinstitute.org/gatk/discussion/9790/haplotypecaller-ploidy-for-pooled-samples-of-bacterial-genomes
https://software.broadinstitute.org/gatk/documentation/article?id=10913

2- after we agreed on GATK germline pipline and whole genome for alignment step, i curated the data set they choose as it was targeted seq against chr21 trisomy only. so I choose another dataset and we agreed on it:
paper:
http://www.bloodjournal.org/content/127/21/2598.long?sso-checked=true#abstract-2
data set:
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP053196 #we agreed on first SRR1839137
```



_ _ _ _ _ _ _ _ _ _ _ _ _ _ 

## Troubleshooting
```java
1- error1: trying to index chr1 & 2 with BWA 
two many levels of symbolic links
```
```java
2- error2:
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.chromosome.1.fa O=Homo_sapiens.GRCh38.dna.chromosome.1.dict
```
```java
3- error3:
Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 138813: Duplicate allele added to VariantContext: CTTTCTTTCTTTCT
```
``` java
4-error4:
A USER ERROR has occurred: Error while trying to create index for /home/emmyneutron/ngs/indexing/homo6_chr1.vcf. Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number 7920778: empty alleles are not permitted in VCF records
```














