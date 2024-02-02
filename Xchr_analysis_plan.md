# X Chr Analysis Plan

## Step 1: Analysis that are helpfull for QC of the X-chromosome

In this section, we will perform analysis that are helpful for quality control of the X-chromosome using PLINK(https://www.cog-genomics.org/plink/2.0/). The set of analysis needed are as follow: 
1. Rate of heterozygote in males only to investigate whether some SNPs have inds with heterozygote alleles (this is an indicator of potential genotyping or imputation errors)
   Plink2a \
   --bfile $PATH_TO_PLINK_FILES \
   --chr X \
   --covar $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate file to avoid auto-conversion auto assigment of sex in male ans female
   --geno-counts --remove-if sex==2 --threads $NUMBER APPROPRIATE THREADS \ 
   --out $PATH_OUTPUT 
3. SNPs allelic frequencies
   in fales
   In females
   test difference in allelic frequencies between males and females
4. HWE test in females only
5. 

* Foo
* Bar
* Baz

## Section 2

### Section 2.1

### Section 2.2

## Section 3
