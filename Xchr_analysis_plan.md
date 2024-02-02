# X Chr Analysis Plan

## Section 1 : Analysis for QC of the X-chromosome

In this section, we will perform analysis that are helpful for quality control of the X-chromosome using PLINK(https://www.cog-genomics.org/plink/2.0/). The set of analysis needed are as follow: 
1. **Rate of heterozygote in males only** this is to investigate whether some SNPs have inds with heterozygote alleles (this is an indicator of potential genotyping or imputation errors)

```
   Plink2a \
   --bfile $PATH_TO_PLINK_FILES \
   --chr X \
   --covar $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate file to avoid auto-conversion auto assigment of sex in male ans female
   --geno-counts --remove-if "$SEX_Variable==2" --threads $NUMBER APPROPRIATE THREADS \ ## keep males only 
   --out $PATH_OUTPUT
```
   
3. **Test for difference in MAF between males and females controls only**

```
Plink2a \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate file to avoid auto-conversion auto assigment of sex in male ans female
--pheno-name $SEX_Variable  --mfilter 12 --assoc fisher \
--out $PATH_OUTPUT
```

4. **Test for differential missingness between males and females**

```
Plink2a \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate file to avoid auto-conversion auto assigment of sex in male ans female
--pheno-name $SEX_Variable  --test-missing \
--out $PATH_OUTPUT
```

5. **HWE test in females only**
   
```
Plink2a \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_PHENOTYPES \ ##
--pheno-name $CAD_Variable \
--keep-if "$CAD_Variable==1" | ## keep control only
--covar $PATH_SEX_OR_COVAR_INCLUDINGSEX
--keep-if "$SEX_Variable==2" \ ## keeep female only
--hardy \
--out $PATH_OUTPUT 
```

## Section 2: In this section we will perfom analyis of CAD for the X chromosome separately for each sex
The analysis can be done using PLINK, REGENIE or SAIGE. REGENIE and SAIGE allow the inclusion of related individuals while PLINK do not not. we will consider 2 different model :
1. Model 1: Activation of the X-chromosome. Here each SNPs is code as  0/1/2 in females and 0/1 in male
1.a. model for females (using PLINK)
1.b. Model for males

2. Model 2: Activation of the X-Chr. This model will be conducted in males only with alleles for each SNP code as 0/2 (assuming that 1 copy of the effect allele in males have the same effect as 2 copy in females)
    

## Section 3: Alternative script for analysis of the X-chr using REGENIE and SAIGE
This section provides examples scripts for X-chr analysis using both REGENIE and SAIGE

1. X-CHR analysis using REGENIE
2. X-chr analysis using SAIGE


## Section 4: description of the summary results to provide for the X-chromosome analysis
This section describe the summary results of the analysis performed for the X-chromosome
1. table content summary results for the QC
2. table content summary results association testing
   

