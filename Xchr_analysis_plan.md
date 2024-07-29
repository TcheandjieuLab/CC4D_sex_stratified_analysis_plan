# X Chr Analysis Plan

## Section 1 : Analysis for QC of the X-chromosome

In this section, we will perform analysis that are helpful for quality control of the X-chromosome using PLINK(https://www.cog-genomics.org/plink/2.0/). The set of analysis needed are as follow: 
1. **Rate of heterozygote in males only** this is to investigate whether some SNPs have inds with heterozygote alleles (this is an indicator of potential genotyping or imputation errors).
   **Important note: If sex is provided as part of the genotype data (e.g., in the .fam file), plink will auto assign male and female. This will then cause the sofware to systematically estimate the rate of heterozygozity for male to be zero and prevent us from catching genotyping errors. To avoid this, the sex should be provided as a covariate file separately, and should be removed put as missing in the .fam file** 

```
   Plink2 \
   --bfile $PATH_TO_PLINK_FILES \
   --chr X \
   --covar $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate covariatefile to avoid auto-conversion auto assigment of sex in male and female
   --geno-counts --remove-if "$SEX_Variable==2" --threads $NUMBER APPROPRIATE THREADS \ ## keep males only 
   --out $PATH_OUTPUT
```

   
3. **Test for difference in MAF between males and females controls only**
**Important note: If sex is provide as part of the genotype data, plink will auto assign male and female. This will then cause the sofware to systematically correct the MAF in male and prevent us from catching genotyping errors. To avoid this, we sex should be provide as a covariate file separately.** 

```
Plink2 \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate covariate file to avoid auto-conversion auto assigment of sex in male and female
--pheno-name $SEX_Variable  --mfilter 12 --assoc fisher \
--out $PATH_OUTPUT
```

4. **Test for differential missingness between males and females**

```
Plink2 \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate file to avoid auto-conversion auto assigment of sex in male ans female
--pheno-name $SEX_Variable  --test-missing \
--out $PATH_OUTPUT
```

5. **HWE test in female controls only**
   
```
Plink2 \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_PHENOTYPES \ ## path to the phenotypes files
--pheno-name $CAD_Variable \
--keep-if "$CAD_Variable==1" | ## keep control only
--covar $PATH_SEX_OR_COVAR_INCLUDINGSEX
--keep-if "$SEX_Variable==2" \ ## keeep female only
--hardy \
--out $PATH_OUTPUT 
```



## Section 2: In this section we will perfom analyis of CAD for the X chromosome separately for each sex
The analysis can be done using PLINK, REGENIE or SAIGE. REGENIE and SAIGE allow the inclusion of related individuals while PLINK do not not. we will consider 2 different model :

### Example script for X-chromosome analysis using PLINK
Model 1: Activation of the X-chromosome. Here each SNPs is code as  0/1/2 in females and 0/1 in male
   
   1.a Model for female-only 

       ```
    Plink2 \
    --pfile $PATH_TO_PLINK_FILES \
    --covar $PATH_TO_FILE_WITH_COVARIATE \
    --covar-name $PATH_SEX_OR_COVAR_INCLUDINGSEX \
    --mac 10 \
    --mach-r2-filter 0.3 \
    --keep-females \
    --xchr-model 1 \ ## model 1 represent the xchr activation
    --pheno $PATH_TO_FILE_WITH_PHENOTYPES \
    --pheno-name $CAD_Variable \
    --threads 6 \
    --glm hide-covar firth-fallback  cols=+a1countcc,+a1freqcc,+machr2,+totallelecc,+nobs \ ## firth-fallback  glm fall on firth regression if low case number 
    --remove $PATH_TO_SUBJECT_to_exclude  \  ## this can be a list of related ind that should be excluded from the model
    --out $PATH_OUTPUT_FEMALE ## path to the output summary statistics 
    ```
   1.b male-only analysis
   
  ```
    Plink2 \
    --pfile $PATH_TO_PLINK_FILES \
    --covar $PATH_TO_FILE_WITH_COVARIATE \
    --covar-name $PATH_SEX_OR_COVAR_INCLUDINGSEX \
    --mac 10 \
    --mach-r2-filter 0.3 \
    --keep-males \
    --xchr-model 1 \ ## model 1 represent the xchr activation
    --pheno $PATH_TO_FILE_WITH_PHENOTYPES \
    --pheno-name $CAD_Variable \
    --threads $number_of_thread \ ## the number of thread here should be adapted to the computing system used (6 is often well tolerated)
    --glm hide-covar firth-fallback  cols=+a1countcc,+a1freqcc,+machr2,+totallelecc,+nobs \ ## firth-fallback  glm fall on firth regression if low case number 
    --remove $PATH_TO_SUBJECT_to_exclude  \  ## this can be a list of related ind that should be excluded from the model
    --out $PATH_OUTPUT_MALES ## path to the output summary statistics 
    ```
Model 2: Inactivation of the X-Chr. This model will be conducted in females only with alleles for each SNP code as 0/2 (assuming that 1 copy of the effect allele in males have the same effect as 2 copy in females)

  ```
    Plink2 \
    --pfile $PATH_TO_PLINK_FILES \
    --covar $PATH_TO_FILE_WITH_COVARIATE \
    --covar-name $PATH_SEX_OR_COVAR_INCLUDINGSEX \
    --mac 10 \
    --mach-r2-filter 0.3 \
    --keep-females \
    --pheno $PATH_TO_FILE_WITH_PHENOTYPES \
    --pheno-name $CAD_Variable \
    --threads $number_of_thread \ ## the number of thread here should be adapted to the computing system used (6 is often well tolerated)
    --xchr-model 2 \ ## model 2 represent the xchr inactivation
    --glm hide-covar firth-fallback  cols=+a1countcc,+a1freqcc,+machr2,+totallelecc,+nobs \ ## firth-fallback  glm fall on firth regression if low case number 
    --remove $PATH_TO_SUBJECT_to_exclude  \  ## this can be a list of related ind that should be excluded from the model
    --out $PATH_OUTPUT_MALES ## path to the output summary statistics 
    ```

## Alternative script for analysis of the X-chr using REGENIE and SAIGE
This section provides examples scripts for X-chr analysis using both REGENIE and SAIGE

1. X-CHR analysis using REGENIE: 
Since REGENIE do not currently support X chromosome activation model, it can only be used to run model 2 for both Male and female.  

***note from REGENIE sofware: To include X chromosome genotypes in step 1 and/or step 2, males should be coded as diploid so that their genotypes are 0/2 (this is done automatically for BED and PGEN file formats with haploid genotypes). Chromosome values of 23 (for human analyses), X, Y, XY, PAR1 and PAR2 are all acceptable and will be collapsed into a single chromosome*** 
   

3. X-chr analysis using SAIGE (if cohort are using SAIGE already)


## Section 3: description of the summary results to provide for the X-chromosome analysis
This section describe the summary results of the analysis performed for the X-chromosome

1. table content summary results for the QC

2. table content summary results association testing
   

