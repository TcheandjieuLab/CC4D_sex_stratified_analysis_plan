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
The analysis can be done using PLINK, REGENIE or SAIGE. REGENIE and SAIGE allow the inclusion of related individuals while PLINK do not not. 

### Example script for X-chromosome analysis using PLINK

#### Female only analysis  

Giving that females can be homozygote or heterozygote for each SNPS, female are analyzed in th esimilar fassion as heterozygote with a dosage model (0/1/2). the analysis here is similar to the autosomal analysis and can be achoeved with the model 1 or model 2 for x-chr analysis in plink.

  ```
    Plink2 \
    --pfile $PATH_TO_PLINK_FILES \
    --covar $PATH_TO_FILE_WITH_COVARIATE \
    --covar-name $PATH_SEX_OR_COVAR_INCLUDINGSEX \
    --mac 10 \
    --mach-r2-filter 0.3 \
    --keep-females \
    --xchr-model 2 \ # model 2 is the default model in plink2
    --pheno $PATH_TO_FILE_WITH_PHENOTYPES \
    --pheno-name $CAD_Variable \
    --threads 6 \
    --glm hide-covar firth-fallback  cols=+a1countcc,+a1freqcc,+machr2,+totallelecc,+nobs \ ## firth-fallback  glm fall on firth regression if low case number 
    --remove $PATH_TO_SUBJECT_to_exclude  \  ## this can be a list of related ind that should be excluded from the model
    --out $PATH_OUTPUT_FEMALE ## path to the output summary statistics 
  ```
#### Male only analysis  
For the X-chromosome, given that males can only be homozygous for either the risk allele or the non-risk allele, it is crucial to determine which model (activation or inactivation of the X-chr) accurately represents the disease risk at the gene/SNP level. In our pipeline, we will be testing both models.
  
##### Model 1: Activation of the X-chromosome. Here each SNPs 0/1 in male assuming that the effect of the SNPS on the disease is equivalent to what observed in a female heterozygote
   
  ```
    Plink2 \
    --pfile $PATH_TO_PLINK_FILES \
    --covar $PATH_TO_FILE_WITH_COVARIATE \
    --covar-name $PATH_SEX_OR_COVAR_INCLUDINGSEX \
    --mac 10 \
    --mach-r2-filter 0.3 \
    --keep-males \
    --xchr-model 1 \ ## For model 1 in PLINK2, Male dosages are on a 0..1 scale on chrX, while females are 0.1.2 scale
    --pheno $PATH_TO_FILE_WITH_PHENOTYPES \
    --pheno-name $CAD_Variable \
    --threads $number_of_thread \ ## the number of thread here should be adapted to the computing system used (6 is often well tolerated)
    --glm hide-covar firth-fallback  cols=+a1countcc,+a1freqcc,+machr2,+totallelecc,+nobs \ ## firth-fallback  glm fall on firth regression if low case number 
    --remove $PATH_TO_SUBJECT_to_exclude  \  ## this can be a list of related ind that should be excluded from the model
    --out $PATH_OUTPUT_MALES ## path to the output summary statistics 
   ```

##### Model 2: Inactivation of the X-Chr. Here each SNPs 0/2 in male assuming that the effect of the SNPS on the disease in male are equivalent to what observed in a female homozygot at risk (this assume 2 copy of the effect allele in males)
  
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
    --xchr-model 2 \ ## for model 2 in PLink Males and females are both on a 0..2 scale on chrX 
    --glm hide-covar firth-fallback  cols=+a1countcc,+a1freqcc,+machr2,+totallelecc,+nobs \ ## firth-fallback  glm fall on firth regression if low case number 
    --remove $PATH_TO_SUBJECT_to_exclude  \  ## this can be a list of related ind that should be excluded from the model
    --out $PATH_OUTPUT_MALES ## path to the output summary statistics 
  ```

***Note: Model 2 is the default choice for X-chromosome analysis using PLINK2(https://www.cog-genomics.org/plink/2.0/assoc#glm). However, in a sex-stratified analysis focusing on females, Models 1 and 2 will yield identical results. The distinction in output results is relevant only for males, who are coded as 0/1 in Model 1 and 0/2 in Model 2.*** 

### Alternative script for analysis of the X-chr using REGENIE 

The current version of Regenie is structured similarly to the default version of the X-chromosome analysis in PLINK2. In this model (Model 2), males are coded as carrying either 0 or 2 copies of the effect allele, while females are coded as carrying 0, 1, or 2 copies of the effect allele. Consequently, only the inactivation of the X-chr model will be considered when analyzing males using Regenie. bolow is an example script for running regenie for x-chr.

#### Step 1: genome wide regression model is fit to the traits. The This step the same step H0 when performed for the sex stratified autosomal chromosome. Thus output generated for this step during sex stratified autosomal analysis can be used here instead
##### H0 in Females

```
regenie \
--step 1 --bed 100K   \  ## bed file with 100K independent SNVs on autosomal chromosomes
--covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_females --bsize 10000 --bt --lowmem --lowmem-prefix tmp_rg1 \
--out H0_females
```
##### H0 in Males

```
regenie \
--step 1 --bed 100K   \ ## bed file with 100K independent SNVs on autosomal chromosomes
--covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_males --bsize 10000 --bt --lowmem --lowmem-prefix tmp_rg1 \
--out H0_males
```

#### Step 2 Association testing with WGS or Imputed data

##### Female-only analysis (coded as 0/1/2)
```
regenie \
--step 2 --bed Xchr --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_females \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_females_pred.list \  ## H0_females_pred.list is obtained in regenie step 1 from the sex stratified autosomal analysis
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out xchr_results_females
```
##### Male-only analysis (coded as 0/2) corresponding to x-chr inactivation
```
regenie \
--step 2 --bed Xchr --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_males \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_males_pred.list \   ## H0_males_pred.list is obtained in regenie step 1 from the sex stratified autosomal analysis
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out xchr_results_females
```

***note from REGENIE sofware: To include X chromosome genotypes in step 1 and/or step 2, males should be coded as diploid so that their genotypes are 0/2 (this is done automatically for BED and PGEN file formats with haploid genotypes). Chromosome values of 23 (for human analyses), X, Y, XY, PAR1 and PAR2 are all acceptable and will be collapsed into a single chromosome*** 
   
## Section 3: description of the summary results to provide for the X-chromosome analysis
This section describe the summary results of the analysis performed for the X-chromosome

1. table content summary results for the QC

2. table content summary results association testing for each sex
https://private-user-images.githubusercontent.com/32551968/305191704-8f09b92d-9481-4f83-8a97-18d3047bc5fe.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MjQ5NjMyMjksIm5iZiI6MTcyNDk2MjkyOSwicGF0aCI6Ii8zMjU1MTk2OC8zMDUxOTE3MDQtOGYwOWI5MmQtOTQ4MS00ZjgzLThhOTctMThkMzA0N2JjNWZlLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNDA4MjklMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjQwODI5VDIwMjIwOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWViNTUzMTM5NTAzM2RmMTU3NGIyYjg0N2RmZDk0Yjc0Y2FmMjAxNzUxM2Y3MWViMWM5MzQyNmQxMzFlM2JiOWUmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.jtauyY5pVAA2AwO5AFdart1OgwargkmnV47jJIvVDLY
   

