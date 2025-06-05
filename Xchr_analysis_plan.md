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
--pheno-name $SEX_Variable  --mfilter 12 --glm --xchr-model 1 \ ## --xchr-model 1 will help keep the association test as Male dosages are on a 0..1 scale on chrX, while females are 0..2 which is what  we want here 
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
--keep-females ## keep female only
--hardy \
--out $PATH_OUTPUT 
```

### Table with summary results of the QC step should be structure as follow:
**SNP**	SNP label for the variant	Identifier from the annotation file	rs693
			chr2:7819
			chr:pos:A1:A2

**CHR**	Chromosome on which SNP resides	Numeric for chromosomes 1-22; [current upload of autosomes only]	1

**POS** genome build	Position of SNP on chromosome	Base pairs on human genome build used	34000345

**EFFECT_ALLELE**	Allele at this site to which the effect has been estimated	Capital letter (A,C,G,T)	A

**NON_EFFECT_ALLELE**	Other allele at this site (please check software documentation before label-ling a1,a2 as non-effect allele)	Capital letter (A,C,G,T)	G

**EAF_in_males**	Allele frequency of the EF-FECT_ALLELE in males	Frequency with 3 digits to the right of the decimal	0.354

**EAF_in_females**	Allele frequency of the EF-FECT_ALLELE in females	Frequency with 3 digits to the right of the decimal	0.354

**pvalue_diff_EAF**	pvalue indicating hte difference in allelic frequency between males and females

**pvalue_diff_missingness**	pvalue indicating the difference in missingness between between males and females

**het_male**	value indicating the rate of heterozygote for each SNPs in males only (heterozygote rate should be run in males only)

**HWE**		pvalue for HWE test for each SNPs in females only (HWE shoud be run in female only)


## Section 2: In this section we will perfom analyis of CAD for the X chromosome separately for each sex
The analysis can be done using PLINK, REGENIE or SAIGE. REGENIE and SAIGE allow the inclusion of related individuals while PLINK do not not. 

### Example script for X-chromosome analysis using PLINK
***NOTE: While sharing your summary statistics, please structure the file name as follow:*** 

***Summary statistics from female: Xchr_female.[ancestry].[cohortname].[sex].[first_lastname].[phenotype].[analysissoftware].[date]****

***Summary statistics from male***: 
-	***Xchr_male_model1.plink.[ancestry].[cohortname].[first_lastname].[phenotype].[date]***
-	***Xchr_male_model2.[ancestry].[cohortname].[first_lastname].[phenotype].[analysissoftware].[date]***

#### Female only analysis  

Giving that females can be homozygote or heterozygote for each SNPS, female are analyzed in th esimilar fassion as heterozygote with a dosage model (0/1/2). the analysis here is similar to the autosomal analysis and can be achoeved with the model 1 or model 2 for x-chr analysis in plink.

  ```
    Plink2 \
    --pfile $PATH_TO_PLINK_FILES \
    --covar $PATH_TO_FILE_WITH_COVARIATE \
    --covar-name $PATH_SEX_OR_COVAR \
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
    --covar-name $PATH_SEX_OR_COVAR \
    --mac 10 \
    --mach-r2-filter 0.3 \
    --keep-males \
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

### Step 2 Association testing with WGS or Imputed data

#### Female-only analysis (coded as 0/1/2)
```
regenie \
--step 2 --bed Xchr --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_females \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_females_pred.list \  ## H0_females_pred.list is obtained in regenie step 1 from the sex stratified autosomal analysis
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out xchr_results_females
```
#### Male-only analysis 
For the X-chromosome, given that males can only be homozygous for either the risk allele or the non-risk allele, it is crucial to determine which model (activation or inactivation of the X-chr) accurately represents the disease risk at the gene/SNP level. In our pipeline, we will be testing both models.
  
#### Model 1: Activation of the X-chromosome. Here each SNPs 0/1 in male assuming that the effect of the SNPS on the disease is equivalent to what observed in a female heterozygote 
***Note: this model has just been recently implemented in REGENIE (REGENIE V4.1) https://github.com/rgcgithub/regenie/releases/tag/v4.1. Here to specifically code male as 0-1, the option ```--skip-dosage-comp``` should to be provide.***

```
regenie \
--step 2 --bed Xchr --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_males \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_males_pred.list \   ## H0_males_pred.list is obtained in regenie step 1 from the sex stratified autosomal analysis
--skip-dosage-comp \ ## With this option, male genotypes (or dosages) will be divided by 2 in non-PAR regions (i.e. male genotypes will be on a [0,1])
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out xchr_results_females
```


#### Model 2: Assuming Xchr anactivation (here male genotypes are coded as 0/2). 
***This is the default model for Xchr analysis in REGENIE. By default, REGENIE assumes males are coded as 0/2 in the non-PAR regions of the genotype file***
```
regenie \
--step 2 --bed Xchr --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_males \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_males_pred.list \   ## H0_males_pred.list is obtained in regenie step 1 from the sex stratified autosomal analysis
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out xchr_results_females
```

***Note from REGENIE sofware: To include X chromosome genotypes in step 1 and/or step 2, males should be coded as diploid so that their genotypes are 0/2 (this is done automatically for BED and PGEN file formats with haploid genotypes). Chromosome values of 23 (for human analyses), X, Y, XY, PAR1 and PAR2 are all acceptable and will be collapsed into a single chromosome. However, there is an updated version of the software (RGENIE V4.1) that allow a model in which male are code on 0-1 scale. to allow this model, the option ```--skip-dosage-comp``` shoould be provided. more indo at https://github.com/rgcgithub/regenie/releases/tag/v4.1*** 
   

### 2. table content summary results association testing for each sex: Below is the list of output column to have in the summary statistics for each sex. 

**SNP**	SNP label for the variant	Identifier from the annotation file	rs693
			chr2:7819
			chr:pos:A1:A2

**CHR**	Chromosome on which SNP resides	Numeric for chromosomes 1-22; [current upload of autosomes only]	1

**POS** genome build	Position of SNP on chromosome	Base pairs on human genome build used	34000345

**EFFECT_ALLELE**	Allele at this site to which the effect has been estimated	Capital letter (A,C,G,T)	A

**NON_EFFECT_ALLELE**	Other allele at this site (please check software documentation before label-ling a1,a2 as non-effect allele)	Capital letter (A,C,G,T)	G

**N_TOTAL**	Total number of cases and controls analyzed	Numeric, integer	1243

**N_CASES**	Total number of cases analyzed	Numeric, integer	1243

**N_CONTROLS**	Total number of controls analyzed	Numeric, integer	1243

**EAF_ALL**	Allele frequency of the EF-FECT_ALLELE in all	Frequency with 3 digits to the right of the decimal	0.354

**EAF_CASES**	Allele frequency of the EF-FECT_ALLELE in cases analyzed	Frequency with 3 digits to the right of the decimal	0.354

**EAF_CONTROLS**	Allele frequency of the EF-FECT_ALLELE in controls ana-lyzed	Frequency with 3 digits to the right of the decimal	0.354

**BETA**	Estimate of the allelic effect, defined as the natural logarithm of the odds ratio, ln(OR)	Numeric float with 3 digits to the right of the decimal	0.203

**SE**	Estimated standard error on the es-timate of the allelic effect, uncor-rected for genomic control	Numeric float with 4 digits to the right of the decimal	0.5611

**PVAL**	Significance of the variant associa-tion, uncorrected for genomic con-trol	Scientific E notation with 3 digits to the right of the decimal	3.24E-10

**INFO**	Measure of information content for the imputed SNP result (range 0-1) (autosomal only)	Numeric float with 3 digits to the right of the decimal (set to missing if genotyped)	0.483

   

