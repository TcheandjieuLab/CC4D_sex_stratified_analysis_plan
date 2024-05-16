# X Chr Analysis Plan

## Section 1 : Analysis for QC of the X-chromosome

In this section, we will perform analysis that are helpful for quality control of the X-chromosome using PLINK(https://www.cog-genomics.org/plink/2.0/). The set of analysis needed are as follow: 
1. **Rate of heterozygote in males only** this is to investigate whether some SNPs have inds with heterozygote alleles (this is an indicator of potential genotyping or imputation errors). #MS: what is the threshold we should use to filter out bad SNPs?
   **Important note: If sex is provided as part of the genotype data (e.g., in the .fam file), plink will auto assign male and female. This will then cause the sofware to systematically estimate the rate of heterozygozity for male to be zero and prevent us from catching genotyping errors. To avoid this, the sex should be provided as a covariate file separately, and should be removed put as missing in the .fam file** 

```
   Plink2a \
   --bfile $PATH_TO_PLINK_FILES \
   --chr X \
   --covar $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate covariatefile to avoid auto-conversion auto assigment of sex in male and female
   --geno-counts --remove-if "$SEX_Variable==2" --threads $NUMBER APPROPRIATE THREADS \ ## keep males only 
   --out $PATH_OUTPUT
```
#MS: 
The rate of heterozygosity for males can be computed as: HAP_ALT_CTS / (HAP_REF_CT + HAP_ALT_CTS) in the .gcount plink output file.
An alternative way to obtain heterozygosity in males with plink 1.9: use the "--hardy" feature combined with "--keep" to use males only
```
plink --bfile genotype_data_xchr --hardy --keep id_males --out hwe_males   ## id_males is a 2-col file that has FID and IID of the male individuals with a header. 
sed -i 's/\// /g' hwe_males.hwe                                            ## to split the three genotypes separated by "/"
sed -i 's/GENO/AA AB BB/g' hwe_males.hwe                                   ## To create three headers (AA, AB, BB) instead of one (GENO)
```
Here, heterozygosity rate can be calculated as AB / (AA + AB + BB)

   
3. **Test for difference in MAF between males and females controls only**
    **Important note: If sex is provide as part of the genotype data, plink will auto assign male and female. This will then cause the sofware to systematically correct the MAF in male and prevent us from catching genotyping errors. To avoid this, we sex should be provide as a covariate file separately.** 

```
Plink2a \
--bfile  $PATH_TO_PLINK_FILES \
--pheno $PATH_TO_FILE_WITH_SEX \ ## sex should be provide as a separate covariate file to avoid auto-conversion auto assigment of sex in male and female
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

5. **HWE test in female controls only**
   
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
#MS: An alternative way using plink1 1.9:
```
plink --bfile genotype_data_xchr --hardy --keep id_females_controls --out hwe_females_controls   ## id_females_controls contains female controls only
```
What is the threshold we should use to filter SNPs? P<10^-8? 



## Section 2: In this section we will perfom analyis of CAD for the X chromosome separately for each sex
The analysis can be done using PLINK, REGENIE or SAIGE. REGENIE and SAIGE allow the inclusion of related individuals while PLINK do not not. we will consider 2 different model :
1. Model 1: Activation of the X-chromosome. Here each SNPs is code as  0/1/2 in females and 0/1 in male
1.a. model for females (using PLINK)
1.b. Model for males

2. Model 2: Inactivation of the X-Chr. This model will be conducted in males only with alleles for each SNP code as 0/2 (assuming that 1 copy of the effect allele in males have the same effect as 2 copy in females)
    

## Section 3: Alternative script for analysis of the X-chr using REGENIE and SAIGE
This section provides examples scripts for X-chr analysis using both REGENIE and SAIGE

1. X-CHR analysis using REGENIE
2. X-chr analysis using SAIGE (if cohort are using SAIGE already)


## Section 4: description of the summary results to provide for the X-chromosome analysis
This section describe the summary results of the analysis performed for the X-chromosome
1. table content summary results for the QC
2. table content summary results association testing
   

