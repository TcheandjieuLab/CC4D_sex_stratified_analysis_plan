#  Analysis plan for sex stratify autosomal chromosomes

the purpose of this pipeline is to provide detail on sex stratif analysis for CAD . 
Please test additive models using logistic regression, accounting for genotype imputation uncertainty (i.e. SNP probability or dosage).

## Section 1: Analysis pipeline
Primary model here will be:
a. Males: CAD(all CAD cases vs. all non-CAD controls) = $SNP + Study covariates 
b. Females: CAD (all CAD cases vs. all non-CAD controls) = $SNP + Study covariates 
c. The should include the following filtering criteria: SNPs with imputation quality (INFO) > 0.3 and Mininum Allele Count (MAC)=10![image](https://github.com/TcheandjieuLab/CC4D_sex_stratified_analysis_plan/assets/32551968/c795e1f7-6cae-4fe5-acbe-466a5c6c2cb3)


### 1. Example script using PLINK

This is a single example command line to run GWAS with PLINK when chromosome are strore in seprate files

1.a. Analysis in males

```
for chr in {1..22}; do
 plink2 \
  --pfile Imputed_files.diresctory/filesmale_imp_chr${chr} \
  --pheno phenotypes_files \
  --covar covariate_files.txt \
  --covar-name PC1-PC5 age \
  --glm hide-covar cols=+a1freq,+machr2  \
  --mac 10 \
  --keep-if Sex==1 \ # Analysis for males
  --threads 20  \
  --mach-r2-filter 0.3  --mac 10 --hwe 1e-10 \
  --memory 42000 --covar-variance-standardize \
  --pheno-name CAD \
  --out male.output.fileschr${chr}.txt
  ```
1.b. Analysis in females

```
for chr in {1..22}; do
 plink2 \
  --pfile Imputed_files.diresctory/filesmale_imp_chr${chr} \
  --pheno phenotypes_files \
  --covar covariate_files.txt \
  --covar-name PC1-PC5 age \
  --glm hide-covar cols=+a1freq,+machr2  \
  --mac 10 \
  --keep-if Sex==2 \ # Analysis for males
  --threads 20  \
  --mach-r2-filter 0.3  --mac 10 --hwe 1e-10 \
  --memory 42000 --covar-variance-standardize \
  --pheno-name CAD \
  --out females.output.fileschr${chr}.txt
  ```
### 2. Example script using REGENIE [module load Regenie/v2.0.1 (users should have regenie installed. I am using v2.0.1)]

Regenie runs in two steps. Before starting the first step, one needs to generate a plink file containing a subset of independent SNPs (e.g., 100K as recommended by Regenie). This file will be used to run H0 and account for relatedness between subjects.
In the first step, 100K.bed 100K.bim 100K.fam will be the files we will use (file containing all inds and 100k SNPs).

#### Step 1: Perform H0. Two runs will be performed, one for females and one for males. We need to create two files that contain FID and IID of females and males (i.e., id_females and id_males, with header (FID IID)).
```
## H0 in Females
regenie \
--step 1 --bed 100K   \
--covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_females --bsize 10000 --bt --lowmem --lowmem-prefix tmp_rg1 \
--out H0_females

## H0 in Males
regenie \
--step 1 --bed 100K   --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_males --bsize 10000 --bt --lowmem --lowmem-prefix tmp_rg1 \
--out H0_males
```

#### Step 2: Perform H1. This step runs the association test for SNVs. In the following case, chromosomes were aggregated in one plink file. If users have data split by chromosome, a loop needs to be used to run all chromosomes.
```
# H1 in Females
regenie \
--step 2 --bed AllChr  --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_females \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_females_pred.list \
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out results_females

# H1 in Males
regenie \
--step 2 --bed AllChr  --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age\
--phenoFile pheno.txt --phenoCol CAD --keep id_males \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_males_pred.list \
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out results_males
```

**If your genomic data is split by chromosome, you can create use a loop style function to run all chrs with a single script as follow:**
```
for i in `seq 1 22`
do
# H1 in Females
regenie \
--step 2 --bed Chr"$i"  --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_females \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_females_pred.list \
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out results_females_chr"$i"

# H1 in Males
regenie \
--step 2 --bed Chr"$i"  --covarFile pheno.txt --covarCol PC1,PC2,PC3,PC4,PC5,Age \
--phenoFile pheno.txt --phenoCol CAD --keep id_males \
--bsize 10000 --bt --firth --approx --pThresh 0.05 --pred H0_males_pred.list \
--threads $i  #multithreads with a number appropriate to your cluster \
--minMAC 10 --minINFO 0.3 --af-cc \
--out results_males_chr"$i"
done
```
#### File formats #### 
pheno.txt: FID IID CAD PC1 PC2 PC3 PC4 PC5 Sex Age

id_females: FID IID

id_males: FID IID

***The file is space-delimited and can have many more covariates. Users have the flexibility to use the set of covariates in the regenie command-line above.***

### 3. Example script using SAIGE

### 4. List of variables to provide in the summary statistics for each sex specific GWAS
**SNP**:	SNP label for the variant	Identifier from the annotation file. Example: rs693 (chr2:7819; chr:pos:A1:A2)

**CHR**:	Chromosome on which SNP resides	Numeric for chromosomes 1-22. Example:	1

**POS with genome build**:	Position of SNP on chromosome	Base pairs on human genome build used. Example:	34000345

**EFFECT_ALLELE**:	Allele at this site to which the effect has been estimated	Capital letter (A,C,G,T). Example:	A

**NON_EFFECT_ALLELE**:	Other allele at this site (please check software documentation before label-ling a1,a2 as non-effect allele)	Capital letter (A,C,G,T). Example: G

**N_TOTAL**:	Total number of cases and controls analyzed	Numeric, integer	Example: 1243

**N_CASES**:	Total number of cases analyzed	(Numeric value)

**N_CONTROLS**:	Total number of controls analyzed	(Numeric value)

**EAF_ALL**:	Allele frequency of the EF-FECT_ALLELE in all;	frequency with 3 digits to the right of the decimal. Example: 0.354

**EAF_CASES**:	Allele frequency of the EF-FECT_ALLELE in cases analyzed;	frequency with 3 digits to the right of the decimal. Example: 0.354

**EAF_CONTROLS**:	Allele frequency of the EF-FECT_ALLELE in controls ana-lyzed;	frequency with 3 digits to the right of the decimal. Example: 0.354

**BETA**:	Estimate of the allelic effect, defined as the natural logarithm of the odds ratio, ln(OR);	numeric float with 3 digits to the right of the decimal. Example: 0.203

**SE**:	Estimated standard error on the es-timate of the allelic effect, uncor-rected for genomic control;	numeric float with 4 digits to the right of the decimal. Example:	0.5611

**PVAL**:	Significance of the variant associa-tion, uncorrected for genomic con-trol;	scientific E notation with 3 digits to the right of the decimal.	Example: 3.24E-10

**INFO**:	Measure of information content for the imputed SNP result (range 0-1) (autosomal only);	numeric float with 3 digits to the right of the decimal (set to missing if genotyped)	Example: 0.483

**HWE**: (should be used on the best guest genotypes, can be obtain from Plink independently):	P-value of the HWE test. Example: 3.00E-06
![image](https://github.com/TcheandjieuLab/CC4D_sex_stratified_analysis_plan/assets/32551968/8f09b92d-9481-4f83-8a97-18d3047bc5fe)

## Section 2: sex interaction model (will be following the analytical pipeline from GLGC)

the sex-interaction test will be dome using GEM (https://github.com/large-scale-gxe-methods/GEM). 
1. This analysis should be conducted using **Unrealted individual only**.
2. For consistency in effect estimate report, We will have Male and Female with the same codage across all studies:
   
    **a. Male coded as 1**
   
    **b. Female coded as 2**

**Example of input/script for sex interaction test using GEM**

```
IID      FID      CAD  age  array  sex   PC1     PC2      PC3      PC4      PC5     ...   PC10
1000001  1000001  0    44   0      0    -0.2283  0.0986  -0.1185  -0.0881  -1.0568  ...   1.2153
1000002  1000002  0    69   0      1    -0.2655  0.1417  -0.1286   0.5837  -0.3744  ...  -0.4781
1000003  1000003  1    66   0      0    -0.2033  0.2145  -0.1035  -0.3426  -0.7164  ...  -0.2142
1000004  1000004  1    54   0      1    -0.1658  0.1421   0.0258  -0.1435   0.2450  ...   0.3132
```

```
./GEM_1.5.2_Intel \
    --bgen ukb_imp_chr${chr}_v3.bgen \
    --sample ${sample} \
    --pheno-file UKBB_CAD.${pop}.pheno.txt \
    --delim '\t' \
    --pheno-name CAD \
    --sampleid-name IID \
    --exposure-name sex \
    --categorical-names sex \
    --covar-names age array PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
    --robust 1 \
    --center 0 \
    --scale 0 \
    --threads 8 \
    --output-style full \
    --out UKBB_CAD.${pop}.${chr}.gem
```

### list of variables to provide in the summary statistic 


***Note: For studies that have already conducted analyses, please discuss the models used (for example, some studies may have already adjusted for age and this is not considered a material deviation from the analysis plan). Use study appropriate software to account for (or exclude as appropriate) relatedness.***
***If undertaking new analyses to contribute then the following approach is advised. If providing previously generated results, please ensure the approach is documented and discuss any major deviations***

