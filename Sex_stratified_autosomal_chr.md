#  Analysis Plan sex stratify Autosomal chr

the purpose of this pipeline is to provide detail on sex stratif analysis for CAD . 
Please test additive models using logistic regression, accounting for genotype imputation uncertainty (i.e. SNP probability or dosage).

## Section 1: Analysis pipeline
Primary model here will be:
a. Males: CAD(all CAD cases vs. all non-CAD controls) = $SNP + Study covariates 
b. Females: CAD (all CAD cases vs. all non-CAD controls) = $SNP + Study covariates 
c. The should include the following filtering criteria: SNPs with imputation quality (INFO) > 0.3 and Mininum Allele Count (MAC)=10![image](https://github.com/TcheandjieuLab/CC4D_sex_stratified_analysis_plan/assets/32551968/c795e1f7-6cae-4fe5-acbe-466a5c6c2cb3)


### 1. Example script using PLINK
### 2. Example script using REGENIE
### 3. Example script using SAIGE

### list of variables to provide in the summary statistics for each sex specific GWAS
**SNP**:	SNP label for the variant	Identifier from the annotation file. Example: rs693 (chr2:7819; chr:pos:A1:A2)

**CHR**:	Chromosome on which SNP resides	Numeric for chromosomes 1-22. Example:	1

**POS with genome build**:	Position of SNP on chromosome	Base pairs on human genome build used. Example:	34000345

**EFFECT_ALLELE**:	Allele at this site to which the effect has been estimated	Capital letter (A,C,G,T). Example:	A

**NON_EFFECT_ALLELE**:	Other allele at this site (please check software documentation before label-ling a1,a2 as non-effect allele)	Capital letter (A,C,G,T). Example: G

**N_TOTAL**:	Total number of cases and controls analyzed	Numeric, integer	Example: 1243

**N_CASES**:	Total number of cases analyzed	(Numeric value)

**N_CONTROLS**:	Total number of controls analyzed	(Numeric value)

**EAF_ALL**:	Allele frequency of the EF-FECT_ALLELE in all	Frequency with 3 digits to the right of the decimal. Example: 0.354

**EAF_CASES**:	Allele frequency of the EF-FECT_ALLELE in cases analyzed	Frequency with 3 digits to the right of the decimal. Example: 0.354

**EAF_CONTROLS**:	Allele frequency of the EF-FECT_ALLELE in controls ana-lyzed	Frequency with 3 digits to the right of the decimal. Example: 0.354

**BETA**:	Estimate of the allelic effect, defined as the natural logarithm of the odds ratio, ln(OR)	Numeric float with 3 digits to the right of the decimal. Example: 0.203

**SE**:	Estimated standard error on the es-timate of the allelic effect, uncor-rected for genomic control	Numeric float with 4 digits to the right of the decimal. Example:	0.5611

**PVAL**:	Significance of the variant associa-tion, uncorrected for genomic con-trol	Scientific E notation with 3 digits to the right of the decimal.	Example: 3.24E-10

**INFO**:	Measure of information content for the imputed SNP result (range 0-1) (autosomal only)	Numeric float with 3 digits to the right of the decimal (set to missing if genotyped)	Example: 0.483

**HWE**: (should be used on the best guest genotypes, can be obtain from Plink independently)	P-value of the HWE test 	numerical	3.00E-06
![image](https://github.com/TcheandjieuLab/CC4D_sex_stratified_analysis_plan/assets/32551968/8f09b92d-9481-4f83-8a97-18d3047bc5fe)

## Section 2: sex interaction model (coming up soon. will be using the analytical pipeline as GLGC)

###list of variables to provide in the summary statistic 


***Note: For studies that have already conducted analyses, please discuss the models used (for example, some studies may have already adjusted for age and this is not considered a material deviation from the analysis plan). Use study appropriate software to account for (or exclude as appropriate) relatedness.***
***If undertaking new analyses to contribute then the following approach is advised. If providing previously generated results, please ensure the approach is documented and discuss any major deviations***

