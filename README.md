![image](https://github.com/user-attachments/assets/1e85c53a-f530-4236-b162-93a264a5dc5a)


# R functions to analysis the Whole Genome Sequencing results of ACORN-WGS 

## About ACORN
ACORN is a Wellcome funded human health clinical AMR surveillance project led by [the Mahidol-Oxford Tropical Medicine Research Unit (MORU)](https://www.tropmedres.ac/) and [the Oxford University of Oxford Clinical Research Unit (OUCRU)](https://www.oucru.org/).

**Why is ACORN needed?**
Existing AMR surveillance systems are based mostly on diagnostic microbiology laboratory antimicrobial susceptibility testing results alone, which limits interpretability of resistant proportions. Resulting data fail to give relevant feedback for treatment decisions for local clinicians and do not allow for direct assessment and subsequent modelling of the clinically relevant impacts and burden of drug resistant infections (DRI). Tools to capture and analyse AMR data in low- and middle-income countries (LMIC) are scarce, which hinders engagement with and use of available data.

To fill these gaps, the major aim of ACORN is to develop and test a comprehensive data capture system for patient-focussed AMR surveillance in LMIC settings. Surveillance will include diagnostic stewardship activities. Data collected will harmonise with and expand on the pathogen-focussed WHO Global Antimicrobial Resistance Surveillance System to enable accurate classification of infection syndromes and patient outcomes. These data will be of critical importance to estimate syndromic and/or pathogen outcomes and associated costs: i.e. how many people die from DRIs and how much does AMR cost?

For more information of ACORN, please visit [the website](https://acornamr.net/#/).

## Overview 
These functions are written to transform the data and visualise the result of ACORN-WGS. 
* `plot_heatmap_vf`
* `plot_patient_timeline`
* `plot_corr_pheno_genotype`
* `amr_transform`

## Dependencies R packages 
*  `dplyr`
*  `ggplot`
*  `stringr`
*  `ComplexHeatmap`
*  `tidyr`

## Author
[Tung Trinh- Oxford University Clinical Research Unit](https://www.oucru.org/people/trinh-son-tung/)
