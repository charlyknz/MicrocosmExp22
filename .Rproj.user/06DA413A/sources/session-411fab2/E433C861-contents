# README NBES 
R script to quantify the net biodiversity effect on ecological stability in a phytoplankton microcosm experiment

## To run the code..
Please store the contents of the repository on your device & run Main.R


## R Session info
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.1.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8





## Rcode:

### Main.R – creates all plots displayed in the MS, script accessess several other R scripts

### 00_NBE_MergeRawData: merges raw data, i.e. count data for each sampling with biovolume measurements

### 01_NBE_temperature_planktotrons.R: creates temperature curves for each treatment

### 02_NBE_BiomassFigure_supplement.R: Total and species-specific biomass for mono- and multi-species cultures

### 03_NBES_calculation.R: Analysis of net biodiversity effect on stability (NBES) using the Overall Ecological Vulnerability metric (OEV, Urrutia-Cordero et al. 2021), resistance and temporal stability, measured as Coefficient of Variation (CV)

### 04_NBE_HectorLoreau_NetBiodivEffect.R: calculation of the net biodiversity effect on ecosystem functioning after Loreau and Hector (2001). 

### 03_NBE_Statistics_Contrasts.R: Statistics introduced to analyse the influence of temperature and species composition on the NBES as well as the net biodiversity effect on functioning. 


## Files:
### AllRawData_InclBV.csv: Created file containing merged raw data including calculated cell Volume from Biovolume of species. 

no: unique ID of N x temp combination 
 
speciesID: species-specific information
 
 combination: species combinations, see MS for abbreviations
 
 temp: temperature treatment, i.e. fluctuation, increase, increase and fluctuation 
 
  rep: replicate no 
  
  sampling: sampling number - 1,3,6,9,12,15. We sampled every 4th day.
  

### NBES.csv:  Data on Net Biodiversity effect on stability

combination: species combinations, see MS for abbreviations

rep: replicate number (1-3)

NBE: Net biodiversity effect 

AUC.RR_exp: expected Stability 

AUC.RR_obs: observed Stability

N: Species richness level (ranging from 2-5)

temp: temperature treatment, i.e. fluctuation, increase, increase and fluctuation 



  

### NBEonFunctioning.csv: Data on Net Biodiversity effect on functioning

combination: species combinations, see MS for abbreviations.

rep: replicate number (1-3)

NetEffect: Net biodiversity effect on functioning

N: Species richness level (ranging from 2-5)

temp: temperature treatment, i.e. control, fluctuation, increase, increase and fluctuation 
  
