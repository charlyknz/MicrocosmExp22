# README  
R script to quantify the net biodiversity effect on ecological stability in a phytoplankton microcosm experiment. Data are stored on Figshare: 10.6084/m9.figshare.25568490


## To run the code..
Please store the contents of the repository on your device, open the Rproject & run Main.R. Please follow the instructions.


## R Session info
R version R version 4.4.3
Running under: macOS Sequoia 15.3.2

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

### 05_NBE_Statistics_Contrasts.R: Statistics introduced to analyse the influence of temperature and species composition on the NBES as well as the net biodiversity effect on functioning. 

### 06_NBES_slopes: calculates influence of single species on NBES.

