#### Main Script ####

# Main script for creating plots and perform analysis 

## Before running: 
# 1. Check your working directory, if not done, please set to the main folder 
getwd()

# 2. load package here
library(here)

# 3. create an output folder in the MicrocosmExp22 folder, if you havent retrieved it from Github
dir.create(here('output')) # all output will be stored here



## Run to create and save plots ##

# Create Biomass Fig. 1 + Fig S2&S3
source(here("code/02_NBE_BiomassFigure_supplement.R"))


# Create NBES Fig 2 & 3 + Fig S1
source(here("code/03_NBES_calculation.R"))
# Note: Within the code above the following is run to create the NBES/NBE plots 
# source(here("code/04_NBE_HectorLoreau_NetBiodivEffect.R"))

# calculate influence on grand mean
source(here('code/06_NBES_slopes.R'))

# Create temperature curves for supplements
source(here('code/01_NBE_temperature_planktotrons.R'))

### Run Stats ####

# Please run separately by actively open the R script
source(here("code/05_NBE_Statistics_Contrasts.R"))



