#### Main Script ####

# Main script for creating plots and perform analysis 

## Before running: 
# 1. Set your working directory to the main folder 
setwd("~/Desktop/Exp22/MicrocosmExp22")

# 2. create an output folder in the MicrocosmExp22 folder 
dir.create(here('output'))



## Create plots ##

# Create biomass Fig. 1
source(here("~/Desktop/Exp22/MicrocosmExp22/code/02_NBE_BiomassFigure_supplement.R"))


# Create Figures 2 & 3 
source(here("~/Desktop/Exp22/MicrocosmExp22/code/03_NBES_calculation.R"))
# Note: Within the code above the following is run to create the NBES/NBE plots 
# source(here("~/Desktop/Exp22/MicrocosmExp22/code/05_NBE_HectorLoreau_NetBiodivEffect.R"))


# Create temperature curves for supplements
source(here('~/Desktop/Exp22/MicrocosmExp22/code/01_NBE_temperature_planktotrons.R'))

### Run Stats ####
source(here("~/Desktop/Exp22/MicrocosmExp22/code/04_NBE_Statistics_Contrasts.R"))



