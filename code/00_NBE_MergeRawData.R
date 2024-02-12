#### R script to transform counting data for Microcosm Exp 22####
#by Charlotte Kunze 28.03.22

## load packages
library(tidyverse)
library(readxl)
library(scales)
library(here)

#### Biovolume ####

##### import data #####
counting_ID_sampling1 <- read_excel("MicrocosmExp22/Data/counting_ID.xlsx", 
                                     sheet = "S0") %>%
  select(c(1:12)) %>%
  mutate(sampling = paste(1))
counting_ID_sampling3 <- read_excel("MicrocosmExp22/Data/counting_ID.xlsx", 
                                    sheet = "S3") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(3))
counting_ID_sampling6 <- read_excel("MicrocosmExp22/Data/counting_ID.xlsx", 
                                    sheet = "S6") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(6))
counting_ID_sampling9 <- read_excel("MicrocosmExp22/Data/counting_ID.xlsx", 
                                     sheet = "S9") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(9))
counting_ID_sampling12 <- read_excel("MicrocosmExp22/Data/counting_ID.xlsx", 
                                     sheet = "S12") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(12))
counting_ID_sampling15 <- read_excel("MicrocosmExp22/Data/counting_ID.xlsx", 
                                     sheet = "S15") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(15))

str(counting_ID_sampling15)

##### import meta data #####
ID_herma <- read_excel("MicrocosmExp22/Data/counting_ID_sampling9.xlsx",  sheet = "overview") %>%
  select(-plate, -V_ml, - magn, -sampling) 
ID_herma[ ID_herma == 'Asterio' ] <- 'A'
ID_herma[ ID_herma == 'DityCux' ] <- 'D'
ID_herma[ ID_herma == 'Guido' ] <- 'G'
ID_herma[ ID_herma == 'ThalaCux' ] <- 'T'
ID_herma[ ID_herma == 'Rhizo' ] <- 'R'


names(ID_herma)
names(counting_ID_sampling9)
names(counting_ID_sampling15)

counting_ID_sampling<- rbind(counting_ID_sampling1, counting_ID_sampling3,counting_ID_sampling6,counting_ID_sampling9, counting_ID_sampling12,counting_ID_sampling15)

counting_ID_sampling <-left_join(counting_ID_sampling, ID_herma, by = c('no', 'temp', 'nut', 'rep', 'species'))

##### calculate cells_ml #####
cells_mL_calculation_Utermo_hl_05chamber <- read_excel("MicrocosmExp22/Data/cells_mL_calculation_Utermöhl_05chamber.xlsx")
str(cells_mL_calculation_Utermo_hl_05chamber)
data <- left_join(counting_ID_sampling,cells_mL_calculation_Utermo_hl_05chamber, by = c('magn', 'V_ml'))

data <- data %>%
  mutate(cells_ml  = counts*At_mm/Afield*GF*V_ml) %>%
  mutate(sampling = as.numeric(sampling)) %>%
  filter(nut == 'Low')
names(data)


##### upload Cell Measurements for BioVolume calculation #####

BioVolume_exp22 <- read_excel("MicrocosmExp22/Data/BioVolume_exp22.xlsx",sheet = "Sheet1") %>%
  rename(speciesID = species) %>%
  select(-sampling,-comment, -Ansatz)
names(BioVolume_exp22)

BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Ditylum']<-'DityCux'
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Thala']<-'ThalaCux'

AllData <- data%>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = cells_ml * BV) %>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = ''), ifelse(species == 'duo', paste(A1, A2, sep = ''), 
                                                                            ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))))%>%
  select(no, speciesID,species, combination, temp, rep, sampling, cellVolume)
summary(AllData)

### create csv to work with ###
write.csv(AllData, file = 'MicrocosmExp22/Data/AllRawData_InclBV.csv')


#### Temperature data ####

## R Script to gather temperature data 

### required packages ###
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)
library(here)

setwd("~/Desktop/Exp22/Temperature_Charly_Planktotrons_2022")

## import excel files using a Loop
directory <-  paste(getwd(),sep = '') # hier sind meine Daten

#read file (from [here][1])
file.list<- list.files(full.names = TRUE, pattern = "*.xls")

df.list <- lapply(file.list, read_excel)
df <- bind_rows(df.list, .id = "id")
df <- filter(df, !unit %in% c(3, 4,9))
str(df)
write.csv(x = df, file = here('~/Desktop/Exp22/MicrocosmExp22/all_temperatures.csv'))
