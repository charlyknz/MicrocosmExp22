#### R script to analyse counting data for Microcosm Exp 22####
#by Charlotte Kunze 28.03.22

## load packages
library(tidyverse)
library(readxl)
library(scales)
library(here)
#### import data ####
counting_ID_sampling1 <- read_excel("counting_ID.xlsx", 
                                     sheet = "S0") %>%
  select(c(1:12)) %>%
  mutate(sampling = paste(1))
counting_ID_sampling3 <- read_excel("counting_ID.xlsx", 
                                    sheet = "S3") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(3))
counting_ID_sampling6 <- read_excel("counting_ID.xlsx", 
                                    sheet = "S6") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(6))
counting_ID_sampling9 <- read_excel("counting_ID.xlsx", 
                                     sheet = "S9") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(9))
counting_ID_sampling12 <- read_excel("counting_ID.xlsx", 
                                     sheet = "S12") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(12))
counting_ID_sampling15 <- read_excel("counting_ID.xlsx", 
                                     sheet = "S15") %>%
  select(c(1:12))%>%
  mutate(sampling = paste(15))

str(counting_ID_sampling15)

#### import meta data ####
ID_herma <- read_excel("counting_ID_sampling9.xlsx",  sheet = "overview") %>%
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

#### calculate cells_ml ####
cells_mL_calculation_Utermo_hl_05chamber <- read_excel("cells_mL_calculation_Utermöhl_05chamber.xlsx")
str(cells_mL_calculation_Utermo_hl_05chamber)
data <- left_join(counting_ID_sampling,cells_mL_calculation_Utermo_hl_05chamber, by = c('magn', 'V_ml'))

data <- data %>%
  mutate(cells_ml  = counts*At_mm/Afield*GF*V_ml) %>%
  mutate(sampling = as.numeric(sampling)) %>%
  filter(nut == 'Low')
names(data)
dum <- filter(data, is.na(cells_ml))
write.csv(data, file = 'Exp22_wrangledData.csv')

