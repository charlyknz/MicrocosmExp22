#### Statistic Rscript ####
# by Charlotte Kunze

#packages
library(tidyverse)
library(ggpubr)

#### Import data ####

##### NBES data ####
# run R scripts 
#source(here("~/Desktop/Exp22/MicrocosmExp22/code/03_NBE_calculatonNBE.R"))
#str(d3)

# or import as csv
d3 <- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/netEffect.csv')%>%
  select(combination, temp, rep, N, speciesID, NBE)
str(d3)

##### NBE on Funct data ####
#source(here("~/Desktop/Exp22/MicrocosmExp22/code/05_NBE_HectorLoreau_NetBiodivEffect.R"))
#str(data)

#or

# import as csv
data <- read.csv('NBEonFunctioning.csv') %>%
  select(combination, temp, rep, N, NetEffect) %>%
  rename(NBE_Func = NetEffect)
str(data)

allNBE <- data %>%
  left_join(., d3, by = c('combination', 'temp', 'rep','N'))

#### plot ####
allNBE %>%
  filter(temp != 'Constant')%>%
ggplot(., aes(x=NBE_Func, y = NBE, color = temp ))+
  geom_point()+
  facet_wrap(~N, scales = 'free_x')+
  theme_bw()

allNBE %>%
  filter(temp != 'Constant')%>%
ggscatter(., x = 'NBE_Func', y = 'NBE', add = 'reg.line')+
  stat_cor()
ggsave(plot = last_plot(), file = '~/Desktop/Exp22/MicrocosmExp22/output/correlationNBEs_temp.png', width = 4, height = 3)
