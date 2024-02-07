#### NBE ####
library(tidyverse)
library(ggpubr)

#### import data ####
d3 <- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/netEffect.csv')%>%
  select(combination, temp, rep, N, speciesID, NBE)
str(d3)

data <- read.csv('NBEonFunctioning.csv') %>%
  select(combination, temp, rep, N, speciesID, NetEffect) %>%
  rename(NBE_Func = NetEffect)
str(data)

allNBE <- data %>%
  left_join(., d3, by = c('combination', 'temp', 'rep','N','speciesID'))

#### plot ####
allNBE %>%
  filter(temp != 'Constant')%>%
ggplot(., aes(x=NBE_Func, y = NBE, color = temp ))+
  geom_point()+
  facet_wrap(~N, scales = 'free_x')+
  theme_bw()
ggsave(plot = last_plot(), file = '~/Desktop/Exp22/MicrocosmExp22/output/FacetNBEs.png', width = 10, height = 3)

allNBE %>%
  filter(temp != 'Constant')%>%
ggscatter(., x = 'NBE_Func', y = 'NBE', add = 'reg.line', color = 'temp')+
  stat_cor(aes(color = temp))
ggsave(plot = last_plot(), file = '~/Desktop/Exp22/MicrocosmExp22/output/correlationNBEs_temp.png', width = 4, height = 3)
