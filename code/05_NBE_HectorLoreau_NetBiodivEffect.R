#### R script to calculate Hector and Loreau's NBE ####
# by Charlotte Kunze

# packages
library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(RColorBrewer)

setwd("~/Desktop/Exp22/MicrocosmExp22/Data")

#### import data ####
expData<-read.csv('AllRawData_InclBV.csv') %>%
  select(-X)
names(expData)

#### colors and shapes ####
temp1Palette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )
shape_values <- c(23,16, 17, 15)

#### Net Biodiversity effect on Functioning (Hector & Loreau, 2001) ####

# Formula: NetEffect = sum(Yoi, na.rm = T)- sum(relBV*Mi)
# with Mi   : Biomass of species at tmax in Monoculture
#      Yoi  : observed biomass of species i in mixture at tmax 
#      relBV: relative Biovolume/ Proportion of species in mixture at t0



# relBV 
# relative biomass of species in mixture at t0
RelBV_t0_ <- expData %>%
  filter(sampling == 1) %>%
  group_by(combination, temp, rep) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(combination, temp, speciesID, rep, relBV)
names(RelBV_t0_)


# Yoi
# observed biomass of species in mixture and monoculture  
counts_duo_tmax_ <- expData %>%
  filter(sampling == 15) %>%
  filter(species == 'duo') %>%
  rename(Yoi = cellVolume)%>%
  select(combination, temp,rep, speciesID, Yoi)
names(counts_duo_tmax_)

# Mi
# biomass of species i at tmax in monoculture
counts_mono_tmax_ <- expData %>%
  filter(sampling == 15) %>%
  filter(species == 'mono') %>%
  rename(Mi = cellVolume) %>%
  select( temp, speciesID, rep,Mi)
names(counts_mono_tmax_)
# NBE
#combine df RelAbund_t0 and OD_counts_tmax
NBE_duo <- left_join(RelBV_t0_,counts_mono_tmax_, by = c( 'temp', 'speciesID', 'rep')) %>%
  left_join(., counts_duo_tmax_, by = c( 'temp', 'speciesID', 'combination', 'rep')) %>%
  drop_na(Yoi) %>%
  drop_na(Mi) %>%
  group_by(temp, combination, rep) %>%
  mutate(NetEffect = sum(Yoi, na.rm = T)- sum(relBV*Mi)) %>%
  group_by(temp, combination) %>%
  mutate(N = as.numeric(paste(2)),
        # compEffect = N*mean(NetEffect)*mean(Mi),
         #selecEffect = N*cov(NetEffect, Mi)
         )

#### Plot 2 species ####
temp1Palette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )

unique(NBE_duo$temp)

NBE_duo$temp[NBE_duo$temp=='CS'] <- 'Constant'
NBE_duo$temp[NBE_duo$temp=='fluct'] <- 'Fluctuation'
NBE_duo$temp[NBE_duo$temp=='inc'] <- 'Increase'
NBE_duo$temp[NBE_duo$temp=='inc+fluc'] <- 'IncreaseFluctuation'


# Net effect plot
NBE_duo %>%
  group_by(combination,  temp) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n())) %>%
  ggplot(., aes(x = combination, y = mean.effect, color = temp, shape = temp))+
  geom_hline(yintercept = 0, color ='darkgrey')+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect + se.effect), width = .2, alpha = .7 )+
  geom_point(size = 2)+
  scale_color_manual(values= temp1Palette)+
  scale_shape_manual(values= shape_values)+
  labs(y = 'Net Biodiversity Effect', title = 'Two species in combination', color = 'Treatment', shape = 'Treatment')+
  #facet_grid(~temp, scales = 'free_y')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'bottom')

ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NBEonFunctioning_onefacet.png'), width = 8, height = 4)

#### Netto Biodiversity Effect Mix ####

# Yoi
# observed biomass of species in mixture and monoculture  
counts_tmax <- expData %>%
  filter(sampling == 15) %>%
  filter(species%in% c('quattro', 'MIX')) %>%
  rename(Yoi = cellVolume)%>%
  select(combination, temp,rep, speciesID, Yoi)
names(counts_tmax)

# Mi
# biomass of species i at tmax in monoculture
names(counts_mono_tmax_)

# NBE
#combine df RelAbund_t0 and OD_counts_tmax
NBE_mix <- left_join(RelBV_t0_,counts_mono_tmax_, by = c( 'temp', 'speciesID', 'rep')) %>%
  right_join(., counts_tmax, by = c( 'temp', 'speciesID', 'combination', 'rep')) %>%
  drop_na(Yoi) %>%
  drop_na(Mi) %>%
  group_by(temp, combination, rep) %>%
  mutate(NetEffect = sum(Yoi, na.rm = T)- sum(relBV*Mi)) %>%
  group_by(temp, combination) %>%
  mutate(N = as.numeric(paste(2)),
         # compEffect = N*mean(NetEffect)*mean(Mi),
         #selecEffect = N*cov(NetEffect, Mi)
  )

#### Plot 4-5 species ####
unique(NBE_mix$temp)

NBE_mix$temp[NBE_mix$temp=='CS'] <- 'Constant'
NBE_mix$temp[NBE_mix$temp=='fluct'] <- 'Fluctuation'
NBE_mix$temp[NBE_mix$temp=='inc'] <- 'Increase'
NBE_mix$temp[NBE_mix$temp=='inc+fluc'] <- 'IncreaseFluctuation'

NBE_mix$combination[NBE_mix$combination == 'ADGRT'] <- '5spec'

# Net effect plot
NBE_mix %>%
  group_by(combination,  temp) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n())) %>%
  ggplot(., aes(x = combination, y = mean.effect, color = temp, shape = temp))+
  geom_hline(yintercept = 0, color ='darkgrey')+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect + se.effect), width = .2,alpha = .7 )+
  geom_point(size = 2)+
  scale_color_manual(values= temp1Palette)+
  scale_shape_manual(values= shape_values)+
  labs(y = 'Net Biodiversity Effect', title = '4 and 5 species in combination', color = 'Treatment', shape = 'Treatment')+
 # facet_grid(~temp, scales = 'free_y')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'bottom')

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NBEonFunctioning_45_onefacet.png'), width = 8, height = 4)

#### Overall NBE ####

allNetBiodiv <- bind_rows(NBE_duo, NBE_mix)%>%
  mutate(N = str_length(combination)) 
#write.csv(allNetBiodiv, file =  'NBEonFunctioning.csv')

NBE<-allNetBiodiv%>%
  group_by( temp, N) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n()))%>%
  ggplot(., aes(x = N, y = mean.effect,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect +se.effect), width = .3, alpha = .5 )+
  geom_point(size = 2.5, alpha = 0.7)+
  labs(x = 'Species Richness', y = 'Net Biodiversity Effect on Functioning', color = 'Treatment', shape = 'Treatment')+
   scale_x_continuous(limits = c(1.5,5.5), breaks = c(2,4,5))+
  scale_color_manual(values= temp1Palette)+
  scale_shape_manual(values= shape_values)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
NBE

### with raw data ###

NBE1 <- allNetBiodiv%>%
  group_by( temp, N) %>%
  mutate(mean.effect = mean(NetEffect, na.rm = T),
         sd.effect = sd(NetEffect,na.rm = T),
         se.effect = sd.effect/sqrt(n()))%>%
  ggplot(. )+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(x = N, y = NetEffect,color = temp), size = 1, shape = 1, alpha = 0.2)+
  geom_errorbar(aes(x = N, y = mean.effect,ymin = mean.effect - se.effect, ymax = mean.effect +se.effect,color = temp), width = .3, alpha = .5 )+
  geom_point(aes(x = N, y = mean.effect,color = temp, shape = temp), size = 2.5, alpha = 0.7)+
  labs(x = 'Species Richness', y = 'Net Biodiversity effect', color = 'Treatment', shape = 'Treatment')+
  scale_x_continuous(limits = c(1,5.5), breaks = c(2,4,5))+
  scale_color_manual(values= temp1Palette)+
  scale_shape_manual(values= shape_values)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
NBE1

####

allNetBiodiv$N[allNetBiodiv$N == 2] <- '2 species'
allNetBiodiv$N[allNetBiodiv$N == 4] <- '4 species'
allNetBiodiv$N[allNetBiodiv$N == 5] <- '5 species'

allNetBiodiv$combination[allNetBiodiv$combination == '5spec'] <- 'ADGRT'

allNetBiodiv$N <- factor(as.factor(allNetBiodiv$N), levels= c( '2 species','4 species','5 species'))
plot245<- allNetBiodiv %>%
  group_by( combination, temp, N) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n()))%>%
  ggplot(., aes(x = combination, y = mean.effect,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect +se.effect), width = .3, alpha = .5 )+
  geom_point(size = 2.5, alpha = 0.7)+
  labs(x = 'Species Combinations', y = 'Net Biodiversity Effect on Functioning', color = 'Treatment', shape = 'Treatment')+
  scale_color_manual(values= temp1Palette)+
  scale_shape_manual(values= shape_values)+
  facet_wrap(~N, scales = 'free_x')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12),
        legend.position = 'none')+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
plot245

library(cowplot)
blank <- ggplot()+
  geom_blank()+
  theme_void() 


plot_grid(NBE+theme(legend.position = 'none'), plot245, ncol = 2, labels = c('(c)', '(d)'), rel_widths= c(1/3, 2/3), hjust = -1.5)
plot_grid( NBE+theme(legend.position = 'none'),plot245,blank,blank, hjust = -0.05, labels = c('(c)', '(d)'), ncol = 2,rel_widths = c( 1/3,2/3), rel_heights = c(1,0.3))
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NBEonFunctioning_all.png'), width = 13, height = 6)




