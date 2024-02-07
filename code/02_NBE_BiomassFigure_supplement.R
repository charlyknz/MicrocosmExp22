#### R script to prepare data for analysis with Dome ####
# 02.10.2022
# by Charlotte Kunze

# packages
library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(RColorBrewer)
library(cowplot)

setwd("~/Desktop/Exp22/MicrocosmExp22/Data")

# create color palettes
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" ) # temp treatments
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7") #species
shape_values <- c(23,16, 17, 15)

#### import species-specific biomass ####
rawData <- read.csv('AllRawData_InclBV.csv')
summary(rawData)


biomass<- rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by( species,combination, speciesID, temp, day) %>%
  mutate(transV = cellVolume/10^9)%>%
  summarise(meanV = mean(transV, na.rm = T),
            sdV = sd(transV, na.rm = T),
            seV = sdV/sqrt(n()))

biomass$temp[biomass$temp=='CS'] <- 'Constant'
biomass$temp[biomass$temp=='fluct'] <- 'Fluctuation'
biomass$temp[biomass$temp=='inc'] <- 'Increase'
biomass$temp[biomass$temp=='inc+fluc'] <- 'IncreaseFluctuation'

biomass$speciesID[biomass$speciesID=='Asterio'] <- 'Asterionellopsis'
biomass$speciesID[biomass$speciesID=='DityCux'] <- 'Ditylum'
biomass$speciesID[biomass$speciesID=='Guido'] <- 'Guinardia'
biomass$speciesID[biomass$speciesID=='Rhizo'] <- 'Rhizosolenia'
biomass$speciesID[biomass$speciesID=='ThalaCux'] <- 'Thalassionema'


###+ single plots 

Duoplot <- biomass %>%
  filter(species == 'duo') %>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Time [days]', color = 'Species')+
  facet_grid(~combination~temp, scales = 'free')+
  scale_colour_manual(values = cbbPalette)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'bold'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'bottom')
Duoplot
ggsave(plot = Duoplot, file = here('MicrocosmExp22/output/DuoBiomass.png'), width = 12, height = 12)

Quattroplot <- biomass %>%
  filter(species == 'quattro') %>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature')+
  facet_grid(~combination~temp, scales = 'free')+
  scale_colour_manual(values = cbbPalette)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 12,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 12, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'bold'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'none')
Quattroplot

#ggsave(plot = Quattroplot, file = here('MicrocosmExp22/output/QuattroBiomass.png'), width = 12, height = 8)

mixplot <- biomass %>%
  filter(!species%in%c('duo','mono', 'quattro') )%>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Time [days]', color = 'Species')+
  facet_grid(~combination~temp, scales = 'free_y')+
  scale_colour_manual(values = cbbPalette)+
 # scale_x_continuous(limits = c(0,16), breaks = c(1,3,6,9,12,15))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 12,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 12, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'bold'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'bottom')
mixplot

#ggsave(plot = mixplot, file = here('MicrocosmExp22/output/MixBiomass.png'), width = 10, height = 3)

plot_grid(Quattroplot, mixplot, ncol = 1, rel_heights = c(6.5/9, 2.5/9),labels = c('(a)', '(b)'))
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/QuattroMixBiomass.png'), width = 12, height = 12)


Mono <- biomass %>%
  filter(species == 'mono') %>%
  ggplot(., aes( x = day, y = meanV, color = temp, group = temp, shape = temp))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~BioV~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature', shape =  'Temperature')+
  facet_wrap(~speciesID,nrow = 1, scales = 'free_y')+
  scale_shape_manual(values = shape_values)+
  scale_colour_manual(values = tempPalette)+
  # scale_x_continuous(limits = c(0,30), breaks = c(1,3,6,9,12,15))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'italic'),
        legend.text = element_text(size = 13), legend.title = element_text(size = 14,  face = 'bold'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme(legend.position = 'bottomright')
Mono
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/MonocultureBiomass.png'), width = 10, height = 3)

#### sum of biomass ####
biomass1<- rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by(species, combination, temp, day) %>%
  mutate(transV = cellVolume/10^9)%>%
  summarise(meanV = sum(transV, na.rm = T),
            sdV = sd(transV, na.rm = T),
            seV = sdV/sqrt(n()))


brewer.pal(n = 8, name = "Set1")
biomass1$temp[biomass1$temp=='CS'] <- 'Constant'
biomass1$temp[biomass1$temp=='fluct'] <- 'Fluctuation'
biomass1$temp[biomass1$temp=='inc'] <- 'Increase'
biomass1$temp[biomass1$temp=='inc+fluc'] <- 'IncreaseFluctuation'


### all ###
biomass1$combination <- factor(as.factor(biomass1$combination) , 
                               levels = c('A','D','G','R','T',
                                          "AD" ,"AG" ,"AR","AT" ,"DG","DR","DT", "GR" , "GT" ,"RT",
                                          "ADGR" ,"ADGT", "ADRT" ,"AGRT", "DGRT" ,"ADGRT"))

biomass1 %>%
  ggplot(., aes( x = day, y = meanV, color = temp, group = temp, shape = temp))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature', shape =  'Temperature')+
  facet_wrap(~combination, ncol = 5, scales = 'free')+
  scale_colour_manual(values = tempPalette)+
  scale_shape_manual(values = shape_values)+
  # scale_x_continuous(limits = c(0,16), breaks = c(1,3,6,9,12,15))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'bold'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  #theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme(legend.position = 'bottom')

ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/Fig2Biomass.png'), width = 15, height = 15)



### CV ####
names(biomass)
rawData$species <- factor(as.factor(rawData$species),levels=c("mono","duo",   "quattro" ,  "MIX"))
rawData$N <- NA
rawData$N[rawData$species == 'mono']<-'1'
rawData$N[rawData$species == 'duo']<-'2'
rawData$N[rawData$species == 'quattro']<-'4'
rawData$N[rawData$species == 'MIX']<-'5'

rawData$temp[rawData$temp=='fluct'] <- 'Fluctuation'
rawData$temp[rawData$temp=='inc'] <- 'Increase'
rawData$temp[rawData$temp=='inc+fluc'] <- 'IncreaseFluctuation'
rawData$temp[rawData$temp=='CS'] <- 'Control'


rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
   group_by(species,temp, combination,N, rep,sampling)%>%
   summarise(sum = sum(cellVolume)) %>%
   group_by(species, temp, combination,rep,N) %>%
   summarise(mean = mean(sum),
         sd = sd(sum),
         CV = mean/sd)%>%
ggplot(.,aes(x = N, y = CV, color = temp))+
  geom_point( alpha = 0.5,size=2)+
  facet_grid(~temp)+
 geom_boxplot(alpha = 0.5)+
  labs(y = 'Temporal stability (Mean/SD)')+
  scale_colour_manual(values = tempPalette)+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/CV.png'), width = 8, height = 3)
