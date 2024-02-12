#### R script to prepare data for analysis with Dome ####
# 02.10.2022
# by Charlotte Kunze

# packages
library(tidyverse)
library(readxl)
library(here)
library(cowplot)

#setwd("~/Desktop/Exp22/MicrocosmExp22")

#### create color palettes ####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" ) # temp treatments
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7") #species
shape_values <- c(23,16, 17, 15)


#### import species-specific biomass ####
rawData <- read.csv('Data/AllRawData_InclBV.csv')
summary(rawData)


biomass<- rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by( species,combination, speciesID, temp, day) %>%
  mutate(cellV_mm_ml = cellVolume/10^9)%>%
  summarise(meanV = mean(cellV_mm_ml, na.rm = T),
            sdV = sd(cellV_mm_ml, na.rm = T),
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


#### Indivdiual species plots ####

Duoplot <- biomass %>%
  filter(species == 'duo') %>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~'['~mm^3~'/'~ml~']'), x = 'Time [days]', color = 'Species')+
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
ggsave(plot = Duoplot, file = here('output/FigS3_Duo_Biomass.png'), width = 12, height = 12)

Quattroplot <- biomass %>%
  filter(species == 'quattro') %>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~'['~mm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature')+
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


mixplot <- biomass %>%
  filter(!species%in%c('duo','mono', 'quattro') )%>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~'['~mm^3~'/'~ml~']'), x = 'Time [days]', color = 'Species')+
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

plot_grid(Quattroplot, mixplot, ncol = 1, rel_heights = c(6.5/9, 2.5/9),labels = c('(a)', '(b)'))
ggsave(plot = last_plot(), file = here('output/FigS4_QuattroMixBiomass.png'), width = 12, height = 12.5)


Mono <- biomass %>%
  filter(species == 'mono') %>%
  ggplot(., aes( x = day, y = meanV, color = temp, group = temp, shape = temp))+
  geom_line(alpha = 0.8, linetype = 'dashed')+
  geom_point(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
  labs(y =  expression(Total~Biovolume~'['~mm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature', shape =  'Temperature')+
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

#### sum of biomass ####
biomass1<- rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by(species, combination, temp, day) %>%
  mutate(cellV_mm_ml = cellVolume/10^9)%>%
  summarise(meanV = sum(cellV_mm_ml, na.rm = T),
            sdV = sd(cellV_mm_ml, na.rm = T),
            seV = sdV/sqrt(n()))


#brewer.pal(n = 8, name = "Set1")
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
  labs(y =  expression(Total~Biovolume~'['~mm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature', shape =  'Temperature')+
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

ggsave(plot = last_plot(), file = here('output/Fig1Biomass.png'), width = 15, height = 15)



