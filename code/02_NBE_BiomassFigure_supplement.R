#### R script to prepare data for analysis with Dome ####
# 02.10.2022
# by Charlotte Kunze

# packages
library(tidyverse)
library(readxl)
library(here)
library(cowplot)


#### create color and shape palettes ####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,'#c7b514' ) # temp treatments
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7") #species
shape_values <- c(23,16, 17, 15)


#### import species-specific biomass ####
rawData <- read.csv('Data/AllRawData_InclBV.csv')
summary(rawData)


# calculate mean biomass and transform biovolume in um3 to mm3 per ml
biomass<- rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by( species,combination, speciesID, temp, day) %>%
  mutate(cellV_mm_ml = cellVolume/10^9)%>%
  summarise(meanV = mean(cellV_mm_ml, na.rm = T),
            sdV = sd(cellV_mm_ml, na.rm = T),
            seV = sdV/sqrt(n()))

# change labels for temperature tratments and species
biomass$temp[biomass$temp=='CS'] <- 'Constant'
biomass$temp[biomass$temp=='fluct'] <- 'Fluctuation'
biomass$temp[biomass$temp=='inc'] <- 'Increase'
biomass$temp[biomass$temp=='inc+fluc'] <- 'Increase + Fluctuation'

biomass$speciesID[biomass$speciesID=='Asterio'] <- 'Asterionellopsis'
biomass$speciesID[biomass$speciesID=='DityCux'] <- 'Ditylum'
biomass$speciesID[biomass$speciesID=='Guido'] <- 'Guinardia'
biomass$speciesID[biomass$speciesID=='Rhizo'] <- 'Rhizosolenia'
biomass$speciesID[biomass$speciesID=='ThalaCux'] <- 'Thalassionema'


#### Plots with species-specific responses ####

#Two species 
Duoplot <- biomass %>%
  filter(species == 'duo') %>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 1, linetype = 'dashed')+
  geom_point(alpha = 1, size = 1.5)+
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
  theme(legend.position = 'bottom',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
Duoplot
ggsave(plot = Duoplot, file = here('output/ExtendedData_FigureS3_Duo_Biomass.pdf'), width = 12, height = 12)

#4 species
Quattroplot <- biomass %>%
  filter(species == 'quattro') %>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 1, linetype = 'dashed')+
  geom_point(alpha = 1, size = 1.5)+
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

# 5 species
mixplot <- biomass %>%
  filter(!species%in%c('duo','mono', 'quattro') )%>%
  ggplot(., aes( x = day, y = meanV, color = speciesID, group = speciesID))+
  geom_line(alpha = 1, linetype = 'dashed')+
  geom_point(alpha = 1, , size = 1.5)+
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
  theme(legend.position = 'bottom',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
mixplot

plot_grid(Quattroplot, mixplot, ncol = 1, rel_heights = c(7/10, 3/10),labels = c('(a)', '(b)'))
ggsave(plot = last_plot(), file = here('output/ExtendedData_FigureS4_QuattroMixBiomass.pdf'), width = 12, height = 12.5)

# monocultures
Mono <- biomass %>%
  filter(species == 'mono') %>%
  ggplot(., aes( x = day, y = meanV, color = temp, group = temp, shape = temp))+
  geom_line(alpha = 1, linetype = 'dashed')+
  geom_point(alpha = 1, size = 1.5)+
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
  theme(legend.position = 'bottomright',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
Mono

#### Total biomass ####

biomass1<- rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by(species, combination, temp, day) %>%
  mutate(cellV_mm_ml = cellVolume/10^9)%>%
  summarise(meanV = sum(cellV_mm_ml, na.rm = T),
            sdV = sd(cellV_mm_ml, na.rm = T),
            seV = sdV/sqrt(n()))


# rename labels
biomass1$temp[biomass1$temp=='CS'] <- 'Constant'
biomass1$temp[biomass1$temp=='fluct'] <- 'Fluctuation'
biomass1$temp[biomass1$temp=='inc'] <- 'Increase'
biomass1$temp[biomass1$temp=='inc+fluc'] <- 'Increase + Fluctuation'


### Create Figure on total biomass over time ###
biomass1$combination <- factor(as.factor(biomass1$combination) , 
                               levels = c('A','D','G','R','T',
                                          "AD" ,"AG" ,"AR","AT" ,"DG","DR","DT", "GR" , "GT" ,"RT",
                                          "ADGR" ,"ADGT", "ADRT" ,"AGRT", "DGRT" ,"ADGRT"))

biomass1 %>%
  ggplot(., aes( x = day, y = meanV, color = temp, group = temp, shape = temp))+
  geom_line(alpha = 1, linetype = 'dashed')+
  geom_point(alpha = 1, size = 1.5)+
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
  theme(legend.position = 'bottom',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))

ggsave(plot = last_plot(), file = here('output/ExtendedData_FigureS2_Biomass.pdf'), width = 15, height = 15)



