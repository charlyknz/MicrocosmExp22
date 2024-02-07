### R script to calculate NBE on Stability ###
# 11.07.2022
# by Charlotte Kunze

#### packages ####
library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(RColorBrewer)
library(cowplot)
library(patchwork)

setwd("~/Desktop/Exp22/MicrocosmExp22/Data")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
#### HI nutrients ####

#calculate mean cells_ml in control dataset 

countData_HI<- expData%>%
  filter(nut == 'HI' & sampling %in% c(1,9,15)) %>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = ''), ifelse(species == 'duo', paste(A1, A2, sep = ''), 
                                                                            ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))))%>%
  group_by(temp, nut, species, Alg, speciesID, sampling, A1, A2, A3, A4, A5, combination) %>%
  summarise(mean.cells = mean(cells_ml, na.rm = T),
            sd.cells = sd(cells_ml, na.rm = T),
            se.cells= sd.cells/sqrt(n())) %>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5) 
names(countData_HI)

#### abundance Net biodiviserty effect ####
# formula:  ∑Yoi -  ∑relAb_t0 * Mi
# Mi - bioimass of species i at tmax in mono
# Yoi - observed biomass of species in mix at tmax
# relAb - relative biomass of species in mix at t0


# only mono and duo cultures
OD_counts_HI <- filter(countData_HI, species %in% c('mono', 'duo')) 

# relative biomass of species in mxiture at t0
RelAbund_t0_HI <- OD_counts_HI %>%
  filter(sampling == 1) %>%
  filter(species == 'duo') %>%
  group_by(nut, combination, temp) %>%
  mutate(sum = sum(mean.cells)) %>%
  ungroup() %>%
  mutate(relAb = mean.cells/sum) %>%
  select(nut, combination, temp, speciesID, relAb)

#observed biomass of species in mixture and monoculture  
counts_duo_tmax_HI <- OD_counts_HI %>%
  filter(sampling == 15) %>%
  filter(species == 'duo') %>%
  rename(Yoi = mean.cells)%>%
  select(nut, combination, temp, speciesID, Yoi)

counts_mono_tmax_HI <- OD_counts_HI %>%
  filter(sampling == 15) %>%
  filter(species == 'mono') %>%
  rename(Mi = mean.cells) %>%
  select(nut, temp, speciesID, Mi)

#combine df RelAbund_t0 and OD_counts_tmax
NetBiodiversity_HI <- left_join(RelAbund_t0_HI,counts_mono_tmax_HI, by = c('nut',  'temp', 'speciesID')) %>%
  left_join(., counts_duo_tmax_HI, by = c('nut',  'temp', 'speciesID', 'combination')) %>%
  drop_na(Yoi) %>%
  drop_na(Mi) %>%
  group_by(temp, nut, combination) %>%
  mutate(NetEffect = sum(Yoi, na.rm = T)- sum(relAb*Mi))

#plot
ggplot(NetBiodiversity_HI, aes(x = combination, y = NetEffect, color = nut))+
  geom_point()+
  geom_hline(yintercept = 0)+
  facet_grid(~temp, scales = 'free_y')+
  theme_bw()
#----------------------------------------------------------------------------------------------------------------#

## BioVolume ##
BioVolume_exp22 
unique(OD_counts_HI$speciesID)
unique(BioVolume_exp22$speciesID)

#no NAs
expectedStab_HI <- OD_counts_HI%>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = mean.cells * BV)
names(expectedStab_HI)

expectedStab_HI %>%
  filter(sampling == 1) %>%
  filter(species == 'duo') %>%
  ggplot(., aes(x = combination, y = mean.cells, color = speciesID))+
  geom_point()
# relative biomass of species in mxiture at t0
RelBV_t0_HI <- expectedStab_HI %>%
  filter(sampling == 1) %>%
  filter(species == 'duo') %>%
  group_by(nut, combination, temp) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(nut, combination, temp, speciesID, relBV)
names(RelBV_t0_HI)
#observed biomass of species in mixture and monoculture  
counts_duo_tmax_HI <- expectedStab_HI %>%
  filter(sampling == 15) %>%
  filter(species == 'duo') %>%
  rename(Yoi = cellVolume)%>%
  select(nut, combination, temp, speciesID, Yoi)
names(counts_duo_tmax_HI)

counts_mono_tmax_HI <- expectedStab_HI %>%
  filter(sampling == 15) %>%
  filter(species == 'mono') %>%
  rename(Mi = cellVolume) %>%
  select(nut, temp, speciesID, Mi)
names(counts_mono_tmax_HI)
#combine df RelAbund_t0 and OD_counts_tmax
NetBiodiversityBV_HI  <- left_join(RelBV_t0_HI,counts_mono_tmax_HI, by = c('nut',  'temp', 'speciesID')) %>%
  left_join(., counts_duo_tmax_HI, by = c('nut',  'temp', 'speciesID', 'combination')) %>%
  drop_na(Yoi) %>%
  drop_na(Mi) %>%
  group_by(temp, nut, combination) %>%
  mutate(NetEffect = sum(Yoi, na.rm = T)- sum(relBV*Mi))

#plot
ggplot(NetBiodiversityBV_HI, aes(x = combination, y = NetEffect, color = nut))+
  geom_point()+
  geom_hline(yintercept = 0)+
  facet_grid(~temp, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), file = 'NetBiodivEffect_HI.png', width = 10, height = 4)
#----------------------------------------------------------------------------------------------------------------#

ggplot(NetBiodiversityBV_, aes(x = combination, y = compEffect))+
  geom_hline(yintercept = 0, color ='darkgrey')+
  geom_point(size = 2)+
  labs(y = 'Complementarity effect', title = 'Two species in combination')+
  facet_grid(~temp, scales = 'free_y')+
  theme_bw()+  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'bottom')


ggplot(NetBiodiversityBV_, aes(x = combination, y = selecEffect))+
  geom_hline(yintercept = 0, color ='darkgrey')+
  geom_point(size = 2)+
  geom_hline(yintercept = 0)+
  labs(y = 'Selection effect', title = 'Two species in combination')+
  facet_grid(~temp, scales = 'free_y')+
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


#----------------------------------------------------------------------------------------------------------------#


#### Expected Stability HI nut ####

## variables
# M - mono culture of species
mono_t_HI <- expectedStab_HI %>%
  filter(species == 'mono' ) %>%
  filter(temp != 'CS') %>%
  select(-mean.cells,-sd.cells,-se.cells) %>%
  rename(Mono_T = cellVolume)%>%
  select(-species, -Alg, -BV)

mono_control_HI <- expectedStab_HI %>%
  filter(species == 'mono') %>%
  filter(temp == 'CS')%>%
  select(-mean.cells,-sd.cells,-se.cells) %>%
  rename(Mono_C = cellVolume) %>%
  select(-species, -temp, -Alg, -BV)

mono_HI <- mono_t_HI %>% 
  left_join(., mono_control_HI, by = c('speciesID', 'nut', 'sampling', 'combination'))


# Mix 
mix_t_HI <- expectedStab_HI %>%
  filter(species == 'duo') %>%
  filter(temp != 'CS')%>%
  select(-mean.cells,-sd.cells,-se.cells) %>%
  group_by(speciesID,temp, nut, sampling, combination) %>%
  mutate(Mix_T = sum(cellVolume,na.rm = T)) %>%
  ungroup()%>%
  select(-species, -Alg, -BV)

mix_control_HI <- expectedStab_HI %>%
  filter(species == 'duo') %>%
  filter(temp == 'CS')%>%
  select(-mean.cells,-sd.cells,-se.cells) %>%
  group_by( speciesID,nut, sampling, combination) %>%
  mutate(Mix_C = sum(cellVolume,na.rm = T)) %>%
  select(-species, -temp, -Alg, -BV)

mix_HI <- mix_t_HI %>%
  select(-cellVolume) %>%
  left_join(., mix_control_HI, by = c('speciesID', 'nut', 'sampling', 'combination')) %>%
  select(-cellVolume)


### merge all df and calculate response ratios ###
all_HI <- mono_HI %>%
  select(-combination)%>%
  right_join(., mix_HI, by = c('temp','speciesID', 'nut', 'sampling')) %>%
  mutate(ratio = Mono_T/Mono_C) %>%
  group_by(temp, sampling, nut, combination) %>%
  mutate(sumC = sum(Mix_C,na.rm = T),
         sumT = sum(Mix_T,na.rm = T)) %>%
  ungroup() %>%
  left_join(., RelBV_t0_HI, by = c("nut","combination", "temp","speciesID" )) %>%
  mutate(exp = sumC*ratio*relBV)%>%
  mutate(RR_mono_exp = (Mono_T - Mono_C)/(Mono_T+Mono_C),
         RR_obs = (Mix_T - Mix_C)/(Mix_T+Mix_C),
         RR_ges_obs = (sumT-sumC)/(sumT+sumC),
         RR_exp = (exp-relBV*sumC)/(exp+relBV*sumC)) %>%
  group_by(temp, sampling, nut, combination) %>%
  mutate(sumexp = sum(exp, na.rm = T),
         RR_ges_exp = (sumexp-sumC)/(sumexp+sumC))%>%
  ungroup()%>%
  mutate(delta_ges = RR_ges_obs-RR_ges_exp) %>%
  group_by(speciesID, temp, sampling, nut, combination) %>%
  mutate(delta_RR = Mix_T - exp)

all_HI$USI <- paste(all_HI$temp, all_HI$nut,all_HI$combination,all_HI$speciesID, sep = '_')

#### AUC ####
stab.auc <- data.frame()

USI <- unique(all_HI$USI)
for(i in 1:length(USI)){
  temp<-all_HI[all_HI$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.deltaRR<-auc(temp$sampling, temp$delta_RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("spline"),absolutearea = FALSE)
    AUC.ges.RR<-auc(temp$sampling, temp$delta_ges,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("spline"),absolutearea = FALSE)
    stab.auc<-rbind(stab.auc,
                    data.frame(temp[1,c(1:18)],
                               AUC.deltaRR,
                               AUC.ges.RR))
    rm(temp)
  }
}

summary(stab.auc)
stab.auc 
unique(stab.auc$combination)

#### plots####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")

ggplot(subset(stab.auc, nut == 'HI'))+
  geom_col( aes(x = combination, y = AUC.deltaRR, fill = speciesID), alpha=0.8)+
  facet_wrap(~temp)+
  ggtitle('HI nutrients')+
  scale_fill_manual(values=cbbPalette)+
  theme_bw()
ggsave(plot = last_plot(), file = 'delatRR_HINut.png', width = 8,height = 3)

# total effect
ggplot(stab.auc)+
  geom_point(aes(x = combination, y = AUC.ges.RR, color = nut),  alpha=0.8, size = 2)+
  facet_grid(~temp, scales = 'free')+
  geom_hline(yintercept = 0)+
  theme_bw()
ggsave(plot = last_plot(), file = 'TotDeltaRR_HI.png', width = 8, height = 3)

# rel BV at to 
all_HI %>%
  filter(sampling == 1 ) %>%
  distinct(speciesID, combination,temp, relBV) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = relBV, fill = speciesID))+
  geom_hline(yintercept = 0.5, color = 'darkgrey', size = 0.2)+
  scale_fill_manual(values=cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  theme_bw()+
  theme(legend.position = 'bottom')
ggsave(plot = last_plot(), file = 'relBVto_exp22_HINut.png', width = 8, height = 4)


