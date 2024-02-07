#### R script to prepare data for analysis with Dome ####
# 02.10.2022
# by Charlotte Kunze

# packages
library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(RColorBrewer)

setwd("~/Desktop/Exp22/MicrocosmExp22/Data")

# look at the data 
#### import count data ####
expData<- read.csv('Exp22_wrangledData.csv') %>%
  select(-X)
names(expData)
expData[ expData == 'TRUE' ] <- 'T'
str(expData)

#----------------------------------------------------------------------------------------------------------------#

#### Net Biodiversity Effect: BioVolume ####

data<- expData%>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = ''), ifelse(species == 'duo', paste(A1, A2, sep = ''), 
                                                                            ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))))%>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5, -Afield, -At_mm, -plate, -GF, -counts, -magn,-V_ml) %>%
  filter(species %in% c('mono', 'duo')) 

BioVolume_exp22 <- read_excel("BioVolume_exp22.xlsx",sheet = "Sheet1") %>%
  rename(speciesID = species) %>%
  select(-sampling,-comment, -Ansatz)
names(BioVolume_exp22)
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Ditylum']<-'DityCux'
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Thala']<-'ThalaCux'


#no NAs
expectedStab_ <- data%>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = cells_ml * BV) %>%
  select(-cells_ml)

#control data
control_stab <- expectedStab_ %>%
  filter(temp == 'CS') %>%
  group_by(nut, combination, sampling, speciesID ) %>%
  summarise(con.vol = mean(cellVolume, na.rm = T)) 

#treatment data 
expectedStab1 <- expectedStab_ %>%
  filter(temp != 'CS') %>%
  left_join(., control_stab, by = c('nut', 'combination', 'sampling', 'speciesID'))

names(expectedStab1)

expectedStab1 %>%
  filter(sampling == 1) %>%
  filter(species == 'duo') %>%
  ggplot(., aes(x = combination, y = cellVolume, color = speciesID))+
  geom_point()

# relative biomass of species in mxiture at t0
RelBV_t0_ <- expectedStab_ %>%
  filter(sampling == 1) %>%
  group_by(nut, combination, temp, rep) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(nut, combination, temp, speciesID, rep, relBV)
names(RelBV_t0_)

#observed biomass of species in mixture and monoculture  
counts_duo_tmax_ <- expectedStab1 %>%
  filter(sampling == 15) %>%
  filter(species == 'duo') %>%
  rename(Yoi = cellVolume)%>%
  select(nut, combination, temp,rep, speciesID, Yoi)
names(counts_duo_tmax_)

counts_mono_tmax_ <- expectedStab1 %>%
  filter(sampling == 15) %>%
  filter(species == 'mono') %>%
  rename(Mi = cellVolume) %>%
  select(nut, temp, speciesID, rep,Mi)
names(counts_mono_tmax_)

#combine df RelAbund_t0 and OD_counts_tmax
NetBiodiversityBV_  <- left_join(RelBV_t0_,counts_mono_tmax_, by = c('nut',  'temp', 'speciesID', 'rep')) %>%
  left_join(., counts_duo_tmax_, by = c('nut',  'temp', 'speciesID', 'combination', 'rep')) %>%
  filter(nut == 'Low')%>%
  drop_na(Yoi) %>%
  drop_na(Mi) %>%
  group_by(temp, nut, combination, rep) %>%
  mutate(NetEffect = sum(Yoi, na.rm = T)- sum(relBV*Mi)) %>%
  group_by(temp, nut, combination) %>%
  mutate(N = as.numeric(paste(2)),
         compEffect = N*mean(NetEffect)*mean(Mi),
         selecEffect = N*cov(NetEffect, Mi))

# Net effect plot
NetBiodiversityBV_ %>%
  group_by(combination, nut, temp) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n())) %>%
ggplot(., aes(x = combination, y = mean.effect))+
  geom_hline(yintercept = 0, color ='darkgrey')+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect + se.effect), width = .8)+
  geom_point(size = 2)+
  labs(y = 'Mean Net Effect +- SE', title = 'Two species in combination')+
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

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NetBiodivEffect_.png'), width = 10, height = 4)


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

#### Expected Stability Duo ####

## variables
# M - mono culture of species
names(expectedStab1)
mono <- expectedStab1 %>%
  filter(species == 'mono' ) %>%
  filter(temp != 'CS') %>%
  rename(Mono_T = cellVolume)%>%
  rename(Mono_C = con.vol) %>%
  select(-species, -Alg, -BV)

# Mix 
duo <- expectedStab1 %>%
  filter(species == 'duo') %>%
  group_by(speciesID,temp, nut, sampling, combination) %>%
  rename(Mix_T = cellVolume) %>%
  rename(Mix_C = con.vol) %>%
  ungroup()%>%
  select(-species, -Alg, -BV,-no)


### merge all df and calculate response ratios ###
all <- mono %>%
  select(-combination)%>%
  right_join(., duo, by = c('temp','speciesID', 'nut', 'sampling', 'rep')) %>%
  mutate(ratio = Mono_T/Mono_C) %>%
  group_by(temp, sampling, nut, combination,rep) %>%
  mutate(sumC = sum(Mix_C,na.rm = T),
         sumT = sum(Mix_T,na.rm = T)) %>%
  ungroup() %>%
  left_join(., RelBV_t0_, by = c("nut","combination", "temp","speciesID", 'rep' )) %>%
  mutate(exp = sumC*ratio*relBV)%>%
  mutate(RR_mono_exp = (Mono_T - Mono_C)/(Mono_T+Mono_C),
         RR_obs = (Mix_T - Mix_C)/(Mix_T+Mix_C),
         RR_ges_obs = (sumT-sumC)/(sumT+sumC),
         RR_exp = (exp-relBV*sumC)/(exp+relBV*sumC)) %>%
  group_by(temp, sampling, nut, combination,rep) %>%
  mutate(sumexp = sum(exp, na.rm = T),
         RR_ges_exp = (sumexp-sumC)/(sumexp+sumC))%>%
  ungroup()%>%
  mutate(delta = RR_obs - RR_exp) %>%
  mutate(delta_ges = RR_ges_obs-RR_ges_exp) %>%
  group_by(speciesID, temp, sampling, nut, combination,rep) %>%
  mutate(delta_RR = Mix_T - exp) 

all$RR_ges_exp[is.na(all$RR_ges_exp)]<-0
all$RR_ges_obs[is.na(all$RR_ges_obs)]<-0

all$delta[is.na(all$delta)]<-0
all$USI <- paste(all$temp, all$nut,all$combination,all$speciesID, all$rep, sep = '_')

str(all)
#### AUC ####
stab.auc <- data.frame()

USI <- unique(all$USI)
for(i in 1:length(USI)){
  temp<-all[all$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.deltaRR<-auc(temp$sampling, temp$delta_RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("spline"),absolutearea = FALSE)
    AUC.RR_obs<-auc(temp$sampling, temp$RR_ges_obs,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("spline"),absolutearea = FALSE)
    AUC.RR_exp<-auc(temp$sampling, temp$RR_ges_exp,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("spline"),absolutearea = FALSE)
    AUC.delta<-auc(temp$sampling, temp$delta,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("spline"),absolutearea = FALSE)
    AUC.ges.RR<-auc(temp$sampling, temp$delta_ges,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("spline"),absolutearea = FALSE)
    stab.auc<-rbind(stab.auc,
                    data.frame(temp[1,c(1:19)],
                               AUC.deltaRR,                               
                               AUC.delta,
                               AUC.ges.RR,
                               AUC.RR_exp,AUC.RR_obs))
    rm(temp)
  }
}

summary(stab.auc)

stab.auc 
unique(stab.auc$combination)

#### plots####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")
stab.auc$temp[stab.auc$temp=='fluct'] <- 'Fluctuation'
stab.auc$temp[stab.auc$temp=='inc'] <- 'Increase'
stab.auc$temp[stab.auc$temp=='inc+fluc'] <- 'IncreaseFluctuation'

stab.auc %>%
  filter(nut == 'Low') %>%
  group_by(combination, speciesID, temp, nut)%>%
  summarise(mean = mean(AUC.deltaRR,na.rm = T),
            sd = sd(AUC.deltaRR,na.rm = T),
            se = sd/sqrt(n())) %>%
ggplot(., aes(x = combination, y = mean, fill = speciesID, color = speciesID))+
  geom_hline(yintercept = 0, color = 'darkgrey') +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .8)+
  geom_col(alpha=0.8)+
  facet_grid(~temp~speciesID, scales = "free")+
  labs(title = '', y = 'AUC(MixT - exp)')+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'bottom')

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/deltaRR_grid.png'), width = 10,height = 7)

# total effect
stab.auc %>%
  filter(nut == 'Low')%>%
group_by(combination, temp, nut)%>%
  summarise(mean = mean(AUC.delta,na.rm = T),
            sd = sd(AUC.delta,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = combination, y =mean))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .8)+
  geom_point(alpha=0.8, size = 2, color = 'black')+
  facet_grid(~temp, scales = 'free')+
  labs(y = 'Net Biodiversity effect on Stability')+
  scale_colour_brewer(palette = "Set1")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/TotDeltaRR_low.png'), width = 8, height = 3)


### complementarity & selection ###
respVar <- stab.auc%>%
  filter(nut == 'Low')%>%
  drop_na(AUC.deltaRR,RR_mono_exp)%>%
  group_by(temp, combination,nut) %>%
  mutate(N = as.numeric(paste(2)),
         compEffect = N*mean(AUC.deltaRR)*mean(RR_mono_exp))%>%
  group_by(temp, combination,nut) %>%
  mutate(selecEffect = N*cov(AUC.deltaRR,RR_mono_exp))

ggplot(respVar, aes(x = combination, y = compEffect))+
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


ggplot(respVar, aes(x = combination, y = selecEffect))+
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




# rel BV at to 
all %>%
  filter(sampling == 1 & nut == 'Low') %>%
  group_by(speciesID, combination,temp) %>%
  summarise(mean.BV = mean(relBV, na.rm = T)) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = mean.BV, fill = speciesID))+
  geom_hline(yintercept = 0.5, color = 'darkgrey', size = 0.2)+
  scale_fill_manual(values=cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  labs( y= 'Mean relative BV')+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/relBVto_exp22_LowNut.png'), width = 8, height = 4)


#### Expected Stability: Quattro and Mix ####

# M - mono culture of species
names(expectedStab1)
mono <- expectedStab1 %>%
  filter(species == 'mono' ) %>%
  filter(temp != 'CS') %>%
  rename(Mono_T = cellVolume)%>%
  rename(Mono_C = con.vol) %>%
  select(-species, -Alg, -BV)

# Mix 
names(expData)
mix <- expData %>%
  filter(species %in% c('MIX', 'quattro')) %>%
  mutate(combination = ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))%>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5, -Afield, -At_mm, -plate, -GF, -counts, -magn,-V_ml) %>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = cells_ml * BV) %>%
  select(-cells_ml) 

mix_con <- mix %>%
  filter(temp == 'CS') %>%
  rename(Mix_C = cellVolume) %>%
  ungroup()%>%
  select(-species, -Alg, -BV,-no,-temp)
str(mix_con)

mix_all <- mix %>%
  filter(temp != 'CS') %>%
  rename(Mix_T = cellVolume) %>%
  left_join(., mix_con, by = c("speciesID", "nut","rep","sampling","combination"))
  

### merge all df and calculate response ratios ###

relBVt0 <- expData %>%
  filter(species %in% c('MIX', 'quattro')) %>%
  mutate(combination = ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))%>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5, -Afield, -At_mm, -plate, -GF, -counts, -magn,-V_ml) %>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = cells_ml * BV) %>%
  select(-cells_ml) %>%
  filter(sampling == 1) %>%
  group_by(nut, combination, temp, rep) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(nut, combination, temp, speciesID, rep, relBV)

allMix <- mono %>%
  select(-combination,-no)%>%
 # filter(nut == 'Low') %>%
  right_join(., mix_all, by = c('temp','speciesID', 'nut', 'sampling', 'rep')) %>%
  mutate(ratio = Mono_T/Mono_C) %>%
  group_by(temp, sampling, nut, combination,rep) %>%
  mutate(sumC = sum(Mix_C,na.rm = T),
         sumT = sum(Mix_T,na.rm = T)) %>%
  ungroup() %>%
  left_join(., relBVt0, by = c("nut","combination", "temp","speciesID", 'rep')) %>%
  mutate(exp = sumC*ratio*relBV)%>%
  mutate(RR_mono_exp = (Mono_T - Mono_C)/(Mono_T+Mono_C),
         RR_obs = (Mix_T - Mix_C)/(Mix_T+Mix_C),
         RR_ges_obs = (sumT-sumC)/(sumT+sumC),
         RR_exp = (exp-relBV*sumC)/(exp+relBV*sumC)) %>%
  group_by(temp, sampling, nut, combination,rep) %>%
  mutate(sumexp = sum(exp, na.rm = T),
         RR_ges_exp = (sumexp-sumC)/(sumexp+sumC))%>%
  ungroup()%>%
  mutate(delta_ges = RR_ges_obs-RR_ges_exp) %>%
  group_by(speciesID, temp, sampling, nut, combination,rep) %>%
  mutate(delta_RR = Mix_T - exp)

allMix$USI <- paste(allMix$temp, allMix$nut,allMix$combination,allMix$speciesID, allMix$rep, sep = '_')

#### AUC ####
stab.auc.mix <- data.frame()

USI <- unique(allMix$USI)
for(i in 1:length(USI)){
  temp<-allMix[allMix$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.deltaRR<-auc(temp$sampling, temp$delta_RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("spline"),absolutearea = FALSE)
    AUC.ges.RR<-auc(temp$sampling, temp$delta_ges,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("spline"),absolutearea = FALSE)
    AUC.RR_exp<-auc(temp$sampling, temp$RR_ges_exp,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("spline"),absolutearea = FALSE)
    AUC.RR_obs<-auc(temp$sampling, temp$RR_ges_obs,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("spline"),absolutearea = FALSE)
    stab.auc.mix<-rbind(stab.auc.mix,
                    data.frame(temp[1,c(1:18,22)],
                               AUC.deltaRR,
                               AUC.ges.RR,AUC.RR_obs,AUC.RR_exp))
    rm(temp)
  }
}

summary(stab.auc.mix)

stab.auc.mix 

unique(stab.auc.mix$combination)

#### plots####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")
stab.auc.mix$temp[stab.auc.mix$temp=='fluct'] <- 'Fluctuation'
stab.auc.mix$temp[stab.auc.mix$temp=='inc'] <- 'Increase'
stab.auc.mix$temp[stab.auc.mix$temp=='inc+fluc'] <- 'IncreaseFluctuation'

stab.auc.mix$combination[stab.auc.mix$combination == 'ADGRT'] <- 'Mix'
stab.auc.mix %>%
  filter(nut == 'Low') %>%
  group_by(combination, speciesID, temp, nut)%>%
  summarise(mean = mean(AUC.deltaRR,na.rm = T),
            sd = sd(AUC.deltaRR,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = combination, y = mean, fill = speciesID, color = speciesID))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .8)+
  geom_col(alpha=0.8)+
  facet_grid(~temp~speciesID, scales = 'free')+
  labs(title = '', y = 'Selection effect')+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  theme(legend.position = 'bottom')+
  guides(color = guide_legend(override.aes = list(size = 3.5)))

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/deltaRR_quattroSpec.png'), width = 11,height = 7)

#pannels for temp only
stab.auc.mix %>%
  filter(nut == 'Low') %>%
  group_by(combination, speciesID, temp, nut)%>%
  summarise(mean = mean(AUC.deltaRR,na.rm = T),
            sd = sd(AUC.deltaRR,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = combination, y = mean, fill = speciesID, color = speciesID))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .8)+
  geom_col(alpha=0.8)+
  facet_grid(~temp~speciesID, scales = 'free')+
  labs( y = 'AUC(MixT - exp)')+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  theme(legend.position = 'bottom')+
  guides(color = guide_legend(override.aes = list(size = 3.5)))

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/deltaRR_quattroSpec_all.png'), width = 11,height = 8)


# total effect
stab.auc.mix %>%
  filter(nut == 'Low')%>%
  group_by(combination, speciesID, temp, nut)%>%
  summarise(mean = mean(AUC.ges.RR,na.rm = T),
            sd = sd(AUC.ges.RR,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = combination, y =mean, col = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .2, alpha = 0.6)+
  geom_point(alpha=0.6, size = 2.5)+
  facet_grid(~temp, scales = 'free')+
  labs(y = 'Net Interaction effect')+
  theme_bw()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/TotDeltaRR_quattro.png'), width = 8, height = 3)

# rel BV at to 

allMix$combination[allMix$combination == 'ADGRT'] <- 'Mix'
allMix %>%
  filter(sampling == 1 & nut == 'Low') %>%
  group_by(speciesID, combination,temp) %>%
  summarise(mean.BV = mean(relBV, na.rm = T)) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = mean.BV, fill = speciesID))+
  scale_fill_manual(values=cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  labs( y= 'Mean relative BV')+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/QuattroMixrelBVto_exp22_LowNut.png'), width = 9, height = 4)


#### Netto Biodiversity Effect Mix ####
data2<- expData%>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = ''), ifelse(species == 'duo', paste(A1, A2, sep = ''), 
                                                                            ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste('Mix')))))%>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5, -Afield, -At_mm, -plate, -GF, -counts, -magn,-V_ml) %>%
  filter(species %in% c('quattro', 'MIX','mono')) 

#no NAs
expectedStab2 <- data2%>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = cells_ml * BV) %>%
  select(-cells_ml)

#control data
control_stab2 <- expectedStab2 %>%
  filter(temp == 'CS') %>%
  group_by(nut, combination, sampling, speciesID ) %>%
  summarise(con.vol = mean(cellVolume, na.rm = T)) 

#treatment data 
expectedStab2 <- expectedStab2 %>%
  filter(temp != 'CS') %>%
  left_join(., control_stab, by = c('nut', 'combination', 'sampling', 'speciesID'))

names(expectedStab2)

expectedStab2 %>%
  filter(sampling == 1) %>%
  ggplot(., aes(x = combination, y = cellVolume, color = speciesID))+
  geom_point()

# relative biomass of species in mxiture at t0
RelBV_t0_2 <- expectedStab2 %>%
  filter(sampling == 1) %>%
  group_by(nut, combination, temp, rep) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(nut, combination, temp, speciesID, rep, relBV)
names(RelBV_t0_2)
#observed biomass of species in mixture and monoculture  
counts_duo_tmax_2 <- expectedStab2 %>%
  filter(sampling == 15) %>%
  filter(!combination %in% c('A', 'D', 'G','R','T')) %>%
  rename(Yoi = cellVolume)%>%
  select(nut, combination, temp,rep, speciesID, Yoi)
names(counts_duo_tmax_2)

counts_mono_tmax_2 <- expectedStab2 %>%
  filter(sampling == 15) %>%
  filter(species == 'mono') %>%
  rename(Mi = cellVolume) %>%
  select(nut, temp, speciesID, rep,Mi)
str(counts_mono_tmax_2)
str(RelBV_t0_2)
str(counts_duo_tmax_)
#combine df RelAbund_t0 and OD_counts_tmax
NetBiodiversityBV2  <- merge(RelBV_t0_2,counts_mono_tmax_2, by = c('nut',  'temp', 'speciesID', 'rep')) %>%
  left_join(., counts_duo_tmax_2, by = c('nut',  'temp', 'speciesID', 'combination', 'rep')) %>%
  filter(nut == 'Low')%>%
  filter(!combination %in% c('A', 'D', 'G','R','T')) %>%
  drop_na(Yoi) %>%
  drop_na(Mi) %>%
  group_by(temp, nut, combination, rep) %>%
  mutate(NetEffect = sum(Yoi, na.rm = T)- sum(relBV*Mi)) %>%
  group_by(temp, nut, combination) %>%
  mutate(N = as.numeric(paste(ifelse(combination == 'Mix', 5, 4))),
         compEffect = N*mean(NetEffect)*mean(Mi),
         selecEffect = N*cov(NetEffect, Mi))
names(NetBiodiversityBV2)
# Net effect plot
NetBiodiversityBV2 %>%
  group_by(combination, nut, temp) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n())) %>%
  ggplot(., aes(x = combination, y = mean.effect))+
  geom_hline(yintercept = 0, color ='darkgrey')+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect + se.effect), width = .8)+
  geom_point(size = 2)+
  labs(y = 'Mean Net Effect +- SE')+
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

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NetBiodivEffect_MixQuattro.png'), width = 10, height = 4)

names(NetBiodiversityBV_)
allNetBiodiv <- bind_rows(NetBiodiversityBV2, NetBiodiversityBV_)%>%
  group_by( nut, temp, N) %>%
  summarise(mean.effect = mean(NetEffect, na.rm = T),
            sd.effect = sd(NetEffect,na.rm = T),
            se.effect = sd.effect/sqrt(n()))
ggplot(allNetBiodiv, aes(x = N, y = mean.effect,color = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 0.7)+
  geom_errorbar(aes(ymin = mean.effect - se.effect, ymax = mean.effect +se.effect), width = .3)+
  labs(x = 'combinations', y = 'Net Biodiversity effect', color = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
 # facet_grid(~N, scales = 'free')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NetBiodiv_all.png'), width = 5, height = 4)

#### total effect plot ####


d1 <- stab.auc %>%
  filter(nut == 'Low')%>%
  select(combination, speciesID, temp, AUC.ges.RR, AUC.RR_exp,AUC.RR_obs)
d2 <- stab.auc.mix %>%
  filter(nut == 'Low')%>%
  select(combination, speciesID, temp, nut,AUC.ges.RR,AUC.RR_exp,AUC.RR_obs)

d2$combination[d2$combination=='Mix']<-'ADGRT'
d3 <- bind_rows(d1, d2) %>%
  mutate(N = str_length(combination)  ) 

OE <- d3%>%
  distinct(combination, temp,AUC.RR_exp,AUC.RR_obs) %>%
  gather(key = 'AUC', value = 'AUC.value',-combination, -temp) %>%
  mutate(N = str_length(combination)  ) %>%
  group_by(temp, N, AUC,combination) %>%
  mutate(mean = mean(AUC.value, na.rm = T),
         sd = sd(AUC.value, na.rm = T),
         se = sd/sqrt(n()))
  ggplot(OE)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(x = combination, y = AUC.value,color = AUC), size = 1, alpha = 0.4)+
  geom_point(aes(x = combination, y = mean,color = AUC),size = 2.5)+
   geom_errorbar(aes(x = combination, y = mean,color = AUC,ymin = mean - se, ymax = mean +se), width = .3)+
  labs(x = 'combination', y = 'Net Interaction effect', color = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~temp~N, scales = 'free_x')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
ggsave(plot = last_plot(), file = here('~/Desktop/Exp22/MicrocosmExp22/output/expectedObserved.png'), width = 10, height = 8)




p1<-d3%>%
  group_by(N, temp,combination) %>%
  summarise(mean.total.DRR = mean(AUC.ges.RR, na.rm = T),
             sd = sd(AUC.ges.RR,na.rm = T),
            se = sd/sqrt(n()))%>%
ggplot(., aes(x = combination, y = mean.total.DRR,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 0.8)+
  geom_errorbar(aes(ymin = mean.total.DRR - se, ymax = mean.total.DRR +se), width = .3)+
  labs(x = 'Combination', y = 'Net Biodiversity Effect on Stability', color = 'Treatment')+
 # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~N, scales = 'free')+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
p1
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/FacetNumberofConn.png'), width = 10, height = 3)

p2 <- d3%>%
  group_by(N, temp) %>%
  summarise(mean.total.DRR = mean(AUC.ges.RR, na.rm = T),
            sd = sd(AUC.ges.RR,na.rm = T),
            se = sd/sqrt(n()))%>%
  ggplot(., aes(x = as.factor(N), y = mean.total.DRR,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 0.8)+
  geom_errorbar(aes(ymin = mean.total.DRR - se, ymax = mean.total.DRR +se), width = .2)+
  labs(x = 'Species Richness', y = 'Net Biodiversity Effect on Stability', color = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  #facet_grid(~N, scales = 'free')+
  theme_bw()+
  theme(legend.position = 'none',
          panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
p2
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NumberofConn.png'), width = 3, height = 3)

library(cowplot)
plot_grid( p2,p1, labels = c('a)', 'b)'), rel_widths = c( 1/3,2/3))
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/plot-grid_fig2.png'),width = 12, height = 4)

#### Figure RR ####

data <- select(stab.auc,combination, speciesID, temp, nut,rep,AUC.RR_obs) %>%
  bind_rows(., select(stab.auc.mix,combination, speciesID, rep,temp,nut, AUC.RR_obs)) %>%
  filter(nut == 'Low') %>%
  mutate(N = str_length(combination),
         N = paste(ifelse(N == 3, 5, N)) )%>%
  group_by(temp, combination, N)%>%
  mutate(mean = mean(AUC.RR_obs, na.rm = T),
         sd = sd(AUC.RR_obs, na.rm = T),
         se = sd/sqrt(n()))
ggplot(data)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
 # geom_point(aes(x = combination, y = AUC.RR_obs), size = 1, alpha = 0.4)+
  geom_point(aes(x = combination, y = mean, color = temp),size = 2.5)+
  geom_errorbar(aes(x = combination, y = mean,ymin = mean - se, ymax = mean +se, color = temp), width = .3)+
  labs(x = 'combination', y = 'AUC RR observed', color = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~N, scales = 'free_x')+
  theme_bw()+
  theme(legend.position = 'bottom',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
ggsave(plot = last_plot(), file = here('~/Desktop/Exp22/MicrocosmExp22/output/AUCRRObserved.png'), width = 10, height = 8)






#### TestData: 1 combination and nut#### 

test <- all %>%
  filter(combination == 'RT'&temp == 'inc'&nut == 'Low' & rep == 2) 
names(test)
test %>%
  filter(sampling == 1) %>%
  distinct(speciesID, combination, relBV) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = relBV, fill = speciesID))+
  scale_fill_manual(values=c("#009E73","#CC79A7"))+
  theme_bw()
#ggsave(plot = last_plot(), file = 'RT_test.png', width = 3, height = 3)
ggplot(test )+
  geom_point(aes(x = sampling, y = Mix_T, shape = speciesID), position = position_dodge(width = 0.5),color = '#FFDB58', alpha=0.8, size = 2)+
  geom_line(aes(x = sampling, y = Mix_T,  group = speciesID), color = '#FFDB58',alpha=0.8, size = 0.5)+
 # geom_point(aes(x = sampling, y = Mono_T,  shape = speciesID), color = 'black',alpha=0.8, size = 2)+
 # geom_line(aes(x = sampling, y = Mono_T,  group = speciesID), color = 'black',alpha=0.8, size = 0.5)+
  geom_point(aes(x = sampling, y = exp, shape = speciesID), color = 'red', alpha=0.5, size = 2)+
  geom_line(aes(x = sampling, y = exp, group = speciesID), color = 'red', alpha=0.5, size = 0.5)+
  labs(title = 'Low nut, rep 3, inc temp', y = 'BioVolume')+
  theme_bw()
ggsave(plot = last_plot(), file = here('output/RT_2.png'), width = 4, height = 2)


#total RR

ggplot(test )+
   geom_point(aes(x = sampling, y = RR_ges_obs), color = 'black',alpha=0.8, size = 2)+
  geom_line(aes(x = sampling, y = RR_ges_obs), size = 0.5, linetype = 'dashed')+
  geom_point(aes(x = sampling, y = RR_ges_exp), color = 'red', alpha=0.5, size = 2)+
  geom_line(aes(x = sampling, y = RR_ges_exp), color = 'red',size = 0.5, linetype = 'longdash')+
  labs(title = 'Low nut, rep 3, inc temp', y = 'Total observed and expected stability')+
  theme_bw()
ggsave(plot = last_plot(), file = here('output/RRgesExpObs.png'), width = 4, height = 2)

### AUC  ###
stab.auc <- data.frame()

USI <- unique(test$speciesID)
for(i in 1:length(USI)){
  temp<-test[test$speciesID==USI[i], ]#creates a temporary data frame for each case
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

ggplot(stab.auc)+
  geom_col( aes(x = speciesID, y = AUC.deltaRR, fill = speciesID))+
  scale_fill_manual(values=c("#009E73","#CC79A7"))+
  geom_point(aes(x = speciesID, y = AUC.ges.RR), color = 'black')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'AUCs_RT1.png', width = 4,height = 2)


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

#### Expected Stability ####

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




