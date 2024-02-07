### R script to prepare data for analysis with Dome ###
# 02.10.2022
# by Charlotte Kunze

#### packages ####
library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(RColorBrewer)
library(cowplot)

setwd("~/Desktop/Exp22/MicrocosmExp22/Data")

# look at the data 

#### import count data ####
expData<- read.csv('Exp22_wrangledData.csv') %>%
  select(-X)
names(expData)
expData[ expData == 'TRUE' ] <- 'T'
str(expData)


#### BioVolume ####

BioVolume_exp22 <- read_excel("BioVolume_exp22.xlsx",sheet = "Sheet1") %>%
  rename(speciesID = species) %>%
  select(-sampling,-comment, -Ansatz)
names(BioVolume_exp22)
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Ditylum']<-'DityCux'
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Thala']<-'ThalaCux'
#

#----------------------------------------------------------------------------------------------------------------#

data<- expData%>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = ''), ifelse(species == 'duo', paste(A1, A2, sep = ''), 
                                                                            ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))))%>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5, -Afield, -At_mm, -plate, -GF, -counts, -magn,-V_ml) %>%
  filter(species %in% c('mono', 'duo')) 

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
    AUC.RR_mono_exp <- auc(temp$sampling, temp$RR_mono_exp, from = min(temp$sampling, na.rm = T), to = max(temp$sampling, na.rm = T),
                           type = c('spline'), absolutearea = F)
    NBE <- AUC.RR_obs - AUC.RR_exp
    stab.auc<-rbind(stab.auc,
                    data.frame(temp[1,c(1:19)],
                               NBE,
                               AUC.deltaRR,                               
                               AUC.delta,
                               AUC.ges.RR,
                               AUC.RR_exp,AUC.RR_obs, AUC.RR_mono_exp))
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

#### SGC data ####
str(stab.auc)
species12<- stab.auc %>%
  filter(nut == 'Low') %>%
  
  group_by(speciesID,temp) %>%
  summarise(mean.NBE = mean(AUC.RR_obs),
            sd = sd(NBE),
            se.NBE = sd/sqrt(n()))
write.csv2(species12, file = 'species12_SGC.csv')

stab.auc %>%
  filter(nut == 'Low') %>%
  group_by(combination, speciesID, temp, nut)%>%
  summarise(mean = mean(NBE,na.rm = T),
            sd = sd(NBE,na.rm = T),
            se = sd/sqrt(n())) %>%
ggplot(., aes(x = combination, y = mean, fill = speciesID, color = speciesID))+
  geom_hline(yintercept = 0, color = 'darkgrey') +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .5, alpha = .7)+
  geom_point(alpha=0.6, size = 2)+
  facet_grid(~temp, scales = "free")+
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

ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/deltaRR_grid.png'), width = 8,height = 5)

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



# rel BV at to 
BV_t0_duo <- all %>%
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
  theme(legend.position = 'none')
BV_t0_duo
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
    NBE <- AUC.RR_obs-AUC.RR_exp
    stab.auc.mix<-rbind(stab.auc.mix,
                    data.frame(temp[1,c(1:18,22)],
                               NBE,AUC.deltaRR,
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

#pannels for temp only
stab.auc.mix %>%
  filter(nut == 'Low') %>%
  group_by(combination, speciesID, temp, nut)%>%
  summarise(mean = mean(AUC.deltaRR,na.rm = T),
            sd = sd(AUC.deltaRR,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = combination, y = mean, fill = speciesID, color = speciesID))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .3)+
  geom_point(size = 2, alpha=0.8)+
  facet_grid(~temp, scales = 'free')+
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

ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/deltaRR_quattroSpec_chaos.png'), width = 8,height = 5)


#### 2- quattro total effect ####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")
stab.auc.mix$Missing <- stab.auc.mix$combination
stab.auc.mix$Missing[stab.auc.mix$Missing =="AGRT"] <- 'D'
stab.auc.mix$Missing[stab.auc.mix$Missing =="ADGT"] <- 'R'
stab.auc.mix$Missing[stab.auc.mix$Missing =="ADRT"] <- 'G'
stab.auc.mix$Missing[stab.auc.mix$Missing =="ADGR"] <- 'T'
stab.auc.mix$Missing[stab.auc.mix$Missing =="DGRT"] <- 'A'
stab.auc.mix$Missing[stab.auc.mix$Missing =="Mix"] <- 'Xnone'

unique(stab.auc.mix$Missing)

stab.auc.mix %>%
  filter(nut == 'Low')%>%
  filter(Missing != 'Xnone')%>%
  group_by(Missing, speciesID, temp, nut)%>%
  summarise(mean = mean(AUC.ges.RR,na.rm = T),
            sd = sd(AUC.ges.RR,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = Missing, y =mean, col = Missing))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .2, alpha = 0.6)+
  geom_point(alpha=0.6, size = 2.5)+
  scale_colour_manual(values = cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  labs(y = 'Net Biodiv effect')+
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

ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NBE_quattro.png'), width = 8, height = 4)

# rel BV at to 
allMix$combination[allMix$combination == 'ADGRT'] <- 'Mix'
BV_t0_mix <- allMix %>%
  filter(sampling == 1 & nut == 'Low') %>%
  group_by(speciesID, combination,temp) %>%
  summarise(mean.BV = mean(relBV, na.rm = T)) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = mean.BV, fill = speciesID))+
  scale_fill_manual(values=cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  labs( y= 'Mean relative BV')+
  theme_bw()+
  theme(legend.position = 'none')
BV_t0_mix
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/QuattroMixrelBVto_exp22_LowNut.png'), width = 9, height = 4)

plot_grid(BV_t0_duo, BV_t0_mix, labels = c('a)', 'b)'), ncol = 1)
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/RelativeBVatT0.png'), width = 8, height = 8)



#### total effect plot ####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )

d1 <- stab.auc %>%
  filter(nut == 'Low')%>%
  select(combination, speciesID, rep,temp, nut,NBE,AUC.ges.RR, AUC.RR_exp,AUC.RR_obs, AUC.RR_mono_exp)
d2 <- stab.auc.mix %>%
  filter(nut == 'Low')%>%
  select(combination, speciesID, rep,temp, nut,NBE,AUC.ges.RR,AUC.RR_exp,AUC.RR_obs)

d2$combination[d2$combination=='Mix']<-'ADGRT'
d3 <- d1 %>%
  select(-AUC.RR_mono_exp)%>%
  bind_rows(., d2) %>%
  mutate(N = str_length(combination)  )  

unique(d3$nut)

data <- d3 %>%
  distinct(combination, rep, temp, nut,N, NBE)
write.csv(data, file = here('MicrocosmExp22/output/netEffect.csv'))

#### START ####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )

d3 <- read.csv2('~/Exp22/MicrocosmExp22/output/netEffect.csv')




OE <- d3%>%
  distinct(combination, temp,AUC.RR_exp,AUC.RR_obs) %>%
  gather(key = 'AUC', value = 'AUC.value',-combination, -temp) %>%
  mutate(N = str_length(combination)  ) %>%
  group_by(temp, N, AUC,combination) %>%
  mutate(mean = mean(AUC.value, na.rm = T),
         sd = sd(AUC.value, na.rm = T),
         se = sd/sqrt(n())) #%>%
  #filter(AUC == 'AUC.RR_obs')
OE$N[OE$N == 2] <- 'Two species'
OE$N[OE$N == 4] <- 'Four species'
OE$N[OE$N == 5] <- 'Five species'

ggplot(OE)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(x = combination, y = AUC.value,color = AUC),alpha = 0.4)+
  geom_errorbar(aes(x = combination, y = mean,color = AUC,ymin = mean - se, ymax = mean +se), alpha = 0.8,width = .3)+
  geom_point(aes(x = combination, y = mean,color = AUC),alpha = 0.8,size = 2.5)+
  labs(x = 'combination', y = 'AUC.RR', color = 'Treatment', shape='Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_manual(values = tempPalette)+
  facet_grid(~temp~N, scales = 'free_x')+
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
#ggsave(plot = last_plot(), file = here('~/Desktop/Exp22/MicrocosmExp22/output/expectedObserved.png'), width = 8, height = 8)




p1<-data%>%
  group_by(N, temp,combination) %>%
  summarise(mean.total.DRR = mean(NBE, na.rm = T),
             sd = sd(NBE,na.rm = T),
            se = sd/sqrt(n()))%>%
  mutate(label = paste(ifelse( N == 2, '2 species', ifelse(N == 4, '4 species', '5 species'))))%>%
ggplot(., aes(x = combination, y = mean.total.DRR,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 0.8)+
  geom_errorbar(aes(ymin = mean.total.DRR - se, ymax = mean.total.DRR +se), width = .3)+
  labs(x = 'Combination', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
 # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~label, scales = 'free')+
  theme_bw()+
  theme(legend.position = 'bottom',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
p1
legendb<-get_legend(p1)
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/FacetNumberofConn.png'), width = 10, height = 3)

p2 <- data%>%
  group_by(N, temp) %>%
  mutate(mean.total.DRR = mean(NBE, na.rm = T),
            sd = sd(NBE,na.rm = T),
            se = sd/sqrt(n()))%>%
  ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(x = as.factor(N), y = AUC.ges.RR,color = temp, shape = temp),size = 0.4, alpha = 0.3)+
 # geom_boxplot(aes(x = as.factor(N), y = AUC.ges.RR, col = temp),size = 1, alpha = 0.8)+
    geom_point(aes(x = as.factor(N), y = mean.total.DRR, color=temp, shape=temp),size = 2.7, alpha = 0.5)+
  geom_errorbar(aes(x=as.factor(N),ymin = mean.total.DRR - se, ymax = mean.total.DRR +se, col=temp), alpha = 0.5, width = .1)+
  labs(x = 'Species Richness', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
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
blank <- ggplot()+
  geom_blank()+
  theme_void() 
  
plot_grid( p2,p1+theme(legend.position = 'none'),blank,blank, labels = c('a)', 'b)'), ncol = 2,legendb,rel_widths = c( 1/3,2/3), rel_heights = c(1,0.3))
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/plot-grid_fig2.jpeg'),width = 12, height = 9)



#### Figure RR ####
names(d1)
MonoData <- d1%>%
  distinct( speciesID, temp,AUC.RR_mono_exp) %>%
  rename(AUC.RR_obs= AUC.RR_mono_exp)%>%
  mutate(combination = paste(ifelse(speciesID == 'Asterio', 'A', ifelse(speciesID == 'DityCux', 'D', ifelse(speciesID == 'Guido', 'G', ifelse(speciesID == 'Rhizo', 'R', 'T')))))) %>%
  bind_rows(., select(d1,combination, speciesID,temp, AUC.RR_obs)) %>%
  bind_rows(., select(d2,combination, speciesID,temp, AUC.RR_obs)) %>%
  mutate(N = str_length(combination),
         N = paste(ifelse(N == 3, 5, N)),
         label = paste(ifelse(N == 1, 'Monoculture', ifelse( N == 2, '2 species', ifelse(N == 4, '4 species', '5 species')))))%>%
  group_by(temp, combination, N)%>%
  mutate(mean = mean(AUC.RR_obs, na.rm = T),
         sd = sd(AUC.RR_obs, na.rm = T),
         se = sd/sqrt(n()))

#change label to factor to adjust coloring
MonoData$label <- factor(as.factor(MonoData$label),levels=c("Monoculture","2 species",   "4 species" ,  "5 species"))

RR1 <- MonoData %>%
  filter(label == 'Monoculture')%>%
    ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(x = combination, y = mean,ymin = mean - se, ymax = mean +se, color = temp), width = .3, alpha = .7)+
  geom_point(aes(x = combination, y = mean, color = temp, shape = temp),size = 2.5, alpha = .6)+
  labs(x = 'combination', y = 'observed (In-)stability (AUC.RR)', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  facet_wrap(~label,  ncol =2)+
  scale_y_continuous(limits = c(-7,7))+
  theme_bw()+
  theme(legend.position = 'right',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
RR1

RR2 <- MonoData %>%
  filter(label == '2 species')%>%
  ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(x = combination, y = mean,ymin = mean - se, ymax = mean +se, color = temp), width = .3, alpha = .7)+
  geom_point(aes(x = combination, y = mean, color = temp, shape = temp),size = 2.5, alpha = .6)+
  labs(x = 'combination', y = '', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  facet_wrap(~label,  ncol =2)+
  scale_y_continuous(limits = c(-7,7))+
  theme_bw()+
  theme(legend.position = 'right',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
RR2

RR4 <- MonoData %>%
  filter(label == '4 species')%>%
  ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(x = combination, y = mean,ymin = mean - se, ymax = mean +se, color = temp), width = .3, alpha = .7)+
  geom_point(aes(x = combination, y = mean, color = temp, shape = temp),size = 2.5, alpha = .6)+
  labs(x = 'combination', y = '', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  facet_wrap(~label,  ncol =2)+
  scale_y_continuous(limits = c(-7,7))+
  theme_bw()+
  theme(legend.position = 'right',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
RR4

RR5 <- MonoData %>%
  filter(label == '5 species')%>%
  ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(x = combination, y = mean,ymin = mean - se, ymax = mean +se, color = temp), width = .3, alpha = .7)+
  geom_point(aes(x = combination, y = mean, color = temp, shape = temp),size = 2.5, alpha = .6)+
  labs(x = 'combination', y = 'observed (In-)stability (AUC.RR)', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  facet_wrap(~label,  ncol =2)+
  scale_y_continuous(limits = c(-7,7))+
  theme_bw()+
  theme(legend.position = 'right',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
RR5

RR1+RR2+RR5+RR4+
  plot_layout(widths = c( 2,3), nrow = 2)+ plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A") 
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/plot-grid_fig2_patchwork.png'),width = 10, height = 8)
######

d3$combination[d3$combination == 'ADGR'] <- 'T'
d3$combination[d3$combination == 'ADGT'] <- 'R'
d3$combination[d3$combination == 'DGRT'] <- 'A'
d3$combination[d3$combination == 'AGRT'] <- 'D'
d3$combination[d3$combination == 'ADRT'] <- 'G'
d3$combination[d3$combination == 'none'] <- 'Xnone'

d3%>%
  group_by(N, temp,combination) %>%
  summarise(mean.total.DRR = mean(AUC.ges.RR, na.rm = T),
            sd = sd(AUC.ges.RR,na.rm = T),
            se = sd/sqrt(n()))%>%
  mutate(label = paste(ifelse( N == 2, '2 species', ifelse(N == 4, '4 species', '5 species'))))%>%
  filter(N != 2) %>%
  ggplot(., aes(x = combination, y = mean.total.DRR,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 0.8)+
  geom_errorbar(aes(ymin = mean.total.DRR - se, ymax = mean.total.DRR +se), width = .3)+
  labs(x = 'Missing', y = 'Net Biodiversity Effect on Stability', color = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~temp, scales = 'free_x')+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  face = "bold",colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/DeleteMe.png'),width = 9, height =3)

#### Monodata ####

str(MonoData)

duos <- MonoData %>%
  select(-mean, -sd, -se) %>%
  filter(N == '2') %>%
  select(-N)%>%
  rename(obs2Instab = AUC.RR_obs) %>%
  mutate(A2 = ifelse(combination == 'AT', 'T', ifelse(combination == 'AR', 'R', ifelse(combination == 'AG', 'G', ifelse(combination == 'AD', 'D',
                                               ifelse(combination == 'DG', 'G', ifelse(combination == 'DR', 'R', ifelse(combination == 'DT', 'T', 
                                              ifelse(combination == 'GR', 'R', ifelse(combination == 'GT', 'T', ifelse(combination == 'RT', 'T'))))))))) )) %>%
  mutate(A1 = ifelse(combination == 'AT', 'A', ifelse(combination == 'AR', 'A', ifelse(combination == 'AG', 'A', ifelse(combination == 'AD', 'A',
                                               ifelse(combination == 'DG', 'D', ifelse(combination == 'DR', 'D', ifelse(combination == 'DT', 'D', 
                                               ifelse(combination == 'GR', 'G', ifelse(combination == 'GT', 'G', ifelse(combination == 'RT', 'R'))))))))) )) %>%
  group_by(combination, A1, A2, temp) %>%
  summarise(mean2stab = mean(obs2Instab))
unique(duos$A1)
str(MonoData)
names(duos)


monos <- MonoData %>%
  ungroup()%>%
  filter(N %in% c(1)) %>%
  select(-mean, -sd, -se, -label,-combination, -N) %>%
  mutate(spec = speciesID, 
         speciesID = paste(ifelse(speciesID == 'Asterio', 'A', ifelse(speciesID == 'DityCux', 'D', ifelse(speciesID == 'Guido', 'G',
                           ifelse(speciesID == 'Rhizo', 'R', 'T')))))) %>%
  group_by(temp, speciesID, spec) %>%
  summarise(mean.mono = mean(AUC.RR_obs, na.rm = T))


### single df for each species ####
Asterio<- duos %>%
  filter(A1 == 'A'| A2 == 'A')%>%  
  mutate(speciesID = ifelse(A2 != A1, A2, A1),
         present = paste('A')) %>%
  left_join(., monos, by = c('speciesID', 'temp'))  %>%
  mutate(sp2_1 = mean2stab - mean.mono) %>%
  group_by(temp, present) %>%
  summarise(mean.presence = mean(sp2_1),
            sd = sd(sp2_1),
            se.presence = sd/sqrt(n()))
  
Dity<- duos %>%
  filter(A1 == 'D'| A2 == 'D')%>%  
  mutate(speciesID = ifelse(A2 == 'D', A1, A2),
         present = paste('D')) %>%
  left_join(., monos, by = c('speciesID', 'temp')) %>%
  mutate(sp2_1 = mean2stab - mean.mono) %>%
  group_by(temp, present) %>%
  summarise(mean.presence = mean(sp2_1),
            sd = sd(sp2_1),
            se.presence = sd/sqrt(n()))

Guido<- duos %>%
  filter(A1 == 'G'| A2 == 'G')%>%  
  mutate(speciesID = ifelse(A2 == 'G', A1, A2),
         present = paste('G')) %>%
  left_join(., monos, by = c('speciesID', 'temp')) %>%
  mutate(sp2_1 = mean2stab - mean.mono) %>%
  group_by(temp, present) %>%
  summarise(mean.presence = mean(sp2_1),
            sd = sd(sp2_1),
            se.presence = sd/sqrt(n()))

Rhizo<- duos %>%
  filter(A1 == 'R'| A2 == 'R')%>%  
  mutate(speciesID = ifelse(A2 == 'R', A1, A2),
         present = paste('R')) %>%
  left_join(., monos, by = c('speciesID', 'temp')) %>%
  mutate(sp2_1 = mean2stab - mean.mono) %>%
  group_by(temp, present) %>%
  summarise(mean.presence = mean(sp2_1),
            sd = sd(sp2_1),
            se.presence = sd/sqrt(n()))

Thala<- duos %>%
  filter(A1 == 'T'| A2 == 'T')%>%  
  mutate(speciesID = ifelse(A2 == 'T', A1, A2),
         present = paste('T')) %>%
  left_join(., monos, by = c('speciesID', 'temp')) %>%
  mutate(sp2_1 = mean2stab - mean.mono) %>%
  group_by(temp, present) %>%
  summarise(mean.presence = mean(sp2_1),
            sd = sd(sp2_1),
            se.presence = sd/sqrt(n()))

allSpecies <- bind_rows(Asterio, Dity, Guido, Rhizo, Thala)
write.csv(allSpecies, file = here('MicrocosmExp22/Data/allSpecies_MeanPresence2_1.csv'))
