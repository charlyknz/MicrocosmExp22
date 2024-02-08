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

# look at the data 

#### import species-specific biomass data ####
expData<-read.csv('AllRawData_InclBV.csv') %>%
  select(-X)
names(expData)
str(expData)
which(is.na(expData$cells_ml))


#control data
control_stab <- expData %>%
  filter(temp == 'CS') %>%
  select(-no)%>%
  group_by(combination, sampling, species,speciesID ) %>%
  summarise(con.vol = mean(cellVolume, na.rm = T)) 

#treatment data 
treat <- expData %>%
  filter(temp != 'CS') %>%
  select(-no)%>%
  left_join(., control_stab, by = c( 'combination','species', 'sampling', 'speciesID'))


# relative biomass of species in mxiture at t0
RelBV_t0_ <- treat %>%
  filter(sampling == 1) %>%
  group_by(combination, temp, rep) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(combination, temp, speciesID, rep, relBV)
names(RelBV_t0_)


#### NBE on Stability Two species ####

## variables
# M - mono culture of species
names(treat)
mono <- treat %>%
  filter(species == 'mono' ) %>%
  select(-species) %>%
  filter(temp != 'CS') %>%
  rename(Mono_T = cellVolume)%>%
  rename(Mono_C = con.vol)

# Two species Mix 
duo <- treat %>%
  filter(species == 'duo') %>%
  group_by(speciesID,temp, sampling, combination) %>%
  rename(Mix_T = cellVolume) %>%
  rename(Mix_C = con.vol) %>%
  ungroup() 


### merge all df and calculate response ratios ###
all <- mono %>%
  select(-combination)%>%
  right_join(., duo, by = c('temp','speciesID',  'sampling', 'rep')) %>%
  mutate(ratio = Mono_T/Mono_C,
         LRR = log(ratio)) %>%
  group_by(temp, sampling,  combination,rep) %>%
  mutate(sumC = sum(Mix_C,na.rm = T),
         sumT = sum(Mix_T,na.rm = T)) %>%
  ungroup() %>%
  left_join(., RelBV_t0_, by = c("combination", "temp","speciesID", 'rep' )) %>%
  mutate(exp = sumC*ratio*relBV)%>%
  mutate(RR_mono = (Mono_T - Mono_C)/(Mono_T+Mono_C),
         RR_obs = (Mix_T - Mix_C)/(Mix_T+Mix_C),
         RR_ges_obs = (sumT-sumC)/(sumT+sumC),
         RR_exp = (exp-relBV*sumC)/(exp+relBV*sumC)) %>%
  group_by(temp, sampling,  combination,rep) %>%
  mutate(sumexp = sum(exp, na.rm = T),
         RR_ges_exp = (sumexp-sumC)/(sumexp+sumC))%>%
  ungroup()%>%
  mutate(delta = RR_obs - RR_exp) %>%
  mutate(delta_ges = RR_ges_obs-RR_ges_exp) %>%
  group_by(speciesID, temp, sampling,  combination,rep) %>%
  mutate(delta_RR = Mix_T - exp) 

resistance_duo <- filter(all, sampling == 3)
all$delta[is.na(all$delta)]<-0
all$USI <- paste(all$temp, all$combination,all$speciesID, all$rep, sep = '_')

str(all)

all <- all[order(all$sampling),]
ggplot(subset(all, sampling == 1),aes( x = combination, y  = delta_ges))+
  geom_point()+
  facet_wrap(~temp)

##### AUC #####
stab.auc <- tibble()

USI <- unique(all$USI)
for(i in 1:length(USI)){
  temp<-all[all$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.deltaRR<-auc(temp$sampling, temp$delta_RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("linear"),absolutearea = FALSE)
    AUC.RR_obs<-auc(temp$sampling, temp$RR_ges_obs,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.RR_exp<-auc(temp$sampling, temp$RR_ges_exp,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.delta<-auc(temp$sampling, temp$delta,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                   type = c("linear"),absolutearea = FALSE)
    AUC.ges.RR<-auc(temp$sampling, temp$delta_ges,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.RR_mono <- auc(temp$sampling, temp$RR_mono, from = min(temp$sampling, na.rm = T), to = max(temp$sampling, na.rm = T),
                           type = c('linear'), absolutearea = F)
    NBE <- AUC.RR_obs - AUC.RR_exp
    stab.auc<-rbind(stab.auc,
                    tibble(temp[1,c(1:3,7,8,25)],
                               NBE,
                               AUC.deltaRR,                           
                               AUC.delta,
                               AUC.ges.RR,
                               AUC.RR_exp,AUC.RR_obs, AUC.RR_mono))
    rm(temp)
  }
}

summary(stab.auc)


unique(stab.auc$combination)

##### plots#####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")
stab.auc$temp[stab.auc$temp=='fluct'] <- 'Fluctuation'
stab.auc$temp[stab.auc$temp=='inc'] <- 'Increase'
stab.auc$temp[stab.auc$temp=='inc+fluc'] <- 'IncreaseFluctuation'

stab.auc %>%
  group_by(combination,speciesID, temp)%>%
  summarise(mean = mean(AUC.delta,na.rm = T),
            sd = sd(AUC.delta,na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = combination, y =mean, color = speciesID))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .8)+
  geom_point(alpha=0.8, size = 2)+
  facet_grid(~temp, scales = 'free')+
  labs(y = 'Net effect Species Stability')+
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

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/SpeciesStabilityEffect.png'), width = 8,height = 5)



# rel BV at to 
all$temp[all$temp=='fluct'] <- 'Fluctuation'
all$temp[all$temp=='inc'] <- 'Increase'
all$temp[all$temp=='inc+fluc'] <- 'IncreaseFluctuation'

BV_t0_duo <- all %>%
  filter(sampling == 1 ) %>%
  group_by(speciesID, combination,temp) %>%
  summarise(mean.BV = mean(relBV, na.rm = T)) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = mean.BV, fill = speciesID))+
  geom_hline(yintercept = 0.5, color = 'darkgrey', linewidth = 0.2)+
  scale_fill_manual(values=cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  labs(y =  expression(Mean~Relative~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Combination', color = 'Species')+
  theme_bw()+
  theme(legend.position = 'none')
BV_t0_duo
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/relBV_t0_TwoSpecies.png'), width = 8, height = 4)


#### NBE on Stability: 4-species and 5-species Mix ####

# M - monoculture of species - use the one from above
str(mono)

# Mix 
names(treat)
mix <- treat %>%
  filter(species %in% c('MIX', 'quattro')) %>%
  ungroup() 

mix_con <- control_stab %>%
  filter(species %in% c('MIX', 'quattro')) %>%
  rename(Mix_C = con.vol) %>%
  ungroup()
str(mix_con)

mix_all <- mix %>%
  rename(Mix_T = cellVolume) %>%
  left_join(., mix_con, by = c("speciesID", "species","sampling","combination"))


### merge all df and calculate response ratios ###

relBVt0 <- treat %>%
  filter(species %in% c('MIX', 'quattro')) %>%
  filter(sampling == 1) %>%
  group_by( combination, temp, rep) %>%
  mutate(sum = sum(cellVolume, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellVolume/sum) %>%
  select(combination, temp, speciesID, rep, relBV)

allMix <- mono %>%
  select(-combination)%>%
  right_join(., mix_all, by = c('temp','speciesID', 'sampling', 'rep')) %>%
  mutate(ratio = Mono_T/Mono_C,
         LRR = log(ratio)) %>%
  group_by(temp, sampling, combination,rep) %>%
  mutate(sumC = sum(Mix_C,na.rm = T),
         sumT = sum(Mix_T,na.rm = T)) %>%
  ungroup() %>%
  left_join(., relBVt0, by = c("combination", "temp","speciesID", 'rep')) %>%
  mutate(exp = sumC*ratio*relBV)%>%
  mutate(RR_obs = (Mix_T - Mix_C)/(Mix_T+Mix_C),
         RR_ges_obs = (sumT-sumC)/(sumT+sumC),
         RR_exp = (exp-relBV*sumC)/(exp+relBV*sumC)) %>%
  group_by(temp, sampling,  combination,rep) %>%
  mutate(sumexp = sum(exp, na.rm = T),
         RR_ges_exp = (sumexp-sumC)/(sumexp+sumC))%>%
  ungroup()%>%
  mutate(delta_ges = RR_ges_obs-RR_ges_exp) %>%
  group_by(speciesID, temp, sampling, combination,rep) %>%
  mutate(delta_RR = Mix_T - exp)

resistance_mix <- filter(allMix,sampling == 3)
allMix$USI <- paste(allMix$temp, allMix$combination,allMix$speciesID, allMix$rep, sep = '_')

##### AUC #####
stab.auc.mix <- tibble()
names(allMix)
USI <- unique(allMix$USI)
for(i in 1:length(USI)){
  temp<-allMix[allMix$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.deltaRR<-auc(temp$sampling, temp$delta_RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("linear"),absolutearea = FALSE)
        AUC.ges.RR<-auc(temp$sampling, temp$delta_ges,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.RR_exp<-auc(temp$sampling, temp$RR_ges_exp,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.RR_obs<-auc(temp$sampling, temp$RR_ges_obs,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    NBE <- AUC.RR_obs-AUC.RR_exp
    stab.auc.mix<-rbind(stab.auc.mix,
                        data.frame(temp[1,c(1:3,7,8,24)],
                                   NBE,AUC.deltaRR,
                                   AUC.ges.RR,AUC.RR_obs,AUC.RR_exp))
    rm(temp)
  }
}

summary(stab.auc.mix)
unique(stab.auc.mix$combination)

##### plots#####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")
stab.auc.mix$temp[stab.auc.mix$temp=='fluct'] <- 'Fluctuation'
stab.auc.mix$temp[stab.auc.mix$temp=='inc'] <- 'Increase'
stab.auc.mix$temp[stab.auc.mix$temp=='inc+fluc'] <- 'IncreaseFluctuation'


#pannels for temp only
stab.auc.mix %>%
  group_by(combination, speciesID, temp)%>%
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

#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/deltaRR_quattroSpec_chaos.png'), width = 8,height = 5)


# rel BV at to 
allMix$combination[allMix$combination == 'ADGRT'] <- 'Mix'

allMix$temp[allMix$temp=='fluct'] <- 'Fluctuation'
allMix$temp[allMix$temp=='inc'] <- 'Increase'
allMix$temp[allMix$temp=='inc+fluc'] <- 'IncreaseFluctuation'


BV_t0_mix <- allMix %>%
  filter(sampling == 1 ) %>%
  group_by(speciesID, combination,temp) %>%
  summarise(mean.BV = mean(relBV, na.rm = T)) %>%
  ggplot(. )+
  geom_col(aes(x = combination, y = mean.BV, fill = speciesID))+
  scale_fill_manual(values=cbbPalette)+
  facet_grid(~temp, scales = 'free')+
  labs(y =  expression(Mean~Relative~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Combination', color = 'Species')+
  theme_bw()+
  theme(legend.position = 'none')
BV_t0_mix
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/relBV_t0_45species.png'), width = 9, height = 4)

plot_grid(BV_t0_duo, BV_t0_mix, labels = c('(a)', '(b)'), ncol = 1)
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/RelativeBVatT0.png'), width = 8, height = 8)



#### Merge Stabilities ####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )

d1 <- stab.auc %>%
  select(combination, speciesID, rep,temp,NBE, AUC.RR_exp,AUC.RR_obs, AUC.RR_mono)
d2 <- stab.auc.mix %>%
  select(combination, speciesID, rep,temp,NBE,AUC.RR_exp,AUC.RR_obs)
d2$combination[d2$combination=='Mix']<-'ADGRT'

d3 <- d1 %>%
  select(-AUC.RR_mono)%>%
  bind_rows(., d2) %>%
  mutate(N = str_length(combination)  )  
which(is.na(d3))

#write.csv(d3, file = here('MicrocosmExp22/Data/netEffect.csv'))


#### START plots####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )

#d3 <- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/netEffect.csv')


#### NBES plot####

p1<-d3%>%
  group_by(N, temp,combination) %>%
  summarise(mean.total.DRR = mean(NBE, na.rm = T),
            sd = sd(NBE,na.rm = T),
            se = sd/sqrt(n()))%>%
  mutate(label = paste(ifelse( N == 2, '2 species', ifelse(N == 4, '4 species', '5 species'))))%>%
  ggplot(., aes(x = combination, y = mean.total.DRR,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 0.9)+
  geom_errorbar(aes(ymin = mean.total.DRR - se, ymax = mean.total.DRR +se), width = .3, alpha = .7)+
  labs(x = 'Combination', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~label, scales = 'free')+
  theme_bw()+
  theme(legend.position = 'right',
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

p2 <- d3%>%
  group_by(N, temp) %>%
  summarise(mean.total.DRR = mean(NBE, na.rm = T),
         sd = sd(NBE,na.rm = T),
         se = sd/sqrt(n()))%>%
  ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(x = as.factor(N), y = mean.total.DRR, color=temp, shape=temp),size = 3.5, alpha = 0.9)+
  geom_errorbar(aes(x=as.factor(N),ymin = mean.total.DRR - se, ymax = mean.total.DRR +se, col=temp), alpha = 0.7, width = .1)+
  labs(x = 'Species Richness', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
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


nbes <- plot_grid( p2,p1+theme(legend.position = 'none'),legendb,hjust = -0.05, labels = c('(a)', '(b)'), ncol = 3,rel_widths = c( 2/7,4/7,1/7), rel_heights = c(10,0.2))
#ggsave(plot = nbes, file = here('MicrocosmExp22/output/Fig_NBES.png'), width = 14, height = 4.5)


#### Merge NBES with NBE on Functioning ####

#run 05_NBE_HectorLoreau_NetBiodivEffect through
source(here("~/Desktop/Exp22/MicrocosmExp22/code/05_NBE_HectorLoreau_NetBiodivEffect.R"))

(nbes/nbef)
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/Fig3_NBE_NBES.png'), width = 14, height = 9)



#### Figure AUC.RR ####
names(d1)
MonoData <- d1%>%
  distinct( speciesID, temp,AUC.RR_mono) %>%
  rename(AUC.RR_obs= AUC.RR_mono)%>%
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
  labs(x = 'Combination', y = 'Observed OEV', color = 'Treatment', shape = 'Treatment')+
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
  labs(x = 'Combination', y = '', color = 'Treatment', shape = 'Treatment')+
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
  labs(x = 'Combination', y = 'Observed OEV', color = 'Treatment', shape = 'Treatment')+
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
  labs(x = 'Combination', y = '', color = 'Treatment', shape = 'Treatment')+
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

RR1+RR2+RR4+RR5+
  plot_layout(widths = c( 3,3), nrow = 2)+ plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a", tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')')

ggsave(plot = last_plot(), file = here('/Users/charlottekunze/Desktop/Exp22/MicrocosmExp22/output/Fig2_ObservedStab_patchwork.png'),width = 10, height = 8)


  
#### NBES resistance ####
resistance_duo 
resistance_mix

all_resistance <- rbind(resistance_duo, resistance_mix)

all_resistance$N <- NA
all_resistance$N[all_resistance$species == 'mono']<-'1'
all_resistance$N[all_resistance$species == 'duo']<-'2'
all_resistance$N[all_resistance$species == 'quattro']<-'4'
all_resistance$N[all_resistance$species == 'MIX']<-'5'

all_resistance$temp[all_resistance$temp=='fluct'] <- 'Fluctuation'
all_resistance$temp[all_resistance$temp=='inc'] <- 'Increase'
all_resistance$temp[all_resistance$temp=='inc+fluc'] <- 'IncreaseFluctuation'

NBESresistance <- all_resistance  %>%
  group_by(temp, N) %>%
  mutate(Mean = mean(delta_ges),
         sd = sd(delta_ges),
         se= sd/sqrt(n()))%>%
  ggplot(.) +
    geom_hline(yintercept = 0)+
    geom_point( aes( x = N, y = delta_ges),color='black', alpha = 0.2)+
  geom_errorbar(aes(x= N, y = Mean, ymin = Mean-se, ymax = Mean+se,color = temp), width = 0.3, alpha = 0.8)+
  geom_point(aes( x = N, y = Mean, color = temp, shape = temp),  size = 3)+
    facet_wrap(~temp)+
  labs(x = 'Species Richness', y=expression(NBES[CV]))+
  scale_colour_brewer(palette = "Set1")+
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
NBESresistance
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NBES_resistance.png'), width = 8, height = 4)


NBESres1 <- all_resistance  %>%
  group_by(temp, N) %>%
  summarise(Mean = mean(delta_ges),
         sd = sd(delta_ges),
         se= sd/sqrt(n()))%>%
  ggplot(.) +
  geom_hline(yintercept = 0)+
  #geom_point( aes( x = N, y = delta_ges), alpha = 0.3)+
  geom_point(aes( x = N, y = Mean, color = temp, shape = temp), size = 3)+
  geom_errorbar(aes(x= N, y = Mean, ymin = Mean-se, ymax = Mean+se, color = temp), alpha=0.8,width = 0.3)+
  #facet_wrap(~temp)+
  labs(x = ' ', y=expression(NBES[resistance]))+
  scale_colour_brewer(palette = "Set1")+
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
NBESres1

labeln <- c('2'='2 species', '4'='4 species','5'='5 species')
NBESres2 <- all_resistance  %>%
  group_by(temp, N,combination) %>%
  summarise(Mean = mean(delta_ges),
         sd = sd(delta_ges),
         se= sd/sqrt(n()))%>%
  ggplot(.) +
  geom_hline(yintercept = 0)+
  #geom_point( aes( x = N, y = delta_ges), alpha = 0.3)+
  geom_point(aes( x = combination, y = Mean, color = temp, shape = temp), size = 3)+
  geom_errorbar(aes(x= combination, y = Mean, ymin = Mean-se, ymax = Mean+se, color = temp), alpha=0.8,width = 0.3)+
  facet_wrap(~N, labeller = labeller(N = labeln), scales = 'free_x')+
  labs(x = ' ', y=expression(NBES[resistance]), color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  theme_bw()+
  theme(legend.position = 'right',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
NBESres2

legend_nb<-get_legend(NBESres2)
PlotNBES <- plot_grid( NBESres1,NBESres2+theme(legend.position = 'none'),hjust = -0.05, labels = c('(a)', '(b)'), ncol = 3, rel_widths = c(2/7,5/7))
PlotNBES
ggsave(plot = PlotNBES, file = here('MicrocosmExp22/output/Fig_NBESresistance.png'), width = 14, height = 4)


#### NBES CV ####
## calculate how much variability is contributed by diversity 
# therefore we assess the CV observed and expected from monoculture
str(all)
CV_duo <- all %>%
  ungroup()%>%
  distinct(temp, rep, sampling, species,combination,RR_ges_obs,RR_ges_exp ) %>%
  group_by(temp, combination, species,rep) %>%
  summarise(MeanObs = mean(RR_ges_obs),
            MeanExp = mean(RR_ges_exp),
            sdObs = sd(RR_ges_obs),
            sdExp = sd(RR_ges_exp),
            CVobs = MeanObs/sdObs,
            CVExp = MeanExp/sdExp,
            NBES_CV = CVobs-CVExp)%>%
  ungroup() %>%
  select(-MeanObs,-MeanExp,-sdObs,-sdExp)

CV_mix <- allMix %>%
  ungroup()%>%
  distinct(temp, rep, sampling, species,combination,rep,RR_ges_obs,RR_ges_exp ) %>%
  group_by(temp, combination, species,rep) %>%
  summarise(MeanObs = mean(RR_ges_obs),
            MeanExp = mean(RR_ges_exp),
            sdObs = sd(RR_ges_obs),
            sdExp = sd(RR_ges_exp),
            CVobs = MeanObs/sdObs,
            CVExp = MeanExp/sdExp,
            NBES_CV = CVobs-CVExp)%>%
  ungroup() %>%
  select(-MeanObs,-MeanExp,-sdObs,-sdExp)

CV <- bind_rows(CV_duo, CV_mix)

CV$N <- NA
CV$N[CV$species == 'mono']<-'1'
CV$N[CV$species == 'duo']<-'2'
CV$N[CV$species == 'quattro']<-'4'
CV$N[CV$species == 'MIX']<-'5'

CV$combination[CV$combination == 'Mix']<-'ADGRT'

CV$temp[CV$temp=='fluct'] <- 'Fluctuation'
CV$temp[CV$temp=='inc'] <- 'Increase'
CV$temp[CV$temp=='inc+fluc'] <- 'IncreaseFluctuation'

NBESCV <- CV %>%
  group_by(temp, N) %>%
  mutate(mean_NBES_CV = mean(NBES_CV),
         sd_NBES_CV = sd(NBES_CV),
         se_NBES_CV=sd_NBES_CV/sqrt(n())) %>%
ggplot(.) +
  geom_hline(yintercept = 0)+
  geom_point( aes( x = N, y = NBES_CV),color = 'black',  alpha = 0.3)+
  geom_point(aes( x = N, y = mean_NBES_CV,color = temp, shape = temp), size = 3)+
  geom_errorbar(aes(x= N, y = mean_NBES_CV, ymin = mean_NBES_CV-se_NBES_CV, ymax = mean_NBES_CV+se_NBES_CV,color = temp),width = 0.3, alpha = 0.8)+
  facet_wrap(~temp)+
  labs(x = 'Species Richness', y=expression(NBES[CV]))+
  scale_colour_brewer(palette = "Set1")+
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
NBESCV
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/NBES_CV.png'), width = 8, height = 4)



plot_grid( NBESresistance,NBESCV,labels = c('(a)', '(b)'), ncol = 1)
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/Fig_NBESmetrics_overall.png'), width = 8, height = 6)


# 
NBEScv1 <- CV  %>%
  group_by(temp, N) %>%
  summarise(mean_NBES_CV = mean(NBES_CV),
         sd_NBES_CV = sd(NBES_CV),
         se_NBES_CV=sd_NBES_CV/sqrt(n())) %>%
  ggplot(.) +
  geom_hline(yintercept = 0)+
  geom_point(aes( x = N, y = mean_NBES_CV, color = temp, shape = temp), size = 3)+
  geom_errorbar(aes(x= N, y = mean_NBES_CV, ymin = mean_NBES_CV-se_NBES_CV, ymax = mean_NBES_CV+se_NBES_CV, color = temp), alpha=0.8,width = 0.3)+
  #facet_wrap(~temp)+
  labs(x = ' ', y=expression(NBES[CV]),color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  theme_bw()+
  theme(legend.position = 'right',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
NBEScv1

legendcv<-get_legend(NBEScv1)

NBEScv2 <- CV  %>%
  group_by(temp, N, combination) %>%
  summarise(mean_NBES_CV = mean(NBES_CV),
            sd_NBES_CV = sd(NBES_CV),
            se_NBES_CV=sd_NBES_CV/sqrt(n())) %>%
  ggplot(.) +
  geom_hline(yintercept = 0)+
  geom_point(aes( x = combination, y = mean_NBES_CV, color = temp, shape = temp), size = 3)+
  geom_errorbar(aes(x= combination, y = mean_NBES_CV, ymin = mean_NBES_CV-se_NBES_CV, ymax = mean_NBES_CV+se_NBES_CV, color = temp), alpha=0.8,width = 0.3)+
  facet_wrap(~N,labeller = labeller(N = labeln), scales = 'free_x')+
  labs(x = ' ', y=expression(NBES[CV]), color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
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
NBEScv2



PlotCV <- plot_grid( NBEScv1+theme(legend.position = 'none'),NBEScv2+theme(legend.position = 'none'),hjust = -0.05, labels = c('(c)', '(d)'), ncol = 3, rel_widths = c(2/7,5/7))
PlotCV

PlotNBES/PlotCV
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/Fig_NBESmetrics.png'), width = 14, height = 6)

#### Explorative plots and analysis ####

### LRR Resistance ###
         
resistance <- treat %>%
  filter(sampling ==3 ) %>%
  group_by(temp, combination, species,rep) %>%
  summarise(treat.tot = sum(cellVolume),
            con.tot = sum(con.vol)) %>%
  mutate(LRR = log(treat.tot/con.tot))

resistance$species <- factor(as.factor(resistance$species),levels=c("mono","duo",   "quattro" ,  "MIX"))

resistance$N <- NA
resistance$N[resistance$species == 'mono']<-'1'
resistance$N[resistance$species == 'duo']<-'2'
resistance$N[resistance$species == 'quattro']<-'4'
resistance$N[resistance$species == 'MIX']<-'5'

resistance$temp[resistance$temp=='fluct'] <- 'Fluctuation'
resistance$temp[resistance$temp=='inc'] <- 'Increase'
resistance$temp[resistance$temp=='inc+fluc'] <- 'IncreaseFluctuation'

ggplot(resistance, aes( x = N, y = LRR))+
  geom_point(alpha = 0.5)+
  geom_boxplot(alpha = 0.5)+
  facet_wrap(~temp)+
  theme_bw()
#ggsave(plot = last_plot(), file = here('~/Desktop/Exp22/MicrocosmExp22/output/LRR.png'), width = 8, height = 4)
# 21 unique combinations x 3 treatments x 3 rep 

### RR over time ###
RR_merge <- rbind(allMix, all) %>%
  distinct(temp, rep, sampling, combination, RR_ges_obs, RR_ges_exp, delta_ges) %>%
  group_by(temp, sampling, combination)%>%
  summarise(mean.delta.RR = mean(delta_ges, na.rm = T),
            sd = mean(delta_ges),
            se = sd/sqrt(n())) %>%
  mutate(S = str_length(combination)) %>%
  mutate(S = paste(ifelse(S == 3, 5, S))) %>%
  arrange(desc(S))
RR_merge$combination <- factor(as.factor(RR_merge$combination) , 
                                  levels = c("AD" ,"AG" ,"AR","AT" ,"DG","DR","DT", "GR" , "GT" ,"RT",
                                             "ADGR" ,"ADGT", "ADRT" ,"AGRT", "DGRT" ,"Mix"))
ggplot(RR_merge, aes(x = sampling, y = mean.delta.RR, color = temp, shape = temp))+
  geom_hline(yintercept = 0, color ='black')+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean.delta.RR-se, ymax = mean.delta.RR + se), width = .8)+
  geom_point()+
  facet_wrap(~combination, ncol = 5)+
  theme_bw()

##plot 2
RR_merge1 <- rbind(allMix, all) %>%
  distinct(temp, rep, sampling, combination, RR_ges_obs, RR_ges_exp, delta_ges) %>%
  mutate(S = str_length(combination)) %>%
  mutate(S = paste(ifelse(S == 3, 5, S))) %>%
  arrange(desc(S))%>%
  group_by(temp, sampling,S)%>%
  summarise(mean.delta.RR = mean(delta_ges, na.rm = T),
            sd = mean(delta_ges),
            se = sd/sqrt(n())) 

ggplot(RR_merge1, aes(x = sampling, y = mean.delta.RR, color = temp, shape = temp))+
  geom_hline(yintercept = 0, color ='black')+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean.delta.RR-se, ymax = mean.delta.RR + se), width = .8)+
  geom_point()+
  facet_wrap(~S)+
  theme_bw()


### Plots Missing species ###

d3$combination[d3$combination == 'ADGR'] <- 'T'
d3$combination[d3$combination == 'ADGT'] <- 'R'
d3$combination[d3$combination == 'DGRT'] <- 'A'
d3$combination[d3$combination == 'AGRT'] <- 'D'
d3$combination[d3$combination == 'ADRT'] <- 'G'
d3$combination[d3$combination == 'none'] <- 'Xnone'

d3%>%
  group_by(N, temp,combination) %>%
  summarise(mean.total.DRR = mean(NBE, na.rm = T),
            sd = sd(NBE,na.rm = T),
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
#ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/DeleteMe.png'),width = 9, height =3)


### SGC data ###

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


### single df for each species ###
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
#write.csv(allSpecies, file = here('MicrocosmExp22/Data/allSpecies_MeanPresence2_1.csv'))
