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


# look at the data 

#### import species-specific biomass data ####
expData<-read.csv('Data/AllRawData_InclBV.csv') %>%
  select(-X) %>%
  mutate(cellV_mm_ml = cellVolume/10^9) %>%
  select(-cellVolume)

names(expData)
str(expData)
which(is.na(expData$cells_ml))


#control data
control_stab <- expData %>%
  filter(temp == 'CS') %>%
  select(-no)%>%
  group_by(combination, sampling, species,speciesID ) %>%
  summarise(con.vol = mean(cellV_mm_ml, na.rm = T)) 

#treatment data 
treat <- expData %>%
  filter(temp != 'CS') %>%
  select(-no)%>%
  left_join(., control_stab, by = c( 'combination','species', 'sampling', 'speciesID'))


# relative biomass of species in mxiture at t0
RelBV_t0_ <- treat %>%
  filter(sampling == 1) %>%
  group_by(combination, temp, rep) %>%
  mutate(sum = sum(cellV_mm_ml, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellV_mm_ml/sum) %>%
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
  rename(Mono_T = cellV_mm_ml)%>%
  rename(Mono_C = con.vol)

# Two species Mix 
duo <- treat %>%
  filter(species == 'duo') %>%
  group_by(speciesID,temp, sampling, combination) %>%
  rename(Mix_T = cellV_mm_ml) %>%
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
names(all)

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



# relative BV at to 
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
  rename(Mix_T = cellV_mm_ml) %>%
  left_join(., mix_con, by = c("speciesID", "species","sampling","combination"))


### merge all df and calculate response ratios ###

relBVt0 <- treat %>%
  filter(species %in% c('MIX', 'quattro')) %>%
  filter(sampling == 1) %>%
  group_by( combination, temp, rep) %>%
  mutate(sum = sum(cellV_mm_ml, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellV_mm_ml/sum) %>%
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


# relative BV at to 
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

cowplot::plot_grid(BV_t0_duo, BV_t0_mix, labels = c('(a)', '(b)'), ncol = 1)
ggsave(plot = last_plot(), file = here('output/RelativeBVatT0.png'), width = 8, height = 8)



#### Merge Stabilities ####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )

d1 <- stab.auc %>%
  distinct(combination, rep,temp,NBE, AUC.RR_exp,AUC.RR_obs)
d2 <- stab.auc.mix %>%
  distinct(combination, rep,temp,NBE,AUC.RR_exp,AUC.RR_obs)
d2$combination[d2$combination=='Mix']<-'ADGRT'

d3 <- d1 %>%
  bind_rows(., d2) %>%
  mutate(N = str_length(combination)  )  
which(is.na(d3))

## save for R-Stats
write.csv(d3, file = here('Data/NBES.csv'))


#### START plots####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )


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
  labs(x = 'Species Combination', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~label, scales = 'free')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
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


nbes <- cowplot::plot_grid( p2,p1+theme(legend.position = 'none'),legendb,hjust = -0.05, labels = c('(a)', '(b)'), ncol = 3,rel_widths = c( 2/7,4/7,1/7), rel_heights = c(10,0.2))


#### Merge NBES with NBE on Functioning ####

#run 05_NBE_HectorLoreau_NetBiodivEffect through
source(here("code/05_NBE_HectorLoreau_NetBiodivEffect.R"))

(nbes/nbef)
ggsave(plot = last_plot(), file = here('output/Fig3_NBE_NBES.png'), width = 15, height = 9)



#### Figure AUC.RR ####
names(d1)
MonoData <- stab.auc %>%
  ungroup()%>%
  distinct( speciesID, temp,rep,AUC.RR_mono) %>%
  rename(AUC.RR_obs= AUC.RR_mono)%>%
  mutate(combination = paste(ifelse(speciesID == 'Asterio', 'A', ifelse(speciesID == 'DityCux', 'D', ifelse(speciesID == 'Guido', 'G', ifelse(speciesID == 'Rhizo', 'R', 'T')))))) %>%
  select(-speciesID) %>%
  bind_rows(., select(d1,combination,temp, rep,AUC.RR_obs)) %>%
  bind_rows(., select(d2 ,combination,temp, rep,AUC.RR_obs)) %>%
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
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
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
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
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
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
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
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
RR5

RR1+RR2+RR4+RR5+
  plot_layout(widths = c( 3,3), nrow = 2)+ plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a", tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')')

ggsave(plot = last_plot(), file = here('output/Fig2_ObservedStab_patchwork.png'),width = 11, height = 8)


  
#### NBES resistance ####
resistance_duo 
resistance_mix

all_resistance <- rbind(resistance_duo, resistance_mix) %>%
  ungroup()%>%
  distinct(temp, rep, delta_ges, species)

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
  labs(x = 'Species Richness', y=expression(NBES[resistance]))+
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

cowplot::plot_grid( NBESresistance,NBESCV,labels = c('(a)', '(b)'), ncol = 1)
ggsave(plot = last_plot(), file = here('output/FigS1_NBESmetrics_overall.png'), width = 8, height = 6)


