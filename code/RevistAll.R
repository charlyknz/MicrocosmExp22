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


#### import species-specific biomass data ####
expData<-read.csv('Data/AllRawData_InclBV.csv') %>%
  select(-X) %>%
  mutate(cellV_mm_ml = cellVolume/10^9) %>%
  select(-cellVolume)

#check import
names(expData)
str(expData)
which(is.na(expData$cells_ml))


#control data
control_stab <- expData %>%
  filter(temp == 'CS') %>%
  select(-no)%>%
  group_by(combination, sampling, species,speciesID ) %>%
  reframe(con.vol = mean(cellV_mm_ml, na.rm = T)) 

#treatment & control data 
treat <- expData %>%
  filter(temp != 'CS') %>%
  select(-no)%>%
  left_join(., control_stab)


# relative biomass of species in mixture at t0
RelBV_t0_ <- treat %>%
  filter(sampling == 1) %>%
  group_by(combination, temp, rep) %>%
  mutate(sum = sum(cellV_mm_ml, na.rm = T)) %>%
  ungroup() %>%
  mutate(relBV = cellV_mm_ml/sum) %>%
  select(combination, temp, speciesID, rep, relBV)
names(RelBV_t0_)


#### Net Biodiversity Effect on Stability (NBES)####

## variables
# M - mono culture of species
names(treat)
mono <- treat %>%
  filter(species == 'mono' ) %>%
  select(-species) %>%
  filter(temp != 'CS') %>%
 # group_by(temp, speciesID, combination, sampling)%>%
  rename(Mono_T = cellV_mm_ml,
         Mono_C = con.vol)

# Two species Mix 
Mixtures <- treat %>%
  filter(species != 'mono')%>% 
  rename(Mix_T = cellV_mm_ml) %>%
  rename(Mix_C = con.vol) 

### merge all df and calculate response ratios ###
alle <- mono %>%
  select(-combination)%>%
  right_join(., Mixtures) %>%
  mutate(ratio = Mono_T/Mono_C,
         LRR = log(ratio)) %>%
  group_by(temp, sampling, species, combination, rep) %>%
  mutate(sumT = sum(Mix_T,na.rm = T),
         sumC = sum(Mix_C,na.rm = T)) %>%
  ungroup() %>%
  left_join(., RelBV_t0_) %>%
  mutate(exp = sumC*ratio*relBV)%>%
  group_by(temp, sampling, species,  combination,rep) %>%
  mutate(sumexp = sum(exp, na.rm = T))%>%
  group_by(temp, sampling, species,  combination)%>%
  mutate(mean_sumexp = mean(sumexp, na.rm = T))%>%
  ungroup()%>%
  mutate(RR_ges_exp = (mean_sumexp-sumC)/(mean_sumexp+sumC),
         RR_ges_obs = (sumT-sumC)/(sumT+sumC),
         RR_mono = (Mono_T - Mono_C)/(Mono_T+Mono_C)) %>%
  ungroup() %>%
  mutate(delta_ges = RR_ges_obs-RR_ges_exp) 

#subset for resistance calculation
resistance <- alle %>%
  filter(sampling == 3) %>%
  distinct(species, combination,temp, rep, delta_ges)


#create Unique identifier USI
alle$USI <- paste(alle$species,alle$temp, alle$combination,alle$speciesID, alle$rep, sep = '_')
str(alle)

#order after timepoints
alle <- alle[order(alle$sampling),]

ggplot(subset(alle, sampling == 1),aes( x = combination, y  = delta_ges))+
  geom_point()+
  facet_wrap(~temp)

##### Area Under the Curve - AUC #####
stab.auc <- tibble() #empty tibble
names(alle)

USI <- unique(alle$USI) ##create Unique identifier USI

for(i in 1:length(USI)){
  temp<-alle[alle$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR_obs<-auc(temp$sampling, temp$RR_ges_obs,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.RR_exp<-auc(temp$sampling, temp$RR_ges_exp,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.deltaRR<-auc(temp$sampling, temp$delta_ges,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                    type = c("linear"),absolutearea = FALSE)
    AUC.RR_mono<-auc(temp$sampling, temp$RR_mono,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                     type = c("linear"),absolutearea = FALSE)
    NBE <- AUC.RR_obs - AUC.RR_exp
    stab.auc<-rbind(stab.auc,
                    tibble(temp[1,c(1:3,7,8,22)],
                           NBE,
                           AUC.deltaRR,
                           AUC.RR_exp,
                           AUC.RR_obs,
                           AUC.RR_mono))
    rm(temp)
  }
}

#### START plots####
tempPalette <- c('black',"#E41A1C" ,"#377EB8" ,"#4DAF4A" )


#### NBES plot####
stab.auc$combination[stab.auc$combination=='Mix']<-'ADGRT'

stab.auc$temp[stab.auc$temp=='fluct'] <- 'Fluctuation'
stab.auc$temp[stab.auc$temp=='inc'] <- 'Increase'
stab.auc$temp[stab.auc$temp=='inc+fluc'] <- 'Increase + Fluctuation'

d3 <- stab.auc %>%
  mutate(N = str_length(combination)  )  

p1<-d3%>%
  group_by(N, temp,combination, species) %>%
  reframe(mean.total.DRR = mean(NBE, na.rm = T),
            sd = sd(NBE,na.rm = T),
            se = sd/sqrt(n()))%>%
  mutate(label = paste(ifelse( N == 2, '2 species', ifelse(N == 4, '4 species', '5 species'))))%>%
  filter(N!= 5)%>%
  ggplot(., aes(x = combination, y = mean.total.DRR,color = temp, shape = temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(size = 2.5, alpha = 1)+
  geom_errorbar(aes(ymin = mean.total.DRR - se, ymax = mean.total.DRR +se), width = .3, alpha = .7)+
  labs(x = 'Species Combination', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
  # scale_x_continuous(limits = c(0.5,5.5), breaks = c(2,4,5))+
  scale_colour_brewer(palette = "Set1")+
  facet_grid(~label, scales = 'free')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 14))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))
p1
legendb<-get_legend(p1)

p2 <- d3%>%
  group_by(N,species, temp) %>%
  reframe(mean.total.DRR = mean(NBE, na.rm = T),
            sd = sd(NBE,na.rm = T),
            se = sd/sqrt(n()))%>%
  ggplot(.)+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(x = as.factor(N), y = mean.total.DRR, color=temp, shape=temp),size = 3.5, alpha = 0.9)+
  geom_errorbar(aes(x=as.factor(N),ymin = mean.total.DRR - se, ymax = mean.total.DRR +se, col=temp), alpha = 0.7, width = .1)+
  labs(x = 'Species Richness', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
p2

nbes <- cowplot::plot_grid( p2,p1+theme(legend.position = 'none'),
                            legendb,
                            hjust = -1.1, 
                            labels = c('(a)', '(b)'), 
                            ncol = 3,
                            rel_widths = c( 2/7,4/7,1/7), 
                            rel_heights = c(10,0.2))
nbes

p3 <- d3%>%
  group_by(N,species, temp) %>%
  mutate(mean.total.DRR = mean(NBE, na.rm = T),
         sd = sd(NBE,na.rm = T),
         se = sd/sqrt(n()))%>%
  ggplot(.,aes(x = as.factor(N), y = mean.total.DRR, color=temp, shape=temp))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_point(aes(y=NBE),color ='black',size = 1.5, alpha = 0.2)+
  geom_errorbar(aes(x=as.factor(N),ymin = mean.total.DRR - se, ymax = mean.total.DRR +se, col=temp), alpha = 0.7, width = .1)+
  geom_point(size = 3.5, alpha = 0.9)+
  labs(x = 'Species Richness', y = 'Net Biodiversity Effect on Stability', color = 'Treatment', shape = 'Treatment')+
  scale_colour_brewer(palette = "Set1")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))+
  facet_grid(~temp)+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
p3
ggsave(p3, file = here('output/NBES_raw.png'), width = 8, height = 4)

write.csv(d3, file = here('Data/NBES_revisited.csv'))



#### Merge NBES with NBE on Functioning ####

#run 05_NBE_HectorLoreau_NetBiodivEffect through
source(here("code/04_NBE_HectorLoreau_NetBiodivEffect.R"))

(nbes/nbef)
ggsave(plot = last_plot(), file = here('output/Figure5_NBE_NBES.tiff'), width = 15, height = 9)



#### Figure Observed Stability ####
names(d3)
MonoData <- d3 %>%
  select(speciesID, temp,rep,AUC.RR_mono)%>%
  unique()%>%
  rename(AUC.RR_obs= AUC.RR_mono)%>%
  ungroup()%>%
  mutate(N = as.numeric(paste('1')),
         combination = paste(ifelse(speciesID == 'Asterio', 'A', ifelse(speciesID == 'DityCux', 'D', 
                             ifelse(speciesID == 'Guido', 'G', ifelse(speciesID == 'Rhizo', 'R', 'T')))))) %>%
  bind_rows(., select(d3,speciesID,combination,N,temp, rep,AUC.RR_obs)) %>%
 # mutate(combination = paste(ifelse(speciesID == 'Asterio', 'A', ifelse(speciesID == 'DityCux', 'D', ifelse(speciesID == 'Guido', 'G', ifelse(speciesID == 'Rhizo', 'R', 'T')))))) %>%
  select(-speciesID) %>%
   mutate(label = paste(ifelse(N == 1, 'Monoculture', ifelse( N == 2, '2 species', ifelse(N == 4, '4 species', '5 species')))))%>%
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
  scale_y_continuous(limits = c(-7.5,7.5))+
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
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))
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
  scale_y_continuous(limits = c(-7.5,7.5))+
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
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))
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
  scale_y_continuous(limits = c(-7.5,7.5))+
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
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))
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
  scale_y_continuous(limits = c(-7.5,7.5))+
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
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))
RR5

RR1+RR2+RR4+RR5+
  plot_layout(widths = c( 3,3), nrow = 2)+ plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a", tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')')

ggsave(plot = last_plot(), file = here('output/Figure4_ObservedStab.tiff'),width = 11, height = 8)



#### NBES resistance ####
names(resistance)

resistance$N <- NA
resistance$N[resistance$species == 'mono']<-'1'
resistance$N[resistance$species == 'duo']<-'2'
resistance$N[resistance$species == 'quattro']<-'4'
resistance$N[resistance$species == 'MIX']<-'5'

resistance$temp[resistance$temp=='fluct'] <- 'Fluctuation'
resistance$temp[resistance$temp=='inc'] <- 'Increase'
resistance$temp[resistance$temp=='inc+fluc'] <- 'Increase + Fluctuation'


NBESresistance <- resistance  %>%
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
CV <- alle %>%
  distinct(temp, rep, sampling, species,combination, RR_ges_obs,RR_ges_exp ) %>%
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

CV$N <- NA
CV$N[CV$species == 'mono']<-'1'
CV$N[CV$species == 'duo']<-'2'
CV$N[CV$species == 'quattro']<-'4'
CV$N[CV$species == 'MIX']<-'5'

CV$combination[CV$combination == 'Mix']<-'ADGRT'

CV$temp[CV$temp=='fluct'] <- 'Fluctuation'
CV$temp[CV$temp=='inc'] <- 'Increase'
CV$temp[CV$temp=='inc+fluc'] <- 'Increase + Fluctuation'

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
ggsave(plot = last_plot(), file = here('output/ExtendedData_FigureS1_NBESmetrics.tiff'), width = 8, height = 6)


#