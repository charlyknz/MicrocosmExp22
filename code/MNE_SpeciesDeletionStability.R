#### R script start ####
# 09.12.2022
# by Charlotte Kunze

# packages
library(tidyverse)
library(readxl)
library(MESS)
library(here)

#### import OD data ####
allData_wBlankExp22 <- read_csv("~/Desktop/Exp22/MicrocosmExp22/Data/allData_wBlankExp22.csv", 
                                locale = locale())


### Data Wrangling ###

exp <- allData_wBlankExp22 %>%
  select(-fileName, -...1) %>%
  mutate(corrected_value = value - blank) %>%
  separate(ID,into = c('algae', 'nut', 'rep'))


### OD only 
names(exp)
OD_exp <- exp %>%
  group_by(method, treat, sampling,algae, nut, combination ) %>%
  filter(method == '440') %>%
  summarise(mean.OD = mean(corrected_value, na.rm = T),
            sd.OD = sd(corrected_value, na.rm = T),
            se.OD = sd.OD/sqrt(n())) %>%
  ungroup() %>%
  rename(species_com = combination,
         combination = algae, 
         temp=treat)  %>%
  mutate(missing = paste(ifelse(combination == 'ADGR', 'T', ifelse(combination == 'AGRT', 'D', ifelse(combination == 'DGRT', 'A', ifelse(combination == 'ADGT', 'R', ifelse(combination == 'ADRT', 'G', ifelse(combination == 'MIX', '0','NA'))))))))
OD_exp$temp[OD_exp$temp == 'cs'] <- 'CS'
OD_exp$temp[OD_exp$temp == 'flucinc'] <- 'inc+fluc'
OD_exp$temp[OD_exp$temp == 'fluc'] <- 'fluct'
OD_exp$nut[OD_exp$nut == 'Hi'] <- 'HI'


#### Species Deletion Stability SDS ####

#1. Calculate RR (mean-con)/(mean+con)
#trt data
data1 <- OD_exp %>%
  filter(species_com %in% c('mix', 'quattro')) %>%
  filter(temp != 'CS') 

# control data
data1_con <- OD_exp %>%
  filter(species_com %in% c('mix', 'quattro')) %>%
  filter(temp == 'CS')%>%
  rename(con.OD = mean.OD) %>%
  select(species_com, nut, sampling, missing, con.OD)

#join data and calcualte RR
all_data1 <- left_join(data1, data1_con, by = c('species_com', 'nut', 'sampling', 'missing')) %>%
  mutate(RR = (mean.OD-con.OD)/(mean.OD+con.OD)) 

#remove unnecessary df
rm(data1_con)
rm(data1)


## for SDS we need RRquattro and RR mix

#mix Data
mix_data1 <- all_data1 %>%
  filter(species_com == 'mix') %>%
  rename(MixRR = RR) %>%
  select(sampling, nut, temp, MixRR)

#join again and calculate delta SDS and SDS ratio
all_mix_data1 <- all_data1 %>%
  filter(species_com!= 'mix') %>%
  left_join(., mix_data1, by = c('sampling', 'nut', 'temp')) %>%
  mutate(UniqueID = paste(nut, temp, missing, sep = '_'))


#rm(mix_data1)
unique(all_mix_data1$temp )

#### AUC SDS ####
stab.auc.SDS <- data.frame()

#SDS = RR/MixRR,
#deltaSDS = AUC.RR-AUC.MixRR
UniqueID <- unique(all_mix_data1$UniqueID)
 for(i in 1:length(UniqueID)){
  temp<-all_mix_data1[all_mix_data1$UniqueID==UniqueID[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$sampling, temp$RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                 type = c("spline"),absolutearea = FALSE)
    AUC.MixRR<-auc(temp$sampling, temp$MixRR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                type = c("spline"),absolutearea = FALSE)
    AUC.deltaSDS <- AUC.RR - AUC.MixRR
    stab.auc.SDS<-rbind(stab.auc.SDS,
                        data.frame(temp[1,c(1:14)],
                                   AUC.MixRR, AUC.RR,AUC.deltaSDS))
    rm(temp)
  }
}

str(stab.auc.SDS)

stab.auc.SDS 
unique(stab.auc.SDS$missing)

#### plot AUC ####
ggplot(subset(stab.auc.SDS, nut == 'Low'), aes(x = missing, y = AUC.deltaSDS))+
  geom_col()+
  facet_wrap(~temp, scales = 'free_y')+
  theme_bw()
#ggsave(plot = last_plot(), file = here('output/AUCSDS.png'), width = 8, height = 3)


#### count data ####
### import  data ###
expData<- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/Exp22_wrangledData.csv') %>%
  select(-X)

BioVolume_exp22 <- read_excel("~/Desktop/Exp22/MicrocosmExp22/Data/BioVolume_exp22.xlsx",sheet = "Sheet1") %>%
  rename(speciesID = species) %>%
  select(-sampling,-comment, -Ansatz)

names(BioVolume_exp22)
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Ditylum']<-'DityCux'
BioVolume_exp22$speciesID[BioVolume_exp22$speciesID == 'Thala']<-'ThalaCux'
names(expData)

expData[ expData == 'TRUE' ] <- 'T'
str(expData)

### Data wrangle ###
countData1<- expData%>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = ''), ifelse(species == 'duo', paste(A1, A2, sep = ''), 
                                                                            ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = ''), paste(A1, A2, A3, A4, A5, sep = '')))))%>%
  ungroup() %>%
  select(-A1, -A2,-A3, -A4, -A5, -Afield, -At_mm, -plate, -GF, -counts, -magn,-V_ml) %>%
  filter(!species %in% c('mono', 'duo') & nut == 'Low') %>%
  merge(., BioVolume_exp22, by = c('speciesID')) %>%
  mutate(cellVolume = cells_ml * BV) %>%
  select(-cells_ml) %>%
  rename(species_com = species)  %>%
  mutate(missing = paste(ifelse(combination == 'ADGR', 'T', ifelse(combination == 'AGRT', 'D', ifelse(combination == 'DGRT', 'A', ifelse(combination == 'ADGT', 'R', ifelse(combination == 'ADRT', 'G', ifelse(combination == 'MIX', '0','NA'))))))))


#### Species Deletion Stability counts ####

#1. Calculate RR (mean-con)/(mean+con)
#trt data
dat1 <- countData1 %>%
  filter(species_com %in% c('MIX', 'quattro')) %>%
  filter(temp != 'CS') %>%
  select(sampling, nut,temp, species_com, missing,rep,cellVolume) %>%
  rename(trt.M = cellVolume)

# control data
names(countData1)
dat1_con <- countData1 %>%
  filter(species_com %in% c('MIX', 'quattro')) %>%
  filter(temp == 'CS')%>%
  group_by(sampling, nut, species_com, missing) %>%
  summarise(con.M = mean(cellVolume, na.rm = T))

#join data and calcualte RR
all_dat1 <- left_join(dat1, dat1_con, by = c('species_com', 'nut', 'sampling', 'missing')) %>%
  mutate(RR = (trt.M-con.M)/(trt.M+con.M)) 

#remove unnecessary df
rm(dat1_con)
rm(dat1)


## for SDS we need RRquattro and RR mix

#mix Data
mix_dat1 <- all_dat1 %>%
  filter(species_com == 'MIX') %>%
  rename(MixRR = RR) %>%
  select(sampling, nut, temp, MixRR,rep)

#join again and calculate delta SDS and SDS ratio
all_mix_dat1 <- all_dat1 %>%
  filter(species_com!= 'MIX') %>%
  left_join(., mix_dat1, by = c('sampling', 'nut', 'temp','rep')) %>%
   mutate(UniqueID = paste(nut, temp, missing,rep,sep = '_'))



#rm(mix_data1)
unique(all_mix_dat1$temp )
names(all_mix_dat1)
#### AUC SDS ####
stab.SDS <- data.frame()

UniqueID <- unique(all_mix_dat1$UniqueID)
for(i in 1:length(UniqueID)){
  temp<-all_mix_dat1[all_mix_dat1$UniqueID==UniqueID[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$sampling, temp$RR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                type = c("spline"),absolutearea = FALSE)
    AUC.MixRR<-auc(temp$sampling, temp$MixRR,  from = min(temp$sampling, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                   type = c("spline"),absolutearea = FALSE)
    AUC.deltaSDS <- AUC.RR-AUC.MixRR
    reverseSDS <- AUC.MixRR - AUC.RR
    stab.SDS<-rbind(stab.SDS,
                        data.frame(temp[1,c(2:11)],
                                   AUC.RR, AUC.MixRR,AUC.deltaSDS,reverseSDS))
    rm(temp)
  }
}

str(stab.SDS)

unique(stab.SDS$missing)

#### plot AUC ####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")
stab.SDS$temp[stab.SDS$temp=='fluct'] <- 'Fluctuation'
stab.SDS$temp[stab.SDS$temp=='inc'] <- 'Increase'
stab.SDS$temp[stab.SDS$temp=='inc+fluc'] <- 'IncreaseFluctuation'

stab.SDS %>%
  group_by(missing, temp) %>%
  summarise(mean = mean(AUC.deltaSDS,na.rm = T),
           sd =  sd(AUC.deltaSDS,na.rm = T),
           se = sd/sqrt(n())) %>%
ggplot(., aes(x = missing, y = mean, fill = missing))+
  geom_hline(yintercept = 0, color = 'darkgrey')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = missing), width = .3)+
  geom_col()+
  facet_wrap(~temp, scales = 'free_y')+
  labs(x = 'Species', y = 'Species Gross Contribution')+
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
  guides(color = guide_legend(override.aes = list(size = 3.5)))

#ggsave(plot = last_plot(), file = here('output/countAUCSDS.png'), width = 8, height = 3)

stab.SDS %>%
  group_by(missing, temp) %>%
  summarise(mean = mean(AUC.SDS,na.rm = T),
            sd =  sd(AUC.SDS,na.rm = T),
            se = sd/sqrt(n())) %>%
ggplot(., aes(x = missing, y = mean, fill = missing))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = missing), width = .3)+
  geom_col()+
  facet_wrap(~temp, scales = 'free_y')+
  labs(x = 'Species', y = 'Mean SDS +-SE', title = 'Log(AUC.RR exclusion/ AUC.RR mix)')+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()
#ggsave(plot = last_plot(), file = here('output/countAUClogSDS.png'), width = 8, height = 3)


#### SGC and difference in observed Instab of 2 species and monocultures ####
species1_2 <- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/allSpecies_MeanPresence2_1.csv')
str(species1_2)

SGC_data <- stab.SDS %>%
  rename(present = missing) %>%
  filter(nut == 'Low') %>%
  group_by(present, temp) %>%
  summarise(mean.SDS = mean(reverseSDS,na.rm = T),
            sd =  sd(reverseSDS,na.rm = T),
            se.SDS = sd/sqrt(n())) %>%
  left_join(., species1_2, by = c('present', 'temp'))

names(SGC_data)
hist(SGC_data$mean.SDS)

####plot SGC species12 ####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")

SGCplot <-ggplot(SGC_data,aes(x = mean.presence, y = mean.SDS, color = present, shape = temp))+
  geom_point(size = 2)+
  geom_errorbarh(aes(xmin = mean.presence-se.presence, xmax = mean.presence + se.presence), width = .5)+
  geom_errorbar(aes(ymin = mean.SDS-se.SDS, ymax = mean.SDS + se.SDS), height = .5)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Mean Species effect (2-1)', y = 'Mean Species Gross Cont (5-4)', color = 'SpeciesID', shape = 'Treatment')+
  scale_x_continuous(limits=c(-9,9), breaks = seq(-4,4,2))+
  scale_color_manual(values=cbbPalette)+
  #facet_wrap(~temp)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 12,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 12, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right')
SGCplot
ggsave(plot= SGCplot, file = here('MicrocosmExp22/output/SGC_Instab_fig.png'), width = 8, height = 6)



#### SGC and NBE of species 12 together ####

species12 <- read.csv2('species12_SGC.csv') %>%
  select(-X)
species12$speciesID[species12$speciesID == 'Asterio'] <- 'A'
species12$speciesID[species12$speciesID == 'DityCux'] <- 'D'
species12$speciesID[species12$speciesID == 'Guido'] <- 'G'
species12$speciesID[species12$speciesID == 'Rhizo'] <- 'R'
species12$speciesID[species12$speciesID == 'ThalaCux'] <- 'T'


unique(species12$speciesID)
str(stab.SDS)
str(species12)

SGC_data <- stab.SDS %>%
  rename(speciesID = missing) %>%
  filter(nut == 'Low') %>%
  group_by(speciesID, temp) %>%
  summarise(mean.SDS = mean(AUC.SDS,na.rm = T),
            sd =  sd(AUC.SDS,na.rm = T),
            se.SDS = sd/sqrt(n())) %>%
  left_join(., species12, by = c('speciesID', 'temp'))

names(SGC_data)

####plot SGC species12 ####
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")

SGCplot <-ggplot(SGC_data,aes(x = mean.NBE, y = mean.SDS, color = speciesID, shape = temp))+
  geom_point(size = 2)+
  geom_errorbarh(aes(xmin = mean.NBE-se.NBE, xmax = mean.NBE + se.NBE), width = .5)+
  geom_errorbar(aes(ymin = mean.SDS-se.SDS, ymax = mean.SDS + se.SDS), height = .5)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'mean NBE (SR = 2)', y = 'Mean Species Gross Cont (4-5)', color = 'SpeciesID', shape = 'Treatment')+
  scale_x_continuous(limits=c(-6,6), breaks = seq(-4,4,2))+
  scale_color_manual(values=cbbPalette)+
  #facet_wrap(~temp)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 12,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 12, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'right')
SGCplot
ggsave(plot= SGCplot, file = here('MicrocosmExp22/output/SGC_figure.png'), width = 8, height = 6)
