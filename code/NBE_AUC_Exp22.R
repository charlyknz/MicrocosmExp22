#### R script to calculate AUC  ####
# by Charlotte Kunze 

#load packages and functions
library(reshape2)
library(tidyverse)
library(MESS) #https://rdrr.io/cran/MESS/man/auc.html
library(nlme)
library(ggpubr)
library(psych)
library(cowplot)


#### import data ####
expData<- read.csv('~/Desktop/Exp22/MicrocosmExp22/Data/Exp22_wrangledData.csv') %>%
  select(-X) %>%
  filter(nut == 'Low')
expData[ expData == 'TRUE' ] <- 'T'
str(expData)

#### data ####

#calculate mean cells_ml in control dataset 
expData_control <- expData%>%
  filter(temp == 'CS') %>%
  rename(con.cells = cells_ml) %>%
  mutate(combination = ifelse(species == 'mono',paste(A1, sep = '_'), ifelse(species == 'duo', paste(A1, A2, sep = '_'), 
                      ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = '_'), paste(A1, A2, A3, A4, A5, sep = '_')))))%>%
  group_by(nut, species, speciesID, sampling, combination) %>%
  summarise(con.bio = mean(con.cells, na.rm = T))

expData<- expData%>%
    filter(temp != 'CS')%>%
    mutate(combination = ifelse(species == 'mono',paste(A1, sep = '_'), ifelse(species == 'duo', paste(A1, A2, sep = '_'), 
                                                                               ifelse(species == 'quattro', paste(A1, A2, A3, A4, sep = '_'), paste(A1, A2, A3, A4, A5, sep = '_')))))
  names(expData)

#bringing in the control data to the main data
exp.data2 <- dplyr::right_join(expData, expData_control, 
                               by = c('nut','species', 'speciesID', 'sampling', 'combination'))
  

# for comparison, calculate the overall response                                    ## data.tot -> Total Response on Forcing
exp.data.tot<-exp.data2%>% 
  dplyr::group_by(nut,sampling, temp, rep,combination) %>%
  dplyr:: summarise(treat.tot = sum(cells_ml, na.rm = T),
                    con.tot = sum(con.bio, na.rm = T))

# LRR (Log Response Ratio) of community 
exp.data.tot$LRR.tot<-log(exp.data.tot$treat.tot/exp.data.tot$con.tot)
exp.data.tot$deltabm.tot<-(exp.data.tot$treat.tot-exp.data.tot$con.tot)/(exp.data.tot$treat.tot+exp.data.tot$con.tot)
hist(exp.data.tot$LRR.tot)
hist(exp.data.tot$deltabm.tot)

#create data set for species specific effect sizes
#take out those rows where cells_ml is 0 in both treatment (cells_ml) + control (con.bio)
names(exp.data.tot)
names(exp.data2)
exp.data3 <- exp.data2 %>%
  dplyr::left_join(.,exp.data.tot, by = c("temp","nut","combination","sampling",  'rep'))
exp.data3$RR <- exp.data3$cells_ml+exp.data3$con.bio
exp.data3 <- dplyr::filter(exp.data3, RR!=0)

# create species specific Log Response Ratio (LRR) for cells_ml
exp.data3$LRR<-log(exp.data3$cells_ml/exp.data3$con.bio)


# create species specific change in cells_ml  (RR) 
exp.data3$RR<-(exp.data3$cells_ml-exp.data3$con.bio)/(exp.data3$cells_ml+exp.data3$con.bio)

# the absence of species in control or treatment creates INF, get rid of these
exp.data3$LRR[exp.data3$LRR=="Inf"]<-NA
exp.data3$LRR[exp.data3$LRR=="-Inf"]<-NA

# create species specific contribution (pi) to cells_ml for control and treat  (treat / tot)         ## contribution 
exp.data3$treat.pi <- exp.data3$cells_ml/exp.data3$treat.tot
exp.data3$con.pi <- exp.data3$con.bio/exp.data3$con.tot

# create species specific LRR and difference for pi           (delta pi)            
exp.data3$delta.pi <- exp.data3$treat.pi-exp.data3$con.pi


# weight LRR  -  mean pi                                                            
exp.data3$mean.pi <- 0.5*(exp.data3$treat.pi+exp.data3$con.pi)
exp.data3$LRR.w <- exp.data3$LRR*exp.data3$mean.pi
exp.data3$RR.w <- exp.data3$RR*exp.data3$mean.pi

# create the deviance between species and community effect sizes                     
exp.data3$LRR.diff <- exp.data3$LRR-exp.data3$LRR.tot                                
exp.data3$LRR.diff.w <- exp.data3$LRR.diff*exp.data3$mean.pi


# differentiate the most dominant species only

exp.mean.dom <- exp.data3 %>%
  dplyr:: group_by(speciesID) %>%
  dplyr:: summarise(mean=mean(mean.pi, na.rm = T))
exp.mean.dom$dom <- exp.mean.dom$mean  #*exp.mean.dom$length                       
exp.data3 <- dplyr::left_join(exp.data3, dplyr::select(exp.mean.dom, -mean), by = c("speciesID"))
hist(exp.data3$dom)
str(exp.data3)

#### LRR ####

exp.data3$LRR
names(exp.data3)
LRR.all <- exp.data3 %>%
  group_by(species,speciesID, sampling, combination, temp) %>%
  summarise(mean = mean(cells_ml, na.rm = T),
            mean.con = mean(con.bio, na.rm = T),
            sd = sd(cells_ml, na.rm = T),
            sd.con = sd(con.bio, na.rm = T),
            N = as.numeric(paste(3)),
            N.con = as.numeric(paste(3)),
            LRR = log(mean/mean.con),
            var.LRR = sd.con^2/(N.con*mean.con^2)+sd^2/(N*mean^2)) 

mono.lrr <- ggplot(filter(LRR.all, species == 'mono'), aes(x = sampling, y = LRR, color = speciesID, group = speciesID))+
  geom_point()+
  geom_errorbar(aes(ymin = LRR-var.LRR, ymax = LRR+var.LRR), width = .8)+
  geom_line(linetype = 'dashed', size = 0.5)+
  geom_hline(yintercept = 0)+
  facet_grid(~temp~combination)+
  theme_bw()
mono.lrr
#ggsave(plot = mono.lrr, file = 'MonoCulture_lrr.png')

duo.lrr <- ggplot(filter(LRR.all, species == 'duo'), aes(x = sampling, y = LRR, color = speciesID, group = speciesID))+
  geom_point()+
  geom_errorbar(aes(ymin = LRR-var.LRR, ymax = LRR+var.LRR), width = .8)+
  geom_line(linetype = 'dashed', size = 0.5)+
  geom_hline(yintercept = 0)+
  facet_grid(~temp~combination)+
  theme_bw()
duo.lrr
#ggsave(plot = duo.lrr, file = 'DuoCulture_lrr.png', width = 8, height = 5)

quattro.lrr <- ggplot(filter(LRR.all, species == 'quattro'), aes(x = sampling, y = LRR, color = speciesID, group = speciesID))+
  geom_point()+
  geom_errorbar(aes(ymin = LRR-var.LRR, ymax = LRR+var.LRR), width = .8)+
  geom_line(linetype = 'dashed', size = 0.5)+
  geom_hline(yintercept = 0)+
  facet_grid(~temp~combination)+
  theme_bw()
quattro.lrr
#ggsave(plot = quattro.lrr, file = 'QuattroCulture_lrr.png', width = 8, height = 5)


mix.lrr <- ggplot(filter(LRR.all, species == 'MIX'), aes(x = sampling, y = LRR, color = speciesID, group = speciesID))+
  geom_point()+
  geom_errorbar(aes(ymin = LRR-var.LRR, ymax = LRR+var.LRR), width = .8)+
  geom_line(linetype = 'dashed', size = 0.5)+
  geom_hline(yintercept = 0)+
  facet_grid(~temp~combination)+
  theme_bw()
mix.lrr
#ggsave(plot = mix.lrr, file = 'MixCulture_lrr.png', width = 3, height = 5)

#######=======================================================================#######




#### AUC Loop ####

names(exp.data3)
exp.data3 <- exp.data3 %>% 
  dplyr::mutate(USI = paste( no,temp, nut,species,combination,rep,speciesID,sep = "_"))     #create unique ID

#order data by day to start with 0 
exp.data3 <- exp.data3[order(exp.data3$sampling),]
names(exp.data3)
#create an empty data frame (then our loop is faster)
exp.stab.auc <- data.frame()


# the following loops cycle through all unique cases (USI)

# Zooplankton
USI <- unique(exp.data3$USI)
for(i in 1:length(USI)){
  temp<-exp.data3[exp.data3$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$sampling, temp$RR,  from = min(1, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                type = c("spline"),absolutearea = FALSE)
    AUC.pi<-auc(temp$sampling, temp$delta.pi, from = min(1, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                type = c("spline"),absolutearea = FALSE)
    AUC.totRR<-auc(temp$sampling, temp$deltabm.tot, from = min(1, na.rm = TRUE), to = max(temp$sampling, na.rm = TRUE),
                   type = c("spline"),absolutearea = FALSE)
    mean.delta.pi<-mean(temp$delta.pi)
    mean.RR<-mean(temp$RR)
    mean.pi<-mean(temp$mean.pi)
    sd.delta.pi<-sd(temp$delta.pi)
    sd.RR<-sd(temp$RR)
    exp.stab.auc<-rbind(exp.stab.auc,
                        data.frame(temp[1,c(1,2,3,4,5,10,13,22,31,37,39)],
                                   AUC.RR,
                                   AUC.pi,
                                   AUC.totRR,
                                   mean.delta.pi,
                                   mean.RR,
                                   mean.pi,
                                   sd.delta.pi,
                                   sd.RR))
    rm(temp)
  }
}

summary(exp.stab.auc)


#### Total RR ####
str(exp.stab.auc)
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")

exp.stab.auc$temp[exp.stab.auc$temp=='fluct'] <- 'Fluctuation'
exp.stab.auc$temp[exp.stab.auc$temp=='inc'] <- 'Increase'
exp.stab.auc$temp[exp.stab.auc$temp=='inc+fluc'] <- 'IncreaseFluctuation'

exp.stab.auc%>%
  distinct(species, temp, combi, no, rep, speciesID, AUC.totRR) %>%
  mutate(N = paste(ifelse(species == 'duo', 2, ifelse(species == 'quattro', 4,ifelse(species == 'mono', 1, 5)))))%>%
  group_by( temp,N)%>%
  mutate(mean.RR= mean(AUC.totRR, na.rm = T),
            sd.RR = sd(AUC.totRR, na.rm=T),
            se.RR = sd.RR/sqrt(n()))  %>%
  ggplot(. )+
  geom_hline(yintercept = 0)+
  geom_point(aes( x = N, y = AUC.totRR, color = temp), alpha = 0.4)+
  geom_errorbar(aes(x = N, y = AUC.totRR, ymin= mean.RR-se.RR, ymax = mean.RR +se.RR), color = 'black',width = .2)+
  geom_point(aes( x = N, y = mean.RR),color = 'black',size = 2.5)+
  scale_colour_brewer(palette = "Set1")+
  labs(y = 'Mean Instability', x = 'No of connections')+
  facet_grid(~temp)+
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
ggsave(plot = last_plot(), file = here('~/Desktop/Exp22/MicrocosmExp22/output/MeanInstability.png'), width = 8, height = 4)



### AUC 
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")

data.plot <- distinct(exp.stab.auc, species, temp, combination, no, rep, speciesID, AUC.RR, AUC.pi) %>%
  group_by(species, temp, combination, speciesID)%>%
  mutate(mean.rr= mean(AUC.RR, na.rm = T),
         sd.rr = sd(AUC.RR, na.rm =T),
         se.rr = sd.rr/sqrt(n()),
         mean.pi= mean(AUC.pi, na.rm = T),
         sd.pi = sd(AUC.pi, na.rm =T),         
         se.pi = sd.pi/sqrt(n()))

#adjust labels
data.plot$species[data.plot$species == 'duo']<- '2 species'
data.plot$species[data.plot$species == 'quattro']<- '4 species'
data.plot$species[data.plot$species == 'MIX']<- '5 species'

data.plot$temp[data.plot$temp=='fluct'] <- 'Fluctuation'
data.plot$temp[data.plot$temp=='inc'] <- 'Increase'
data.plot$temp[data.plot$temp=='inc+fluc'] <- 'IncreaseFluctuation'

levels(as.factor(data.plot$species))

ggplot(filter(data.plot, species != "mono"), aes(mean.pi, mean.rr,color = speciesID )) +
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymin = mean.rr - sd.rr, ymax = mean.rr + sd.rr), width = .2)+
  geom_errorbarh(aes(xmin = mean.pi - sd.pi, xmax = mean.pi + sd.pi), height = .2)+
  geom_vline(xintercept = 0, alpha = 0.25) +                                      
  geom_hline(yintercept = 0, alpha = 0.25) +
  facet_grid(~species~temp, scales = 'free') +
  scale_color_manual(values=cbbPalette)+
  labs(x ='Relative contribution to stability', y = "Absolute contribution to stability", shape = 'Exp') +  
  theme_bw()+
  theme(legend.position = 'right')+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 6)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave(plot = last_plot(), file = here('output/AUC_first.png'), width = 10, height = 8)

data.plot$speciesID[data.plot$speciesID=='ThalaCux'] <-'Thalassionema'
data.plot$speciesID[data.plot$speciesID=='DityCux'] <-'Ditylum'
data.plot$speciesID[data.plot$speciesID=='Guido'] <-'Guinardia'
data.plot$speciesID[data.plot$speciesID=='Rhizo'] <-'Rhizosolenia'
data.plot$speciesID[data.plot$speciesID=='Asterio'] <-'Asterionellopsis'



RR1 <- data.plot %>%
filter( species == '5 species') %>%
  distinct(species, speciesID,temp, mean.rr, mean.pi, se.rr, se.pi) %>%
ggplot(., aes( x = speciesID, y = mean.rr, fill = speciesID))+
  geom_hline(yintercept = 0, col = 'darkgrey')+
  geom_errorbar(aes(ymin = mean.rr - se.rr, ymax = mean.rr + se.rr, color = speciesID), width = .2)+
  geom_col()+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  facet_wrap(~temp)+
  labs(y = 'Absolute Contribution (AUC.RR+-SE)')+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
RR1

pi1 <- data.plot %>%
  filter( species == '5 species') %>%
  distinct(species, speciesID,temp, mean.rr, mean.pi, se.rr, se.pi) %>%
  ggplot(., aes( x = speciesID, y = mean.pi, fill = speciesID))+
  geom_hline(yintercept = 0, col = 'darkgrey')+
  geom_errorbar(aes(ymin = mean.pi - se.pi, ymax = mean.pi + se.pi, color = speciesID), width = .2)+
  geom_col()+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  facet_wrap(~temp)+
  labs(y = 'Relative Contribution (AUC.RR+-SE)', x = ' ')+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
pi1
plot_grid(pi1,RR1, ncol = 1)
ggsave(plot = last_plot(), file = here('output/AUC_speciesID_mix.png'), width =6, height =7)
