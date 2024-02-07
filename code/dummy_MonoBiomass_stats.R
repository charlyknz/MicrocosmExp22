rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by(species, combination, temp) %>%
  filter(species == 'mono')%>%
  mutate(transV = cellVolume/10^9)%>%
  summarise(meanV = sum(transV, na.rm = T),
            sdV = sd(transV, na.rm = T),
            seV = sdV/sqrt(n())) %>%
  ggplot(., aes( x = temp, y = meanV, fill = temp, group = temp))+
  geom_col(alpha = 0.8)+
  geom_errorbar(aes(ymin = meanV-seV, ymax = meanV+seV), width = .8)+
 # labs(y =  expression(Mean~Biovolume~x~10^9~'['~µm^3~'/'~ml~']'), x = 'Time [days]', color = 'Temperature')+
  facet_grid(~combination, scales = 'free_y')+
  scale_fill_manual(values = tempPalette)+
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


monoData<-rawData%>%
  mutate(day = sampling * 2,
         day = as.numeric(paste(ifelse(sampling == 1, 1, sampling*2))))%>%
  group_by(species, combination, temp,rep) %>%
  filter(species == 'mono')%>%
  summarise(sum = sum(cellVolume, na.rm = T))

A <- filter(monoData, combination == 'A')
D <- filter(monoData, combination == 'D')
G <- filter(monoData, combination == 'G')
R <- filter(monoData, combination == 'R')
'T' <- filter(monoData, combination == 'T')

aov1<- aov(sum ~ temp, A)
aov2<- aov(sum ~ temp, D)
aov3<- aov(sum ~ temp, G)
aov4<- aov(sum ~ temp, R)
aov5<- aov(sum ~ temp, T)

summary(aov1)
TukeyHSD(aov1)

summary(aov2)
TukeyHSD(aov2)

summary(aov3)
TukeyHSD(aov3)

summary(aov4)
TukeyHSD(aov4)

summary(aov5)
TukeyHSD(aov5)
 