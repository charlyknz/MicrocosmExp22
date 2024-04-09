#### R script to plot temperature curves ####
### required packages ###
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)
library(here)

#### Import data ####
df <- read_csv("Data/all_temperatures.csv", locale = locale())
str(df)


## data wrangling
temp <- df %>%
  mutate(datetime = timestamp) %>%
  separate(timestamp, into = c('date', 'time'), sep = ' ') %>%
  filter(date > '2022-01-19') %>%
  filter(date < '2022-02-18') %>%
  mutate(time1 = as.character(time)) 


# set color palette
tempPalette <- c('black' ,"#E41A1C","#377EB8","#4DAF4A" ) # temp treatments

# order mesocosms after treatments
dataPlot <- temp%>%
  filter( unit %in%c(8,11,2,5)) %>%
  mutate(treat = ifelse(unit == 11, 'Constant', ifelse(unit == 2, 'Increase', ifelse(unit == 8, 'Fluctuation', 'Increase + Fluctuation'))))
levels(as.factor(dataPlot$treat))

#temperature plots
plot <- ggplot(dataPlot, aes(x = datetime, y = actual_tempmiddle, color = treat))+
  geom_line(linewidth = 1.0)+
  facet_wrap(~treat,  ncol = 1)+
  labs( x = 'Date', y = 'Temperature (in °C)')+
  scale_x_datetime(breaks = '7 days')+
  scale_color_manual(values = tempPalette)+
  theme_classic()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.text = element_text(face = 'bold', size = 18),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank())+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 11,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 11,  colour = "black", angle = 0, hjust = 0.4))
plot
ggsave(plot = plot, file = here('output/FigS1_Temp_curves_exp22.tiff'), height = 12, width = 5)
