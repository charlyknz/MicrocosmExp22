## R Script to gather temperature data 

#### required packages ####
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)
library(here)

setwd("~/Desktop/Exp22/Temperature_Charly_Planktotrons_2022")

## import excel files using a Loop
directory <-  paste(getwd(),sep = '') # hier sind meine Daten

#read file (from [here][1])
file.list<- list.files(full.names = TRUE, pattern = "*.xls")

df.list <- lapply(file.list, read_excel)
df <- bind_rows(df.list, .id = "id")
df <- filter(df, !unit %in% c(3, 4,9))
str(df)
write.csv(x = df, file = 'all_temperatures.csv')
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
#### START ####
df <- read_csv("~/Desktop/Exp22/Temperature_Charly_Planktotrons_2022/all_temperatures.csv", 
                             locale = locale())
str(df)
## data wrangling

temp <- df %>%
  mutate(datetime = timestamp) %>%
  separate(timestamp, into = c('date', 'time'), sep = ' ') %>%
  filter(date > '2022-01-19') %>%
  filter(date < '2022-02-18') %>%
  mutate(time1 = as.character(time)) 

dev.off()

# order MC after treatments
tempPalette <- c('black' ,"#E41A1C","#4DAF4A","#377EB8"  ) # temp treatments

dataPlot <- temp%>%
  filter( unit %in%c(8,11,2,5)) %>%
  mutate(treat = ifelse(unit == 11, 'Constant', ifelse(unit == 2, 'Increase', ifelse(unit == 8, 'Fluctuation', 'Fluctuation + Increase'))))
levels(as.factor(dataPlot$treat))

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
ggsave(plot = plot, file = here('~/Desktop/Exp22/MicrocosmExp22/output/FigSTemp_curves_exp22.tiff'), height = 12, width = 5)
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##



