## R Script to gather temperature data 

#### required packages ####
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)

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
label_temp <- c('11' = 'constant','12' = 'constant', '2' = 'increase','6' = 'increase', '5' = 'fluctuation+increase', '10' = 'fluctuation+increase',
                '7' = 'fluctuation','8' = 'fluctuation')
plot <- ggplot(subset(temp, unit %in%c(11,2,5,8)), aes(x = datetime, y = actual_tempmiddle))+
  geom_line(size = 1.0)+
  facet_wrap(~unit, labeller = labeller(unit = label_temp), ncol = 2)+
  labs( x = 'Date', y = 'Temperature (in °C)')+
  scale_x_datetime(breaks = '5 days')+
  theme_classic()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank())
plot
ggsave(plot = plot, file = 'FigSTemp_curves_exp22.png', width = 10, height = 8)
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##



