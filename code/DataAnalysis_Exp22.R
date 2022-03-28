#### R script to analyse counting data for Microcosm Exp 22####
#by Charlotte Kunze 28.03.22

## load packages
library(tidyverse)
library(readxl)
library(scales)
library(here)
#### import data ####
counting_ID_sampling15 <- read_excel("MicrocosmExp22/Data/counting_ID_sampling15.xlsx", 
                                     sheet = "sampling15") %>%
  select(-cells_mL)

str(counting_ID_sampling15)

#### import meta data ####
ID_herma <- read_excel("MicrocosmExp22/Data/counting_ID_sampling15.xlsx",  sheet = "overview") %>%
  select(-plate, -V_ml, - magn, -sampling) 
ID_herma[ ID_herma == 'Asterio' ] <- 'A'
ID_herma[ ID_herma == 'DityCux' ] <- 'D'
ID_herma[ ID_herma == 'Guido' ] <- 'G'
ID_herma[ ID_herma == 'ThalaCux' ] <- 'T'
ID_herma[ ID_herma == 'Rhizo' ] <- 'R'


names(ID_herma)


counting_ID_sampling15 <- left_join(counting_ID_sampling15, ID_herma, by = c('no', 'temp', 'nut', 'rep', 'species'))
#### calculate cells_ml ####
cells_mL_calculation_Utermo_hl_05chamber <- read_excel("MicrocosmExp22/cells_mL_calculation_Utermöhl_05chamber.xlsx")
str(cells_mL_calculation_Utermo_hl_05chamber)
data <- left_join(counting_ID_sampling15,cells_mL_calculation_Utermo_hl_05chamber, by = c('magn', 'V_ml'))

data <- data %>%
  mutate(cells_ml  = counts*At_mm/Afield*GF*V_ml)
names(data)

data$col  <- NA


mono <- filter(data, species == 'mono')%>%
  mutate(combination = A1)
duo <- filter(data, species == 'duo') %>%
  mutate(combination = paste(A1, A2, sep = ''))
quattro <- filter(data, species == 'quattro')%>%
  mutate(combination = paste(A1, A2, A3, A4, sep = ''))
mix <- filter(data, species == 'MIX')%>%
  mutate(combination = paste(A1, A2, A3, A4, A5, sep = ''))



#### first plots ####
# The palette with black:
cbbPalette <- c("#E69F00", "#000000","#0072B2", "#009E73","#CC79A7")

duo %>%
  group_by(combination, temp, nut, speciesID, col) %>%
  summarise(mean = mean(cells_ml, na.rm = T)) %>%
ggplot(., aes(x = combination, y = mean))+
  geom_col(aes(fill = speciesID))+
  facet_grid(~nut~temp, scales = 'free_x')+
  scale_fill_manual(values=cbbPalette)+
  scale_y_log10("cells [mL]",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs( title = '2 species combination', y = 'cells [ml]', x = 'Nutrient level')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/duo_abosoluteCellsML.png'), width = 10, height = 4 )


quattro %>%
  group_by(combination, temp, nut, speciesID, col) %>%
  summarise(mean = mean(cells_ml, na.rm = T)) %>%
  ggplot(., aes(x = combination, y = mean))+
  geom_col(aes(fill = speciesID))+
  facet_grid(~nut~temp, scales = 'free_x')+
  scale_fill_manual(values=cbbPalette)+
  scale_y_log10("cells [mL]",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs( title = '4 species combination', y = 'cells [ml]', x = 'Nutrient level')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/quattro_abosoluteCellsML.png'), width = 10, height = 4 )

mix %>%
  group_by(combination, temp, nut, speciesID, col) %>%
  summarise(mean = mean(cells_ml, na.rm = T)) %>%
  ggplot(., aes(x = nut, y = mean))+
  geom_col(aes(fill = speciesID))+
  facet_grid(~temp, scales = 'free_x')+
  scale_fill_manual(values=cbbPalette)+
  scale_y_log10("cells [mL]",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs( title = '5 species combination', y = 'cells [ml]', x = 'Nutrient level')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/mix_abosoluteCellsML.png'), width = 10, height = 4 )


#### relative abundance plots ####
duo %>%
  group_by(combination, temp, nut) %>%
  mutate(sum = sum(cells_ml)) %>%
  ungroup()%>%
  mutate(rel_ab = cells_ml/ sum) %>%
  group_by(combination, temp, nut, speciesID, rep) %>%
  summarise(mean = mean(rel_ab, na.rm = T)) %>%
  ggplot(., aes(x = combination, y = mean))+
  geom_col(aes(fill = speciesID))+
  facet_grid(~nut~temp, scales = 'free_x')+
  scale_fill_manual(values=cbbPalette)+
   labs( title = '2 species combination', y = 'relative abundance [cells mL-1]', x = 'Nutrient level')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/duo_relCellsML.png'), width = 10, height = 4 )


quattro %>%
  group_by(combination, temp, nut) %>%
  mutate(sum = sum(cells_ml)) %>%
  ungroup()%>%
  mutate(rel_ab = cells_ml/ sum) %>%
  group_by(combination, temp, nut, speciesID, rep) %>%
  summarise(mean = mean(rel_ab, na.rm = T)) %>%
  ggplot(., aes(x = combination, y = mean))+
  geom_col(aes(fill = speciesID))+
  facet_grid(~nut~temp, scales = 'free_x')+
  scale_fill_manual(values=cbbPalette)+
  labs( title = '4 species combination', y = 'relative abundance [cells mL-1]', x = 'Nutrient level')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/quattro_relativeCellsML.png'), width = 10, height = 4 )

mix %>%
  group_by(combination, temp, nut) %>%
  mutate(sum = sum(cells_ml)) %>%
  ungroup()%>%
  mutate(rel_ab = cells_ml/ sum) %>%
  group_by(combination, temp, nut, speciesID, rep) %>%
  summarise(mean = mean(rel_ab, na.rm = T)) %>%
  ggplot(., aes(x = nut, y = mean))+
  geom_col(aes(fill = speciesID))+
  facet_grid(~temp, scales = 'free_x')+
  scale_fill_manual(values=cbbPalette)+
   labs( title = '5 species combination', y = 'rel. abundance [cells mL-1]', x = 'Nutrient level')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MicrocosmExp22/output/mix_relativeCellsML.png'), width = 10, height = 4 )
