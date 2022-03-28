### R script to analyse exp 22 data ####
# by Charlotte Kunze

#load packages
library(tidyverse)

#----------------------------------------------------------------------------------------------------------------#

# import data
allData_wBlankExp22 <- read_csv("~/Desktop/Exp22/output/allData_wBlankExp22.csv", 
                                   locale = locale())

#----------------------------------------------------------------------------------------------------------------#

#### Data Wrangling ####

exp <- allData_wBlankExp22 %>%
  select(-fileName, -X1) %>%
  mutate(corrected_value = value - blank) %>%
  separate(ID,into = c('algae', 'nut', 'rep'))


### OD only 
OD_exp <- exp %>%
  group_by(method, treat, combination, sampling,algae, nut ) %>%
  filter(method == '440') %>%
  summarise(mean = mean(corrected_value, na.rm = T),
            sd = sd(corrected_value, na.rm = T),
            se = sd/sqrt(n()))



ggplot(filter(OD_exp, combination == 'mono'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/mono_OD.png', width = 9, height = 6)

ggplot(filter(OD_exp, combination == 'duo1'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/duo1_OD.png', width = 9, height = 6)

ggplot(filter(OD_exp, combination == 'duo2'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/duo2_OD.png', width = 9, height = 6)

ggplot(filter(OD_exp, combination == 'quattro'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/quattro_OD.png', width = 9, height = 6)

ggplot(filter(OD_exp, combination == 'mix' & !sampling %in%c(5,7)), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/mix_OD.png', width = 9, height = 6)



### RFU only ###
rfu_exp <- exp %>%
  group_by(method, treat, combination, sampling,algae, nut ) %>%
  filter(method == '395680') %>%
  summarise(mean = mean(corrected_value, na.rm = T),
            sd = sd(corrected_value, na.rm = T),
            se = sd/sqrt(n()))

### rfu plots ###

ggplot(filter(rfu_exp, combination == 'mono'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/mono_rfu.png', width = 9, height = 6)

ggplot(filter(rfu_exp, combination == 'duo1'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/duo1_rfu.png', width = 9, height = 6)

ggplot(filter(rfu_exp, combination == 'duo2'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/duo2_rfu.png', width = 9, height = 6)

ggplot(filter(rfu_exp, combination == 'quattro'), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(algae~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/quattro_rfu.png', width = 9, height = 6)

ggplot(filter(rfu_exp, combination == 'mix' & !sampling %in%c(5,7)), aes(x = sampling, y = mean, color = nut))+
  geom_point()+
  geom_line(linetype = 'dashed')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + sd), width = .8)+
  facet_grid(~treat, scales = 'free_y')+
  theme_bw()
ggsave(plot = last_plot(), 'output/mix_rfu.png', width = 9, height = 6)

