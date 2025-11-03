source(here("~/Desktop/phD/Exp22/MicrocosmExp22/code/03_NBES_calculation.R"))

str(all_resistance)
str(CV)

remove1<-stab.auc.mix%>%
  distinct(combination, rep,temp, AUC.RR_obs)

AllStabil<-stab.auc %>% 
  distinct(combination, rep,temp, AUC.RR_obs) %>% 
  bind_rows(.,remove1) %>% 
  left_join(., all_resistance) %>% 
  left_join(., CV)

plot.a <- ggpubr::ggscatter(AllStabil, x = "CVobs", y = "AUC.RR_obs", 
                            xlab = "CV", ylab = "OEV", cor.method = "spearman", add = "reg.line", cor.coef = T)
plot.b <- ggpubr::ggscatter(AllStabil, x = "delta_ges", y = "AUC.RR_obs", 
                            ylab = "OEV", xlab = "Resistance", cor.method = "spearman", add = "reg.line", cor.coef = T)

cowplot::plot_grid(plot.a, plot.b, labels = c("a)", "b)"))
ggsave(plot = last_plot(), file = here("output/StabilityCorrelation.png"), width = 8, height = 4)
