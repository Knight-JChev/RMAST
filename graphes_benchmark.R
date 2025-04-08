library(ggplot2)

# Idée de plot 1 ####
metricsframe = data.frame()
for (file in list.files("Benchmark_data")){
  currentMetrics = readRDS(paste0("Benchmark_data/",file))
  
  df = data.frame(measures = c(currentMetrics$allAccuracy,
                                         currentMetrics$allPrecision,
                                         currentMetrics$allRecall),
                            type = as.factor(
                              c(rep("acc", currentMetrics$nbRepeats),
                                rep("pre", currentMetrics$nbRepeats),
                                rep("rec", currentMetrics$nbRepeats))),
                            nbErrors = as.factor(
                              rep(currentMetrics$nb.move, currentMetrics$nbRepeats*3)))
  metricsframe = rbind(metricsframe, df)
}  

supp.labs <- c("Accuracy", "Precision", "Recall")
names(supp.labs) <- c("acc","pre","rec")
  
plot = ggplot(data = metricsframe, aes(x = nbErrors, y = measures, 
                                       fill = type, colour = type, alpha = nbErrors)) +
  geom_violin(position = "dodge") +
  facet_grid(.~type, labeller = as_labeller(supp.labs)) +
  labs(title = "Violin plot des métriques" ,
       x = "Nombre d'erreurs",
       y = "Valeur de métrique") +
  theme_bw() + theme(panel.grid.major = element_line(colour = "grey80"),
                     axis.ticks = element_blank(),
                     panel.grid.minor.x=element_blank(),
                     panel.grid.major.x=element_blank()) +
  scale_fill_manual(values = c("mediumpurple2", "aquamarine2", "sienna2")) +
  scale_colour_manual(values = c("darkorchid1", "darkslategray3", "orangered")) +
  scale_alpha_manual(values = c(1, 0.8,0.6,0.4,0.2)) +
  coord_cartesian(ylim = c(0,1)) +
  guides(alpha = 'none', fill = "none", colour = "none")
plot


## Idée de plot 2 ####
metricsframe = data.frame()
for (file in list.files("Benchmark_data")){
  currentMetrics = readRDS(paste0("Benchmark_data/",file))
  
  df = data.frame(Accuracy = currentMetrics$allAccuracy,
                  Precision = currentMetrics$allPrecision,
                  Recall = currentMetrics$allRecall,
                  nbErrors = as.factor(
                    rep(currentMetrics$nb.move, currentMetrics$nbRepeats)))
  metricsframe = rbind(metricsframe, df)
} 

for (colname in colnames(metricsframe)[1:3]){
  plot = ggplot(data = metricsframe, aes(x = nbErrors, y = colname)) +
    geom_smooth(aes(fill = "red")) +
    labs(title = "Violin plot des métriques" ,
         x = "Nombre d'erreurs",
         y = "Valeur de métrique") +
    theme_bw() + theme(panel.grid.major = element_line(colour = "grey80"),
                       axis.ticks = element_blank(),
                       panel.grid.minor.x=element_blank(),
                       panel.grid.major.x=element_blank()) +
    coord_cartesian(ylim = c(0,1)) +
  print(plot)
}


# Tests
d <- ggplot(mpg, aes(fl))

df <- data.frame(grp = c("A", "B"), fit = 4:5, se = 1:2)
j <- ggplot(df, aes(grp, fit, ymin = fit - se, ymax = fit + se))

# patchwork package pour les insets 
n <- d + geom_bar(aes(fill = fl))
n + scale_fill_manual(
  values = c("skyblue", "royalblue", "blue", "navy"), 
  limits = c("d", "e", "p", "r"), breaks =c("d", "e", "p", "r"), 
  name = "fuel", labels = c("D", "E", "P", "R"))
n
