library(ggplot2)
library(dplyr)

# Charger l'exemple ####
metricsframe2 = readRDS("./Example_benchmark_data")

# OU Formatage des données du benchmark ####
dir = "Benchmark_data/200rep_3a7dist_1a25moves/"
metricsframe = data.frame()
for (file in list.files(dir)){
  if (!(dir.exists(paste0(dir,file)))){
    df = readRDS(paste0(dir,file))
    metricsframe = rbind(metricsframe, df)
  }
}
newmoves = c("1","3", "5", "10", "15", "20", "25")
metricsframe2 = arrange(mutate(metricsframe, nb.move = factor(nb.move, levels = newmoves)), nb.move)

## Matrice de couleur "normalisée" ####
mieux = metricsframe2 %>% group_by(nbLeaves, nbRepeats, nb.move, dist) %>%
  count(perfect, FN) %>%
  reframe(matval = FN*n) %>%
  group_by(nbLeaves, nbRepeats, nb.move, dist) %>%
  tally(matval) %>%
  mutate(normVal = round(n/(as.numeric(levels(nb.move)[nb.move])*nbRepeats), digits = 3)) %>%
  mutate(compte_faux = paste0("n = ",n))

ggplot(mieux) +
  geom_tile(aes(x = dist, y = nb.move, fill = normVal)) +
  geom_text(aes(x = dist, y = nb.move, label = normVal), col = "mistyrose", size = 6) +
  #geom_text(aes(x = dist, y = nb.move, label = compte_faux), col = "white", nudge_y = -0.3, size = 4) + 
  scale_fill_viridis_c(begin = 0.2, end = 0.6, direction = -1, option = "plasma")+
  labs (subtitle = paste0("200 arbres de 100 feuilles par case"),
        x = "Distance de déplacement",
        y = "Nombre de feuilles déplacées",
        fill = "Proportion \n d'erreur")+
  theme_bw() + theme(panel.background = element_blank(),
                     panel.border = element_blank(),
                     panel.grid = element_blank(),
                     axis.title = element_text(color = "black", size = 12),
                     axis.ticks = element_line(),
                     axis.text = element_text(color = "black", size = 10, face = "bold"),
                     axis.line = element_blank(),
                     plot.title=element_text(size=20),
                     plot.subtitle=element_text(size=15),
                     legend.text=element_text(size=12),
                     legend.title = element_text(size = 12),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15))+
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0))
  
saveRDS (metricsframe2, "Example_benchmark_data")
