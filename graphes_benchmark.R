library(ggplot2)
library(dplyr)

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

# New plots fixed dist ####
  ## Changer les labels pour les distances
  dist.labs <- c("Dist 3", "Dist 4", "Dist 5", "Dist 7")
  names(dist.labs) <- c("3", "4", "5", "7")
  
  # Changer les labels pour le nb de feuilles
  move.labs <- c("1 Leaf Moved", "3 Leaves Moved","5 Leaves Moved","10 Leaves Moved",
                 "15 Leaves Moved", "20 Leaves moved", "25 Leaves Moved")
  names(move.labs) <- c("1","3", "5", "10", "15", "20", "25")
  
  ## Nombre "d'erreurs" par arbre  ----
  plot3b = ggplot(data = metricsframe2) +
    geom_bar(aes(x = factor(FN))) +
    facet_grid(dist~nb.move, labeller=labeller(dist = dist.labs, nb.move = move.labs)) +
    labs(title = "Nombre d'abre en fonction du nombre de FN" ,
         x = "Nombre de feuilles déplacées retenues",
         y = "Nombre d'arbres observés")+
    scale_x_discrete(labels= c("0","","","","","5","","","","","10",
                               "","","","","15","","","","","20",
                               "","","","","25"))
  plot3b

  ## Matrice de couleur
  pamal = metricsframe2 %>% group_by(nbLeaves, nb.move, dist) %>%
    count(perfect, FN) %>%
    reframe(matval = FN*n) %>%
    group_by(nbLeaves, nb.move, dist) %>%
    tally(matval)
  
  ggplot(pamal, aes(x = dist, y = nb.move, fill = n, label = n)) +
    geom_tile() +
    geom_text(col = "rosybrown1") +
    scale_fill_viridis_c(begin = 0.2, end = 0.6, direction = -1, option = "plasma")+
    labs (title = "Feuilles erronnées retenues en fonction de la distance de déplacement",
          x = "Distance de déplacement",
          y = "Nombre de feuilles déplacées",
          subtitle = paste0("20 arbres de 40 feuilles par case"))+
    theme_bw() + theme(panel.background = element_blank(),
                       panel.border = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_text(color = "black", size = 12),
                       axis.ticks = element_line(),
                       axis.text = element_text(color = "black", size = 10, face = "bold"),
                       axis.line = element_blank())+
    scale_y_discrete(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0))
  
  ## Matrice de couleur "normalisée"
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
    geom_text(aes(x = dist, y = nb.move, label = compte_faux), col = "white", nudge_y = -0.3, size = 4) + 
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
  
  ## Nombre d'erreurs par distance entre feuilles déplacées ----
  tb = metricsframe2 %>% filter (FN >= 1) # Garder que les cas avec erreurs retenues
  
  # Déplier le vecteur de distance pour le mettre dans un tableau
  plot5bframe = data.frame()
  for (i in 1:nrow(tb)){
    tmp = data.frame()
    # Récup distances de branche entre les feuilles erronées retenues
    edgeofWrong = tb$edgeDist[[i]][match(tb$FNtips[[i]], tb$wrongTips[[i]])]
    
    for (j in 1:length(edgeofWrong)){
      tmp[j,1] = tb$id[i]
      tmp[j,2] = tb$dist[i]
      tmp[j,3] = tb$nb.move[i]
      tmp[j,4] = tb$FN[i]
      tmp[j,5] = edgeofWrong[j]
    }
    plot5bframe = rbind(plot5bframe, tmp)
  }
  
  colnames(plot5bframe) <- c("id","dist","nb.move","FN","edgeofWrong")
  
  plot5b = ggplot(data = plot5bframe) +
    geom_bar(aes(x = factor(FN), fill = factor(edgeofWrong))) +
    facet_grid(dist~nb.move, labeller=labeller(dist = dist.labs, nb.move = move.labs)) +
    labs(title = "Nombre de 'mauvaises' branche, colorées par distance",
         x = "Nombre de feuilles déplacées retenues",
         y = "#branches des mauvaises feuilles retenues") +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank())
  plot5b
  

# Tests ####


#' *patchwork package pour les insets*
n <- d + geom_bar(aes(fill = fl))
n + scale_fill_manual(
  values = c("skyblue", "royalblue", "blue", "navy"), 
  limits = c("d", "e", "p", "r"), breaks =c("d", "e", "p", "r"), 
  name = "fuel", labels = c("D", "E", "P", "R"))
n
