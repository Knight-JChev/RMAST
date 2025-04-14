library(ggplot2)
library(dplyr)

dir = "Benchmark_data/Data_20rep_40leaves_dist/"
metricsframe = data.frame()
for (file in list.files(dir)){
  if (!(dir.exists(paste0(dir,file)))){
    df = readRDS(paste0(dir,file))
    metricsframe = rbind(metricsframe, df)
  }
}
# New plots random Tax ####
  ## Distribution des tailles d'arbres ----
  plot2 = ggplot(data = metricsframe) +
    geom_dotplot(aes(x = nb.move, y = obsTreeSize), fill = metricsframe$perfect, 
                 binwidth = 1, dotsize = 0.5, binaxis = "y",
                 stackdir = "center", stackratio = 1.2) +
    labs(title = "Distribution des tailles d'arbres" ,
         x = "Nombre d'erreurs demandées",
         y = "Taille observée de l'arbre") +
    theme_bw() + theme(panel.grid.major = element_line(colour = "grey80", linewidth = 0.5),
                       panel.grid.minor.y = element_line(colour = "grey", linewidth =0.2),
                       panel.grid.major.y = element_line(colour = "grey", linewidth =0.5),
                       panel.grid.minor.x=element_blank()) +
    scale_y_continuous(minor_breaks = seq(15, 25, by = 1)) +
    coord_cartesian(ylim = c(0,25))
  plot2
  
  plot2b = ggplot(data = metricsframe) +
    geom_dotplot(aes(x = nb.move, y = obsTreeSize, fill = perfect),
                 binwidth = 1, dotsize = 0.5, binaxis = "y",
                 stackdir = "center", stackratio = 1, position = "dodge") +
    labs(title = "Distribution des tailles d'arbres" ,
         x = "Nombre d'erreurs demandées",
         y = "Taille observée de l'arbre") +
    theme_bw() + theme(panel.grid.major = element_line(colour = "grey80", linewidth = 0.5),
                       panel.grid.minor.y = element_line(colour = "grey", linewidth =0.2),
                       panel.grid.major.y = element_line(colour = "grey", linewidth =0.5),
                       panel.grid.minor.x=element_blank()) +
    scale_y_continuous(minor_breaks = seq(15, 25, by = 1)) +
    coord_cartesian(ylim = c(0,25))
  plot2b
  
  
  ## Nombre "d'erreurs" par arbre  ----
  plot3 = ggplot(data = metricsframe) +
    geom_bar(aes(x = factor(FN))) +
    facet_grid(.~nb.move) +
    labs(title = "Nombre d'abre en fonction du nombre de FN" ,
         x = "Nombre de feuilles déplacées retenues",
         y = "Nombre d'arbres observés")
  plot3
  
  ## Créer un tableau temporaire en allongeant un vecteur en colonne
  plot4frame = data.frame()
  for (i in 1:nrow(metricsframe)){
    tmp = data.frame()
    for (j in 1:length(metricsframe$edgeDist[[i]])){
      tmp[j,1] = metricsframe$id[i]
      tmp[j,2] = metricsframe$nb.move[i]
      tmp[j,3] = metricsframe$FN[i]
      tmp[j,4] = metricsframe$edgeDist[[i]][j]
    }
    plot4frame = rbind(plot4frame, tmp)
  }
  colnames(plot4frame) <- c("id","nb.move","FN","edgeDist")
  
  plot4 = ggplot(data = plot4frame) +
    geom_bar(aes(x = factor(FN), fill = factor(edgeDist))) +
    facet_grid(.~nb.move) +
    labs(title = "Branches entre feuilles déplacées, coloré par distance" ,
         x = "Nombre de feuilles déplacées retenues",
         y = "Nombre de branches entre les feuilles déplacées")
  plot4
  
  ## ---------
  t = metricsframe %>% filter (FN >= 1)
  
  plot5frame = data.frame()
  for (i in 1:nrow(t)){
    tmp = data.frame()
    edgeofWrong = t$edgeDist[[i]][match(t$FNtips[[i]], t$wrongTips[[i]])]
    
    for (j in 1:length(edgeofWrong)){
      tmp[j,1] = t$id[i]
      tmp[j,2] = t$nb.move[i]
      tmp[j,3] = t$FN[i]
      tmp[j,4] = edgeofWrong[j]
    }
    plot5frame = rbind(plot5frame, tmp)
  }
  
  colnames(plot5frame) <- c("id","nb.move","FN","edgeofWrong")
  
  plot5 = ggplot(data = plot5frame) +
    geom_bar(aes(x = factor(FN), fill = factor(edgeofWrong))) +
    facet_grid(.~nb.move) +
    labs(title = "Nombre de 'mauvaises' branche, colorées par distance",
         x = "Nombre de feuilles déplacées retenues",
         y = "#branches des mauvaises feuilles retenues") +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank())+
    scale_y_continuous(breaks = seq(1,10, by=1))
  plot5

# New plots fixed dist ####
  ## Changer les labels pour les distances
  dist.labs <- c("Dist 3", "Dist 4", "Dist 5")
  names(dist.labs) <- c("3", "4", "5")
  
  # Changer les labels pour le nb de feuilles
  move.labs <- c("1 Leaf Moved", "5 Leaves Moved", "10 Leaves moved")
  names(move.labs) <- c("1", "5", "10")
  
  ## Nombre "d'erreurs" par arbre  ----
  plot3b = ggplot(data = metricsframe) +
    geom_bar(aes(x = factor(FN))) +
    facet_grid(dist~nb.move, labeller=labeller(dist = dist.labs, nb.move = move.labs)) +
    labs(title = "Nombre d'abre en fonction du nombre de FN" ,
         x = "Nombre de feuilles déplacées retenues",
         y = "Nombre d'arbres observés")
  plot3b

  ## Matrice de couleur
  pamal = metricsframe %>% group_by(nbLeaves, nb.move, dist) %>%
    count(perfect, FN) %>%
    reframe(matval = FN*n) %>%
    group_by(nbLeaves, nb.move, dist) %>%
    tally(matval)
  
  ggplot(pamal, aes(x = dist, y = nb.move, fill = n, label = n)) +
    geom_tile() +
    geom_text(col = "pink") +
    scale_fill_viridis_c(begin = 0.2, end = 0.6, direction = -1, option = "plasma")+
    theme_bw()
  
  ## Nombre d'erreurs par distance entre feuilles déplacées ----
  tb = metricsframe %>% filter (FN >= 1) # Garder que les cas avec erreurs retenues
  
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
d <- ggplot(mpg, aes(fl))
d + geom_bar()

df <- data.frame(grp = c("A", "B"), fit = 4:5, se = 1:2)
j <- ggplot(df, aes(grp, fit, ymin = fit - se, ymax = fit + se))

f <- ggplot(mpg, aes(class, hwy))
f + geom_dotplot(binaxis = "y", stackdir = "center") 

g <- ggplot(diamonds, aes(cut, color))
g + geom_count()

#' *patchwork package pour les insets*
n <- d + geom_bar(aes(fill = fl))
n + scale_fill_manual(
  values = c("skyblue", "royalblue", "blue", "navy"), 
  limits = c("d", "e", "p", "r"), breaks =c("d", "e", "p", "r"), 
  name = "fuel", labels = c("D", "E", "P", "R"))
n
