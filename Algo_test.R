library(ape)
library(dplyr)

createTaxPhy <- function(NbLeaf = 10){
  world = rtree(NbLeaf, rooted = FALSE)
  world = makeNodeLabel(world, method = "number", prefix = "Node")
  (q1 = summary(world$edge.length)[2])
  ttt = sample(2:(world$Nnode-1), size = 2, replace = F)
  tip1 = ttt[1]
  tip2 = ttt[2]
  
  
  # Dichotomie vers polytomie ####
  #' [Dichotomie vers polytomie]
  par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
  layout(matrix(c(1,2),1,2)) # Matrice pour tracer les plots
  world = di2multi(world, tol= round(q1, digits = 1)) # Binaire vers polytomies
  
  # Création arbre taxo et phylo ####
  #' [Visu taxo et phylo]
  layout(matrix(c(1,2,3),1,3)) # Matrice pour tracer les plots
  
  plot(world, cex = 1.2, main = "World")
  nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  tiplabels(world$tip.label[tip1], tip1, adj=0)
  tiplabels(world$tip.label[tip2], tip2, adj=0, bg = "lightblue")
  
  taxo = drop.tip(world, tip1) # Enlever une feuille
  tipToChange = which(taxo$tip.label==world$tip.label[tip2])
  plot(taxo, main = "Taxo", cex = 1.2)
  nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  tiplabels(world$tip.label[tip2], tipToChange, adj=0, bg = "lightblue")
  
  phylo = drop.tip(world, tip2) # Enlever une feuille
  tipToChange = which(phylo$tip.label==world$tip.label[tip1])
  phylo$tip.label[tipToChange] = world$tip.label[tip2] # Changer le nom de la feuille qui est "déplacée"
  plot(phylo, main = "Phylo", cex = 1.2)
  nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  tiplabels(world$tip.label[tip2], tipToChange, adj=0, bg = "orchid")

  return(c(taxo, phylo, world))
}

createTaxPhy(11)

 
