library(ape)
library(dplyr)

# Initialisation de l'arbre ####
ori = rtree(11, rooted = FALSE)
(q1 = summary(ori$edge.length)[2])

# Visualisation de l'arbre ####
par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
plot(ori)
nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  
writeNewick <-function(tree){  
  # Métriques et repères sur l'arbre
  nbFeuilles = length(tree$tip.label)
  numNoeuds = Nnode(tree) + nbFeuilles
  
  # Data frame des edges
  edgePaths = as.data.frame(tree$edge)
  edgePaths = edgePaths %>%
    rename("from"="V1","to"="V2")
  
  root = nbFeuilles + 1
  return (newick(root, nbFeuilles, edgePaths))
}

newick <- function (noeud, nbFeuilles, edgePaths){
     current = filter(edgePaths, from == noeud)
     nwk = ""
     
       for (i in current[,2]){

         # Si l'enfant est le premier
         if (i == current[1, 2])  { 
           
           if (i %in% 1:nbFeuilles)  nwk = paste0( "(", nwk, i, ",")
           else nwk = paste0( "(", nwk, newick(i, nbFeuilles, edgePaths), ",")
         
        # Si l'enfant est le dernier
         } else if (i == current[length(current[,2]), 2])  { 
           
           if (i %in% 1:nbFeuilles)  nwk = paste0(nwk,i,"):", noeud)
           else nwk = paste0(nwk,newick(i, nbFeuilles, edgePaths),"):", noeud)
         
         # Si l'enfant est au milieu  
         } else {
           
           if (i %in% 1:nbFeuilles)  nwk = paste0(nwk,i,",")
           else nwk = paste0(nwk,newick(i, nbFeuilles, edgePaths),",")
           
         }  
              }
     return(nwk)
  }

writeNewick(ori)



