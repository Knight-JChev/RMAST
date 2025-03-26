library(ape)
library(dplyr)

# Initialisation de l'arbre ####
ori = rtree(11, rooted = FALSE)
(q1 = summary(ori$edge.length)[2])

# Visualisation de l'arbre
par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
plot(ori)
nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")

# Métriques et repères sur l'arbre
nbFeuilles = length(ori$tip.label)
numNoeuds = Nnode(ori) + nbFeuilles

# Data frame des edges
df = as.data.frame(ori$edge)
df = df %>%
  rename("from"="V1","to"="V2")

root = nbFeuilles + 1

newick(17, df)


newick <- function (noeud, df){
     current = filter(df, from == noeud)
     nwk = "("
     
       for (i in current[,2]){
           print(i)
        
         # Si l'enfant est le dernier
         if (i == current[length(current[,2]), 2]){ 
           
           if (i %in% 1:nbFeuilles)  nwk = paste0(nwk,i,")")
           else nwk = paste0(nwk,"(",newick(i, df),")")
         
         } else {
           
           if (i %in% 1:nbFeuilles)  nwk = paste0(nwk,i,",")
           else nwk = paste0(nwk,"(",newick(i, df),"),")
           
         }  
              }
     return(nwk)
  }

  