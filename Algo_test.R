library(ape)
library(dplyr)
ori = rtree(11, rooted = FALSE)
ori = makeNodeLabel(ori, method = "number", prefix = "Node")
(q1 = summary(ori$edge.length)[2])
ttt = sample(2:(ori$Nnode-1), size = 2, replace = F)
tip1 = ttt[1]
tip2 = ttt[2]


# Avant après polytomie ####
#' [Avant après polytomie]
summary(ori)
par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
layout(matrix(c(1,2),1,2)) # Matrice pour tracer les plots
plot (ori)
ori = di2multi(ori, tol= round(q1, digits = 1)) # Binaire vers polytomies
plot(ori)

# Création arbre taxo et phylo ####
#' [Visu taxo et phylo]
layout(matrix(c(1,2,3),1,3)) # Matrice pour tracer les plots

plot(ori, cex = 1.2, main = "Arbre monde")
nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
tiplabels(ori$tip.label[tip1], tip1, adj=0)
tiplabels(ori$tip.label[tip2], tip2, adj=0, bg = "lightblue")

taxo = drop.tip(ori, tip1) # Enlever une feuille
plot(taxo, main = "Arbre Taxo", cex = 1.2)
nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
tipToChange = which(taxo$tip.label==ori$tip.label[tip2])
tiplabels(ori$tip.label[tip2], tipToChange, adj=0, bg = "lightblue")

phylo = drop.tip(ori, tip2) # Enlever une feuille
plot(phylo, main = "Arbre Phylo", cex = 1.2)
nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
tipToChange = which(phylo$tip.label==ori$tip.label[tip1])
phylo$tip.label[tipToChange] = ori$tip.label[tip2] # Changer le nom de la feuille qui est "déplacée"
tiplabels(ori$tip.label[tip2], tipToChange, adj=0, bg = "orchid")


# Prendre et manipuler des sous arbres ####
test = extract.clade (taxo,) # choisir un sous arbre

phyloNodes = phylo$node.label
phyloNodesNum = unique(phylo$edge[,1]) # liste des numéro des noeuds internes
phyloLeaf = phylo$tip.label # liste des noms des feuilles

df = as.data.frame(phylo$edge)

df = df %>%
  rename("from"="V1","to"="V2") %>%
  mutate(from = recode(from, !!!setNames(phyloNodes, phyloNodesNum))) %>%
  mutate(to = recode(to, !!!setNames(phyloNodes, phyloNodesNum))) %>%
  mutate(to = recode(to, !!!setNames(phyloLeaf, phyloLeafNum)))

node.depth(phylo)

parcours <- function(arbre){
  
}

taxoNodes = taxo$node.label
taxoNodesNum = unique(taxo$edge[,1]) 
taxoLinks = 

tablesave = data.frame()

keyPhylo <- makeNodeLabel(phylo, "md5sum")$node.label
keyTaxo <- makeNodeLabel(taxo, "md5sum")$node.label
matchPhyTax <- match(keyPhylo, keyTaxo)
matchTaxPhy <- match(keyTaxo, keyPhylo)







for (i in phyloNodes ){
  PhyloSubtree = extract.clade (phylo, i)
  
  for (j in taxoNodes){
    TaxoSubtree = extract.clade (taxo, j)
    
    if (PhyloSubtree == TaxoSubtree){
      print("MATCH : Phylo", i, "et Taxo", j)
    } else {
      print("MISS : Phylo", i, "et Taxo", j)
    }
  }
}
data.
