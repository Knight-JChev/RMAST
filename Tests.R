library(ape)
ori = rtree(250, rooted = FALSE)
(q1 = summary(ori$edge.length)[2])
ttt = sample(1:(ori$Nnode-1), size = 2, replace = F)
tip1=ttt[1]
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
tiplabels(taxo$tip.label[tip2-1], tip2-1, adj=0, bg = "lightblue")

phylo = drop.tip(ori, tip2) # Enlever une feuille
phylo$tip.label[tip1] = ori$tip.label[tip2] # Changer le nom de la feuille qui est "déplacée"
plot(phylo, main = "Arbre Phylo", cex = 1.2)
nodelabels(adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
tiplabels(adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
tiplabels(phylo$tip.label[tip1],tip1, adj=0, bg = "orchid")

# Comparaison arbre taxo et phylo
compare = comparePhylo(taxo, phylo, plot = T, force.rooted = TRUE)
a=all.equal.phylo(phylo, taxo, use.edge.length = F, index.return = T)

# Fonctions qui peuvent être intéressantes
nodepath()
node.depth()
mrca()
as.bitsplits()
reorder.phylo ()
which.edge()
write.tree(phylo) # possibilité de sortir la valeur en console
getAnywhere(comparePhylo())


deparse(substitute(phylo))
Ntip(phylo)

# Tests perso ####
key1 <- makeNodeLabel(taxo, "md5sum")$node.label
key2 <- makeNodeLabel(phylo, "md5sum")$node.label
mk12 <- match(key1, key2)
mk21 <- match(key2, key1)
tmp <- is.na(mk12)

## Choper les différences ####
noeud_different <- which(tmp)+Ntip(taxo) # numéros des noeuds différents
MRCA = getMRCA(taxo, noeud_different) # numéro du noeud où la divergence commence
le_plus_a_droite = noeud_different[  
  which.min(node.depth(taxo)[noeud_different])] #numéro du noeud le plus "bas" dans la phylo


def.par <- par(no.readonly = TRUE)
layout(matrix(1:2, 1, 2))
plot(phylo, use.edge.length = F, main = "taxo")
nodelabels(node = which(tmp)+Ntip(taxo), pch = 19, col = "blue", cex = 2)
nodelabels()
tiplabels()
legend(location, legend = paste("Clade absent in", tree2), pch = 19, col = "blue")

nodepath(taxo, from=MRCA, to=le_plus_a_droite)

# Prendre et manipuler des sous arbres ####
phyloNodes = phylo$node.label
phyloNodesNum = unique(phylo$edge[,1]) # liste des numéro des noeuds internes
phyloLeaf = phylo$tip.label # liste des noms des feuilles

df = as.data.frame(phylo$edge)
df = df %>%
  rename("from"="V1","to"="V2") %>%
  mutate(from = recode(from, !!!setNames(phyloNodes, phyloNodesNum))) %>%
  mutate(to = recode(to, !!!setNames(phyloNodes, phyloNodesNum))) %>%
  mutate(to = recode(to, !!!setNames(phyloLeaf, phyloLeafNum)))
