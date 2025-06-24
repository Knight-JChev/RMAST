library(ape)
library(dplyr)
par(xpd = TRUE, mar = c(2, 2, 0.3, 0))


# Prend en compte la distance de déplacement des feuilles
createTaxPhyAlt <- function(nbLeaves = 10, nb.move = 2, dist = 4){
  
  # Vérif 
  if (dist <= 2 || dist > nbLeaves/2){
    return(print("Argument 'dist', has to be between 2 and half the number of leaves"))
  }
  # Création de l'arbre taxo
  taxo = rtree(nbLeaves, rooted = FALSE)
  taxo$edge.length[] <- 1 # Toutes les arrêtes font 1 de long
  
  # Matrice de distance et choix des feuilles ####
  # Récupérer des couples à une certaine distance 
  # Peux pas prendre 2 fois la même feuille
  distMat = cophenetic.phylo(taxo) 
  distMat[lower.tri(distMat)] <- NA
  tmpcoords = data.frame(which(distMat == dist, arr.ind = T)) # Coordonnées OK
  # Enlever doubles dans ligne; dans colonnes; entre lignes et colonnes
  tmpcoords = distinct(tmpcoords, row, .keep_all = T) 
  tmpcoords = distinct(tmpcoords, col, .keep_all = T) 
  tmpcoords = tmpcoords[-(na.omit(match(tmpcoords$row, tmpcoords$col))),]
  
  # Choix de le feuille à déplacer et de la feuille où déplacer
  tipMoved = c()
  tipTo = c()
  i = 1
  if (nrow(tmpcoords) >= nb.move){ # Au moins autant de résultats que demandé
    while (i <= nb.move){
      tipTo = append(tipTo, tmpcoords[i,1])
      tipMoved = append(tipMoved, tmpcoords[i,2])
      i = i+1
    }
  } else return(createTaxPhyAlt(nbLeaves, nb.move, dist)) # Sinon, relancer
  
  # Formatage pour benchmark
  tipsNumber = sub(".", "", taxo$tip.label)
  wrongTips = taxo$tip.label[tipMoved]
  
  # Création arbre phylo ####
  #' [Visu taxo et phylo]
  par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
  taxo = makeNodeLabel(taxo, method = "number", prefix = "tax") 
  
  layout(matrix(c(1,2),1,2)) # Matrice pour tracer les plots
  
  # Plot arbre taxo avec les bons nom
  plot(taxo, cex = 1.5, font = 2)
  #nodelabels(taxo$node.label, adj = c(1,-0.2), frame = "n", cex = 0.8, font = 2, col="red")
  tiplabels(taxo$tip.label[tipTo], tipTo, frame = "c",
            adj=0, font = 2, cex = 1.5, bg = "mediumpurple1")
  tiplabels(wrongTips, tipMoved, frame = "c",
            adj=0, bg = "lightblue", font = 2, cex = 1.5)
  
  # Création arbre phylo
  phylo = taxo
  tipToAdjust = c()
  for (i in 1:length(tipMoved)){
    # Ajuste l'indice des feuilles dans phylo comme on en enlève une par une
    tmptree = keep.tip(taxo, taxo$tip.label[tipMoved[i]])
    phylo = drop.tip(phylo, taxo$tip.label[tipMoved[i]])
    tipToAdjust = append(tipToAdjust,
                         which(phylo$tip.label == taxo$tip.label[tipTo[i]])) 
    phylo = bind.tree(phylo, tmptree, where = tipToAdjust[i], position =  0.5)
  }
  phylo$edge.length[] <- 1
  
  # Plot arbre phylo
  tipMovedTo = match(taxo$tip.label[tipMoved], phylo$tip.label)
  phylo = makeNodeLabel(phylo, method = "number", prefix = "phy")
  plot(phylo, cex = 1.5, font = 2)
  #nodelabels(phylo$node.label, adj = c(1,-0.2), frame = "n", cex = 0.8, font = 2, col="red")
  tiplabels(taxo$tip.label[tipTo], tipToAdjust, frame = "c",
            adj=0, bg = "mediumpurple1", font = 2, cex = 1.5)
  tiplabels(taxo$tip.label[tipMoved], tipMovedTo,
            adj=0, frame = "c", bg = "lightblue", font = 2, cex = 1.5)
  
  
  # Nombre d'arrêtes entre les paires de feuilles déplacées ####
  edgeDist = c()
  notWrongTips = c()
  wrongTips = sub(".","",wrongTips)
  
  for (i in 1:nb.move) {
    edges = (length(nodepath(taxo, tipMoved[i], tipTo[i])))-1
    if (edges == 2) {
      notWrongTips = append(notWrongTips, tipsNumber[tipMoved[i]])
    }
    edgeDist = append(edgeDist, edges)
  }
  
  if (!(identical(pmatch(notWrongTips, wrongTips), integer(0)))){
    realWrongTips = wrongTips[-pmatch(notWrongTips, wrongTips)]
  } else realWrongTips = wrongTips
  
  
  taxo[[6]] = "taxo"
  names(taxo)[6] = "name"
  
  phylo[[6]] = "phylo"
  names(phylo)[6] = "name"
  
  return(list(taxo = taxo, phylo = phylo, 
              notWrongTips = notWrongTips, truth = c(tipsNumber[-c(tipMoved,tipTo)], notWrongTips),
              edgeDist = edgeDist, realWrongTips = realWrongTips, wrongTips = wrongTips))
}
createTaxPhyAlt(nbLeaves = 10, nb.move = 3, dist = 3)

# Arbre équilibré ####
layout(matrix(c(1),1,1))
a = stree(8, "b")
a$edge.length = rep(1, length(a$tip.label)*2-2)
a$node.label = c("a","b","d","e","c","f","g")

plot.phylo(a, "c", show.node.label = F, show.tip.label = F)
nodelabels(a$node.label, adj = c(1,-0.2), frame = "n", cex = 0.8, font = 2, col="red")
ape::write.tree(a)

# Exemple humain:plante ####
layout(matrix(c(1,2,3),1,3))

a = read.tree(text = "(((Plante1:1,Plante2:1)d:1,(Plante3:1,Plante4:1)e:1)b:1,((Plante5:1,Plante6:1)f:1,(Plante7:1,Humain:1)g:1)c:1)a;")
b = read.tree(text = "(((Plante1:1,Plante2:1)d:1,(Plante3:1,Plante4:1)e:1)b:1,((Plante5:1,Plante6:1)f:1,(Plante7:1)g:1)c:1,Humain:3)a;")
c = read.tree(text = "(((Plante1:1,Plante2:1)d:1,(Plante3:1,Plante4:1)e:1)b:1,((Plante5:1,Plante6:1)f:1,(Plante7:1)g:1)c:1)a;")

plot.phylo(a, show.node.label = F, cex = 2.5)
plot(b, cex = 2.5)
plot(c, cex = 2.5)

# Exemple MAST KAST ####
# Arbres de départ

layout(matrix(c(1),1,1))
a = read.tree(text = "(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1);")
b = read.tree(text = "(((a:1,c:1):1,b:1):1,((d:1,e:1):1,f:1):1);")

png("./graphes/kast_a1.png", width = 350, height = 550)
plot(a)
tiplabels(a$tip.label, cex = 2.2, frame = "c", col = c("black"), 
          bg = c("cyan2","lemonchiffon","mediumpurple",rep("coral1",3)) , srt=90)
dev.off()

png("./graphes/kast_b1.png", width = 350, height = 550)
plot(b)
tiplabels(b$tip.label, cex = 2.2, frame = "c", col = c("black"), 
          bg = c("cyan2","mediumpurple","lemonchiffon",rep("coral1",3)) , srt=90)
dev.off()

# MASTS
#layout(matrix(c(1,2,3),1,3))
m1 = read.tree(text = "((a:1,b:1):1,((d:1,e:1):1,f:1):1);")
m2 = read.tree(text = "((a:1,c:1):1,((d:1,e:1):1,f:1):1);")
m3 = read.tree(text = "((b:1,c:1):1,((d:1,e:1):1,f:1):1);")

png("./graphes/kast_m1.png", width = 350, height = 550)
plot(m1)
tiplabels(m1$tip.label, cex = 2.2, frame = "c", col = "black", 
          bg = c("cyan2","lemonchiffon",rep("coral1",3)) , srt=90)
dev.off()

png("./graphes/kast_m2.png", width = 350, height = 550)
plot(m2)
tiplabels(m2$tip.label, cex = 2.2, frame = "c", col = "black", 
          bg = c("cyan2","mediumpurple",rep("coral1",3)) , srt=90)
dev.off()

png("./graphes/kast_m3.png", width = 350, height = 550)
plot(m3)
tiplabels(m3$tip.label, cex = 2.2, frame = "c", col = "black", 
          bg = c("lemonchiffon","mediumpurple",rep("coral1",3)) , srt=90)
dev.off()

# KAST
layout(matrix(c(1),1,1))
k1 = read.tree(text = "((d:1,e:1):1,f:1):1;")

png("./graphes/kast_k1.png", width = 350, height = 550)
plot(k1)
tiplabels(k1$tip.label, cex = 2.2, frame = "c", col = c("black"), 
          bg = c("coral1") , srt=90)
dev.off()

# Exemple MAST DIAG et MAST matching ####
layout(matrix(c(1,2),1,2))
a = read.tree(text = "(((t3:0.5,t4:0.5)a1:0.5,(t1:0.5,t6:0.5)b1:0.5)d1:0.5,((t5:0.5,t2:0.5)c1:0.5)f0:0.5)top;")
plot.phylo(a, show.node.label = F, show.tip.label = F)
# Plots high nodes in circles
nodelabels(a$node.label[2], node = 2+6, col="black", 
           adj = c(0.5,0.5), frame = "c", cex = 1.6, font = 2, 
           bg=c("coral1"),srt = 90)
# Plot low nodes
nodelabels(a$node.label[c(3,4)], node = c(3,4)+6, col="black",
           adj = c(0.5,0.5), frame = "c", cex = 1.4, font = 2, 
           bg=c(rep("coral1",2)),srt = 90)
# Plpot leaves
tiplabels(a$tip.label, frame = "c", cex = 1.2, font = 2, srt = 90, bg = "white")
tiplabels(a$tip.label[c(3,5)], tip = c(3,5), frame = "c", cex = 1.2, font = 2, 
          srt = 90, bg = c("springgreen4","olivedrab3"), col = c("grey90","black"))


b = read.tree(text = "(((t1:0.5,t2:0.5)c1:1,((t3:0.5,t4:0.5)a0:0.5,(t5:0.5,t6:0.5)b0:0.5)d0:0.5)f1:0)top;")
plot.phylo(b, show.node.label = F, show.tip.label = F)
# Plots high nodes in circles
nodelabels(b$node.label[4], node = 4+6, col="black", 
           adj = c(0.5,0.5), frame = "c", cex = 1.6, font = 2, 
           bg=c("aquamarine"),srt = 90)
# Plot low nodes
nodelabels(b$node.label[c(5,6)], node = c(5,6)+6, col="black",
           adj = c(0.5,0.5), frame = "c", cex = 1.4, font = 2, 
           bg=c(rep("aquamarine",2)),srt = 90)
# Plot leaves
tiplabels(b$tip.label, frame = "c", cex = 1.2, font = 2, srt = 90,  bg = "white")
tiplabels(b$tip.label[c(5,1)], tip = c(5,1), frame = "c", cex = 1.2, font = 2, 
          srt = 90, bg = c("olivedrab3","springgreen4"), col = c("black","grey90"))

# Exemple arbre et tableau ####
layout(matrix(c(1,2),1,2))
a = stree(4, "b", tip.label = c(1,2,3,4))
a$edge.length = rep(1, length(a$tip.label)*2-2)
a$node.label = c(5,6,7)

plot.phylo(a, "c", show.node.label = F, show.tip.label = F, direction = "downwards")
tiplabels(a$tip.label, cex = 1.5, frame = "c", col = "black", bg = "coral1" )
nodelabels(a$node.label, frame = "c", cex = 1.5, font = 2, col="black", bg = "aquamarine")

# Oui ####
x = rtree(12)
x$edge.length[] <- 1
plot(x, show.tip.label = F)

