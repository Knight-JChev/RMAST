library(ape)
library(dplyr)
par(xpd = TRUE, mar = c(2, 2, 0.3, 0))

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

a = read.tree(text = "(((Mammifère1:1.0):1.0,(Mammifère2:1.0,Mammifère3:1.0):1.0):1.0,((Plante5:1.0,Plante6:1.0):1.0,(Plante7:1.0,Humain: 1.0):1.0):1.0);")
b = read.tree(text = "(((Humain: 1.0,Mammifère1:1.0):1.0,(Mammifère2:1.0,Mammifère3:1.0):1.0):1.0,((Plante5:1.0,Plante6:1.0):1.0,(Plante7:1.0):1.0):1.0);")
c = read.tree(text = "(((Mammifère1:1.0):1.0,(Mammifère2:1.0,Mammifère3:1.0):1.0):1.0,((Plante5:1.0,Plante6:1.0):1.0,(Plante7:1.0):1.0):1.0);")

plot.phylo(a, show.node.label = F, cex = 2.5)
plot(b, cex = 2.5)
plot(c, cex = 2.5)

# Exemple MAST KAST ####
# Arbres de départ

layout(matrix(c(1),1,1))
a = read.tree(text = "(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1);")
b = read.tree(text = "(((a:1,c:1):1,b:1):1,((d:1,e:1):1,f:1):1);")

png("./graphes/kast_a1.png", width = 350, height = 550)
plot(a, direction = "downwards")
makeNodeLabel(phy = a, method = c("1","2","3","4","5"))
tiplabels(a$tip.label, cex = 2.2, frame = "c", col = c("black"), 
          bg = c("cyan2","lemonchiffon","mediumpurple",rep("coral1",3)))
nodelabels(a$node)
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

