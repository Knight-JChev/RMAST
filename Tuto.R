library("ape")
plot(bird.families)
is.binary(bird.families)
plot(multi2di(bird.families))

# Tuto ####
foo <- function() { 
  col <- "green"
  for (i in 1:2){
    axis(i, col = col, col.ticks = col, col.axis = col, las = 1)}
  box(lty = "19")}

## Plotter des arbres ####

tr <- compute.brlen(stree(8, "l"), 0.1) 
tr$tip.label[] <- ""
layout(matrix(1:6, 3, 2, byrow=T))
par(las = 1)
plot(tr, "u"); foo() 
plot(tr, type ="u", use.edge.length = FALSE); foo()
plot(tr, "f"); foo()
plot(tr, "f", FALSE); foo()
plot(tr, "p"); foo()
plot(tr, "c"); foo()

## Pimper des arbres ####
mytr <- read.tree(text = "((Pan:5,Homo:5):2,Gorilla:7);")
plot(mytr); foo()

### Pimper les feuilles ####
par(mfrow=c(1,1), xpd = TRUE)
geo <- factor(c("Africa", "World", "Africa"))
(mycol <- c("blue", "red")[geo])
plot(mytr, tip.color = mycol, cex = c(1, 1, 1.5))

### Pimper les branches ####
(i <- which.edge(mytr, c("Homo", "Pan")))
co = rep("black", Nedge(mytr))
co[i] = "blue"
plot(mytr, edge.col = co)
edgelabels()

### Ajouter des labels ####
plot(mytr)
add.scale.bar()
axisPhylo()
axis(3)
nodelabels()
tiplabels()
edgelabels()
mytr$edge

### Ajouter des plots ####
x <- matrix(1:6, 3, 10) 
dimnames(x) <- list(c("Homo", "Gorilla", "Pan"), LETTERS[1:ncol(x)])
par(mar = c(10, 2, 5, 5)) 
plot(mytr, x.lim = 30) 
phydataplot(x, mytr, "m", border = "white", offset = 3, width = 1) 
phydataplot(x, mytr, "m", border = "white", offset = 15,
            width = 1, continuous = 2)

n <- 100 
p <- 3
set.seed(3) 
tr <- rcoal(n) 
x <- matrix(runif(p * n), n, p)

rownames(x) <- tr$tip.label
COL <- c("red", "green", "blue")
par(xpd = T)
plot(tr, "f", x.lim = c(-5, 5), show.tip.label = FALSE, no.margin = TRUE)
for (j in 1:p) ring(x[, j], tr, offset = j - 1 + 0.1, col = COL[j])

### Modifier des infos dans l'arbre ####
tree$edge.length[8] <- 1000

### Exemple des infos dans l'arbre 
plot(rtree(5))
(pp <- get("last_plot.phylo", envir = .PlotPhyloEnv))
nodelabels()
tiplabels()
mytr$edge
