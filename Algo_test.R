library(ape)
library(dplyr)
library(RcppHungarian)

# Fonction pour faire la phylo et la taxo.
# Retourne un vecteur avec taxo, phylo, arbre d'origine
createTaxPhy <- function(NbLeaves = 10){
  world = rtree(NbLeaves, rooted = FALSE)
  (q1 = summary(world$edge.length)[2])
  ttt = sample(2:(world$Nnode-1), size = 2, replace = F)
  tip1 = ttt[1]
  tip2 = ttt[2]
  
  
  # Dichotomie vers polytomie ####
  #' [Dichotomie vers polytomie]
  par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
  layout(matrix(c(1,2),1,2)) # Matrice pour tracer les plots
  world = di2multi(world, tol= round(q1, digits = 1)) # Binaire vers polytomies
  world = makeNodeLabel(world, method = "number", prefix = "n")
  
  # Création arbre taxo et phylo ####
  #' [Visu taxo et phylo]
  layout(matrix(c(1,2,3),1,3)) # Matrice pour tracer les plots
  
  plot(world, cex = 1.2, main = "World")
  nodelabels(world$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(world$tip.label, adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  tiplabels(world$tip.label[tip1], tip1, adj=0)
  tiplabels(world$tip.label[tip2], tip2, adj=0, bg = "lightblue")
  
  taxo = drop.tip(world, tip1) # Enlever une feuille
  tipToChange = which(taxo$tip.label==world$tip.label[tip2])
  plot(taxo, main = "Taxo", cex = 1.2)
  nodelabels(taxo$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(taxo$tip.label, adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  tiplabels(world$tip.label[tip2], tipToChange, adj=0, bg = "lightblue")
  
  phylo = drop.tip(world, tip2) # Enlever une feuille
  tipToChange = which(phylo$tip.label==world$tip.label[tip1])
  phylo$tip.label[tipToChange] = world$tip.label[tip2] # Changer le nom de la feuille qui est "déplacée"
  plot(phylo, main = "Phylo", cex = 1.2)
  nodelabels(phylo$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(phylo$tip.label, adj = c(-2,0.3), frame = "n", cex = 1.5, font = 2, col="green")
  tiplabels(world$tip.label[tip2], tipToChange, adj=0, bg = "orchid")
  
  return(list(taxo = taxo, phylo = phylo, world = world))
}

# Fonction pour faire une dataframe avec le nom des feuilles
edgesToDf <- function(tree){
  
  nbNodes = (length(tree$tip.label)+Nnode(tree))
  newLabels = c(tree$tip.label, tree$node.label)
  
  edgePaths = as.data.frame(tree$edge)
  edgePaths = edgePaths %>%
    rename("from"="V1","to"="V2") %>%
    mutate(from = recode(from, !!!setNames(newLabels, as.character(1:nbNodes)))) %>%
    mutate(to = recode(to, !!!setNames(newLabels, as.character(1:nbNodes))))
  
  return(edgePaths)
}

# Métriques et repères sur l'arbre
treeMetrics <- function (T1, T2){
  nbT1Leaves = length(T1$tip.label)
  leavesT1 = T1$tip.label # liste des noms des feuilles
  rootT1 = T1$node.label[1]
  edgeT1 = edgesToDf(T1)
  
  nbT2Leaves = length(T2$tip.label)
  leavesT2 = T2$tip.label
  rootT2 = T2$node.label[1]
  edgeT2 = edgesToDf(T2)
  
  metrics = list(nbT1Leaves = nbT1Leaves,
                 leavesT1= leavesT1,
                 rootT1 = rootT1,
                 edgeT1 = edgeT1,
                 nbT2Leaves = nbT2Leaves,
                 leavesT2 = leavesT2,
                 rootT2 = rootT2,
                 edgeT2 = edgeT2)
  
  return(metrics)
}

# Mast en cas d'étages différents
maststep <- function(subrootT1, subrootT2, T1, T2, mastlist){ #On garde les arbres d'origine

  # Utiliser les subroot pour faire des sous arbres si la racine est différente
  if (subrootT1 != T1$node.label[1])
    subT1 = extract.clade(T1, subrootT1)
  if (subrootT2 != T2$node.label[1])
    subT2 = extract.clade(T2, subrootT2)
  
  # Calculer des métriques sur les sous arbres
  subTmetrics = treeMetrics(T1 = subT1, T2 = subT2)
  
  # Tableaux des paths du noeud en cours 
  currentT1Node = filter(subTmetrics$edgeT1, from == subTmetrics$rootT1) %>%
                    arrange(to)
  currentT2Node = filter(subTmetrics$edgeT2, from == subTmetrics$rootT2) %>%
                    arrange(to)
  
  # Boucle T1 sur sous-arbre de T2 ####
  
  # Si premier noeud de l'arbre, séparer les résultats des sous arbres
  if (subTmetrics$rootT2 == "n1"){
    for (firstnode in currentT2Node[,2]){
      if (firstnode %in% subTmetrics$leavesT1){
        mastlist = append(mastlist, firstnode)
      }
      else{
        mastlist = append(mastlist, maststep(subrootT1, firstnode, T1, T2, mastlist))
      }
    }
  }
  
  # Si on est dans un des sous arbre, faire le maststep récursivement
  else { 
   for (subnodeT2 in currentT2Node[,2]){ # i = sous-noeuds du noeud en cours
    print(paste0("Sous-noeud en cours ", subnodeT2))
    
    if (subnodeT2 %in% subTmetrics$leavesT2) {  # si le sous-noeud est une feuille

      if (subnodeT2 %in% subTmetrics$leavesT1){ # si cette feuille appartient à l'autre sous-arbre
        mastlist = append(mastlist, subnodeT2)
        print(paste0("Mastlist in ", (paste0(mastlist, collapse = " "))))
        }
    } else {
      # Ajoute la mastlist des enfants a celle du noeud en cours
      mastlist = append(mastlist, maststep(subrootT1, subnodeT2, T1, T2, mastlist))
    }
   }
    # Concaténer les feuilles
    mastlist = paste0(mastlist, collapse = "")
  }
  
  print(paste0("Mastlist end ", (paste0(mastlist, collapse = " "))))
  
  # Choisir le max de la mastlist
  return (mastlist[which.max(lapply(mastlist, nchar))])
}

maststep(subrootT1 = T1$node.label[1], subrootT2 = T2$node.label[1],
     T1 = T1, T2 = T2, mastlist = c())

# Création des arbres à tester
x = createTaxPhy(10)

T1 = x$taxo
T2 = x$phylo
subT1 = T1
subT2 = T2

mastunion <- function(subT1, subT2){
  # Mast/match en cas d'étages similaires
  subTmetrics = treeMetrics(T1 = subT1, T2 = subT2)
  
  # Tableaux des paths du noeud en cours, "to" ordonné noeuds avant feuilles
  currentT1Node = filter(subTmetrics$edgeT1, from == subTmetrics$rootT1) %>%
    arrange(to)
  currentT2Node = filter(subTmetrics$edgeT2, from == subTmetrics$rootT2) %>%
    arrange(to)
  
  # Matrice des produits cartésiens avec sous-noeuds de T2 en colonne et de T1 en ligne
  associations = matrix(nrow = nrow(currentT1Node), ncol = nrow(currentT2Node), 
                        dimnames = list(c(paste0("T1",currentT1Node[,2])),
                                        c(paste0("T2",currentT2Node[,2]))))
  matrice = associations # matrice décompte longueur similarité
  
  for (i in 1:nrow(currentT1Node)){
    for (j in 1:nrow(currentT2Node)){
      associations[i,j] = maststep(subrootT1 = currentT1Node[i,2], 
                              subrootT2 = currentT1Node[j,2], 
                              T1, T2, mastlist= c())
      matrice[i,j] = nchar(associations[i,j])
    }
  }
  
  # Inverser la matrice pour résoudre maximisation avec Algo Hongrois 
  matrice = abs(matrice - max(matrice))
  bestmatches = HungarianSolver(matrice)$pairs
  
  # Concaténation du meilleur groupe d'associations
  best = c()
  for (i in 1:nrow(bestmatches)){
    best = append(liste, associations[bestmatches[i,1], bestmatches[i,2]])
  }
  best = paste0(best, collapse="")
  
  return(best)
}

mast <- function (T1, T2){
  T1T2 = maststep(subrootT1 = T1$node.label[1], subrootT2 = T2$node.label[1],
                  T1 = T1, T2 = T2, mastlist = c())
  
  T2T1 = maststep(subrootT1 = T2$node.label[1], subrootT2 = T1$node.label[1],
                  T1 = T2, T2 = T1, mastlist = c())
  
  unions = mastunion(T1,T2)
  
  
}


# Réservoir
liste = c("abc", "de", "gklm")
res = liste[which.max(lapply(liste, nchar))]

matrix = matrix(ncol =length(T2$node.label) , nrow = length(T1$node.label),
                dimnames = list(c(paste0("T1",T1$node.label),(c(paste0("T2",T2$node.label))))))

test = c(paste0("T1",T1$node.label))
length(T1$node.label)

length(T1$node.label)
list(list(paste0("T1",T1$node.label),((paste0("T2",T2$node.label)))))


