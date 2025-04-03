library(ape)
library(dplyr)
library(stringr)
library(RcppHungarian)
# Fonction pour faire la phylo et la taxo.
# Retourne un vecteur avec taxo, phylo, arbre d'origine
createTaxPhy <- function(nbLeaves = 10, nb.move = 1){
  world = rtree(nbLeaves, rooted = FALSE)
  (q1 = summary(world$edge.length)[2])
  ttt = sample(2:(world$Nnode-1), size = 2*nb.move, replace = F)
  tip1 = ttt[1:nb.move]
  tip2 = ttt[(nb.move+1):(nb.move*2)]
  
  
  # Dichotomie vers polytomie ####
  #' [Dichotomie vers polytomie]
  par(xpd = TRUE) # Permettre à la légende de sortir du cadre de la figure
  layout(matrix(c(1,2),1,2)) # Matrice pour tracer les plots
  world = di2multi(world, tol= round(q1, digits = 1)) # Binaire vers polytomies
  world = makeNodeLabel(world, method = "number", prefix = "n")
  
  # Création arbre taxo et phylo ####
  #' [Visu taxo et phylo]
  layout(matrix(c(1,2,3),1,3)) # Matrice pour tracer les plots
  
  plot(world, cex = 1.5, main = "World", font = 2)
  nodelabels(world$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(world$tip.label[tip1], tip1, adj=0, font = 2, cex = 1.2)
  tiplabels(world$tip.label[tip2], tip2, adj=0, bg = "lightblue", font = 2, cex = 1.2)
  
  taxo = drop.tip(world, tip1) # Enlever une feuille
  tipToChange = match(world$tip.label[tip2], taxo$tip.label)
  taxo = makeNodeLabel(taxo, method = "number", prefix = "tax")
  plot(taxo, main = "Taxo", cex = 1.2, font = 2)
  nodelabels(taxo$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(world$tip.label[tip2], tipToChange, adj=0, bg = "lightblue", font = 2, cex = 1.2)
  
  phylo = drop.tip(world, tip2) # Enlever une feuille
  tipToChange = match(world$tip.label[tip1], phylo$tip.label)
  phylo = makeNodeLabel(phylo, method = "number", prefix = "phy")
  phylo$tip.label[tipToChange] = world$tip.label[tip2] # Changer le nom de la feuille qui est "déplacée"
  plot(phylo, main = "Phylo", cex = 1.2, font = 2)
  nodelabels(phylo$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(world$tip.label[tip2], tipToChange, adj=0, bg = "orchid", font = 2, cex = 1.2)
  

  
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
  nbT1Nodes = T1$Nnode
  
  nbT2Leaves = length(T2$tip.label)
  leavesT2 = T2$tip.label
  rootT2 = T2$node.label[1]
  edgeT2 = edgesToDf(T2)
  nbT2Nodes = T2$Nnode
  
  metrics = list(nbT1Leaves = nbT1Leaves,
                 nbT1Nodes = nbT1Nodes,
                 leavesT1= leavesT1,
                 rootT1 = rootT1,
                 edgeT1 = edgeT1,
                 nbT2Leaves = nbT2Leaves,
                 nbT2Nodes = nbT2Nodes,
                 leavesT2 = leavesT2,
                 rootT2 = rootT2,
                 edgeT2 = edgeT2)
  
  return(metrics)
}

# Création des arbres à tester 
x = createTaxPhy(nbLeaves = 20, nb.move = 4)

T1 = x$taxo
T2 = x$phylo


# Mast en cas d'étages différents
maststep <- function(subrootT1, subrootT2, T1, T2){ #On garde les arbres d'origine
  # Initialisation des variables ####
  
    mastlist = list(c(),c(),c()) # liste de stockage
  
    # Utiliser les subroot pour faire des sous arbres si la racine est différente
    print(paste0("SBRT 1 :", subrootT1," depuis ", T1$Nnode, " noeuds", ";  SBRT 2 : ", subrootT2, " depuis ", T2$Nnode, " noeuds")) 
    
    if (subrootT1 != T1$node.label[1]){
      subT1 = extract.clade(T1, subrootT1)
    }
    else { subT1 = T1
    }
    
    if (subrootT2 != T2$node.label[1]){
      subT2 = extract.clade(T2, subrootT2)
    }
    else { subT2 = T2
    }
    
    # Calculer des métriques sur les sous arbres
    subTmetrics = treeMetrics(T1 = subT1, T2 = subT2)
    
    # Conditions de raccourci
    ## Arrêter s'il n'y a aucune feuille en commun
    commonLeaves = (subTmetrics$leavesT1 %in% subTmetrics$leavesT2)
    overlap = sum(commonLeaves, na.rm=T)
    if (overlap == 0){
      print(paste0("NOMATCH SBRT 1 :", subrootT1," depuis ", T1$Nnode, " noeuds", ";  SBRT 2 : ", subrootT2, " depuis ", T2$Nnode, " noeuds")) 
      return ("")
    }
    ## Si il y a un overlap de 1, retourner la feuille en commun
    else if (overlap == 1){
      return (subTmetrics$leavesT1[which(commonLeaves)])
    }
    ## Si il y a un overlap de 1, retourner la feuille en commun
    else if (overlap == 2 & length(commonLeaves)>= 2){
      return (paste0(subTmetrics$leavesT1[which(commonLeaves)]
                     , collapse ="" ))
    }
    
    # Tableaux des paths du noeud en cours 
    currentT1Node = filter(subTmetrics$edgeT1, from == subTmetrics$rootT1) %>%
                      arrange(to)
    currentT2Node = filter(subTmetrics$edgeT2, from == subTmetrics$rootT2) %>%
                      arrange(to)
  
  # Boucle T1 sur sous-arbre de T2 ####
    # Si premier noeud de l'arbre, séparer les résultats des sous arbres
    print("1 : ")
    if (subTmetrics$rootT2 == "n1"){
      for (firstnode in currentT2Node[,2]){
        if (firstnode %in% subTmetrics$leavesT2) {  # si le sous-noeud est une feuille
          if (firstnode %in% subTmetrics$leavesT1){ # si cette feuille appartient à l'autre sous-arbre
          mastlist[[1]] = append(mastlist[[1]], firstnode)
          }
        }
        else{
          print(paste0("entrée mastep 1; premier noeud ", firstnode, " de ", T2$Nnode, " noeuds dans ", T1$Nnode, " noeuds "))
          
          mastlist[[1]] = append(mastlist[[1]], maststep(subrootT1, firstnode, T1, T2))
          
          print(paste0("sortie mastep 1; premier noeud ", firstnode, " de ", T2$Nnode, " noeuds dans ", T1$Nnode, " noeuds "))
        }
      }
    }
    
    # Si on est dans un des sous arbre, faire le maststep récursivement
    else { 
     for (subnodeT2 in currentT2Node[,2]){ # i = sous-noeuds du noeud en cours
       print(paste0("1 : Sous-noeud en cours ", subnodeT2, " de ", T2$Nnode," noeuds dans ", subrootT1))
      
      if (subnodeT2 %in% subTmetrics$leavesT2) {  # si le sous-noeud est une feuille
        if (subnodeT2 %in% subTmetrics$leavesT1){ # si cette feuille appartient à l'autre sous-arbre
          mastlist[[1]] = append(mastlist[[1]], subnodeT2)
          print(paste0("1 : Mastlist in ", (paste0(mastlist[[1]], collapse = " "))))
          }
      } else {
        # Ajoute la mastlist des enfants a celle du noeud en cours
        print(paste0("ENTREE mastep 1; sous noeud ", subnodeT2, " de ", T2$Nnode, " noeuds dans ", subrootT1))
        
        mastlist[[1]] = append(mastlist[[1]], maststep(subrootT1, subnodeT2, T1, T2))
        
        print(paste0("SORTIE mastep 1; sous noeud ", subnodeT2, " de ", T2$Nnode, " noeuds dans ", subrootT1))
      }
     }
      # Concaténer les feuilles
      mastlist[[1]] = paste0(mastlist[[1]], collapse = "")
    }
  
  print(paste0("1 : Mastlist end ", (paste0(mastlist[[1]], collapse = " "))))
  
  # Boucle T2 sur sous-arbre de T1 ####
  print("2 : ")
    # Si premier noeud de l'arbre, séparer les résultats des sous arbres
    if (subTmetrics$rootT1 == "n1"){
      for (firstnode in currentT1Node[,2]){
        if (firstnode %in% subTmetrics$leavesT1) {  # si le sous-noeud est une feuille
          if (firstnode %in% subTmetrics$leavesT2){ # si cette feuille appartient à l'autre sous-arbre
          mastlist[[2]] = append(mastlist[[2]], firstnode)
          }
        }
        else {
          print(paste0("ENTREE mastep 2; premier noeud ", firstnode, " de ", T1$Nnode, " noeuds dans ", T2$Nnode, " noeuds "))
          
          mastlist[[2]] = append(mastlist[[2]], maststep(subrootT2, firstnode, T2, T1))
          print(paste0("SORTIE mastep 2; premier noeud ", firstnode, " de ", T1$Nnode, " noeuds dans ", T2$Nnode, " noeuds "))
        }
      }
    }
    
    # Si on est dans un des sous arbre, faire le maststep récursivement
    else { 
      for (subnodeT1 in currentT1Node[,2]){
        print(paste0("2 : Sous-noeud en cours ", subnodeT1, " de ", T1$Nnode," noeuds dans ", subrootT2))
        
        if (subnodeT1 %in% subTmetrics$leavesT1) {  # si le sous-noeud est une feuille
          if (subnodeT1 %in% subTmetrics$leavesT2){ # si cette feuille appartient à l'autre sous-arbre
            mastlist[[2]] = append(mastlist[[2]], subnodeT1)
            print(paste0("2 : Mastlist in ", (paste0(mastlist[[2]], collapse = " "))))
          }
        } else {
          # Ajoute la mastlist des enfants a celle du noeud en cours
          print(paste0("ENTREE mastep 2; sous noeud ", subnodeT1, " de ", T1$Nnode, " noeuds dans ", subrootT2))
          mastlist[[2]] = append(mastlist[[2]], maststep(subrootT2, subnodeT1, T2, T1))
          print(paste0("SORTIE mastep 2; sous noeud ", subnodeT1, " de ", T1$Nnode, " noeuds dans ", subrootT2))
        }
      }
      # Concaténer les feuilles
      mastlist[[2]] = paste0(mastlist[[2]], collapse = "")
    }
  print(paste0("2 : Mastlist end ", (paste0(mastlist[[2]], collapse = " "))))
  
  # Matching des sous-arbres ####
  print("3 : ")
    # Matrice des produits cartésiens avec sous-noeuds de T2 en colonne et de T1 en ligne
    associations = matrix(nrow = nrow(currentT1Node), ncol = nrow(currentT2Node), 
                          dimnames = list(c(paste0("T1",currentT1Node[,2])),
                                          c(paste0("T2",currentT2Node[,2]))))
    countMat = associations # matrice compte longueur similarité
    
    for (i in 1:nrow(currentT1Node)){
      for (j in 1:nrow(currentT2Node)){
        iNode = currentT1Node[i,2]
        jNode = currentT2Node[j,2]
        print (paste0("i : ", iNode, "; j : ", jNode))
        
        # Si l'un des noeuds courant est une feuille
        if ((iNode %in% subTmetrics$leavesT1) && 
            (jNode %in% subTmetrics$leavesT2)){
          if (iNode==jNode){ # sont elles identiques ?
            associations[i,j] = iNode 
            countMat[i,j] = 1
          } else {
            associations[i,j] = NA 
            countMat[i,j] = 0 
          }
        } else if ((iNode %in% subTmetrics$leavesT1) || 
                 (jNode %in% subTmetrics$leavesT2)){
          # Sinon la feuille est-elle dans l'autre sous-arbre ?
          if ((iNode %in% subTmetrics$leavesT1) &&
                     (iNode %in% extract.clade(T2, jNode)$tip.label)){
            associations[i,j] = iNode 
            countMat[i,j] = 1
            
          } else if ((jNode %in% subTmetrics$leavesT2) &&
                     (jNode %in% extract.clade(T1, iNode)$tip.label)){
            associations[i,j] = jNode 
            countMat[i,j] = 1
          
          } else {
            associations[i,j] = NA 
            countMat[i,j] = 0
          }
        } else {
          print(paste0("SORTIE mastep 3; sous noeuds ", iNode, " et ", jNode))
          associations[i,j] = maststep(subrootT1 = iNode, 
                                       subrootT2 = jNode, 
                                       T1, T2)
          print(paste0("SORTIE mastep 3; sous noeuds ", iNode, " et ", jNode))
          countMat[i,j] = str_count(associations[i,j], pattern = "t")
        }
        print(countMat[i,j])
      }  
    }
    
    # Inverser la countMat pour résoudre maximisation avec Algo Hongrois 
    countMat = abs(countMat - max(countMat))
    print("countMatrix :")
    print(countMat)
    bestmatches = HungarianSolver(countMat)$pairs
    
    # Retire les matchs sans optimum dans une matrice rectangle
    if (!(identical(which(bestmatches[,2]==0), integer(0)))) # Vérifie que le which n'es pas vide
      bestmatches = bestmatches[-which(bestmatches[,2]==0),] 
    
    # Concaténation du meilleur groupe d'associations
    for (i in 1:nrow(bestmatches)){
      mastlist[[3]] = append(mastlist[[3]], associations[bestmatches[i,1], bestmatches[i,2]])
    }
    mastlist[[3]] = paste0(mastlist[[3]], collapse="")
    print(paste0("3 : Mastlist end ", (paste0(mastlist[[3]], collapse = " "))))
    
  # Choisir le max de la mastlist ####
    besthit = c()
    for (i in 1:length(mastlist)){
      print(paste0("Mastlist : ",i, " ", (paste0(mastlist[[i]], collapse = " "))))
      
      besthit = append(besthit, 
                       mastlist[[i]][which.max(lapply(mastlist[[i]], str_count, pattern = "t"))])
    }
    
    print(paste0("besthits :",besthit))
    print(besthit[which.max(lapply(besthit, str_count, pattern = "t"))])
    return (besthit[which.max(lapply(besthit, str_count, pattern = "t"))])
}

resultat = maststep(subrootT1 = T1$node.label[1], subrootT2 = T2$node.label[1],
     T1 = T1, T2 = T2)
y = keep.tip(T2, c("t10","t7","t3"))


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
  countMat = associations # matrice compte longueur similarité
  
  for (i in 1:nrow(currentT1Node)){
    for (j in 1:nrow(currentT2Node)){
      associations[i,j] = maststep(subrootT1 = currentT1Node[i,2], 
                              subrootT2 = currentT1Node[j,2], 
                              T1, T2, mastlist)
      countMat[i,j] = nchar(associations[i,j])
    }
  }  
  
  # Inverser la countMat pour résoudre maximisation avec Algo Hongrois 
  countMat = abs(countMat - max(countMat))
  bestmatches = HungarianSolver(matrice)$pairs
  
  # Retire les matchs sans optimum dans une matrice rectangle
  bestmatches = bestmatches[-which(bestmatches[,2]==0),] 
  
  # Concaténation du meilleur groupe d'associations
  for (i in 1:nrow(bestmatches)){
    mastlist[[3]] = append(mastlist[[3]], associations[bestmatches[i,1], bestmatches[i,2]])
  }
  mastlist[[3]] = paste0(mastlist[[3]], collapse="")
}
mastunion(T1,T2)

# test zone
test <- function (a,b){
  if (a==b) {
    return (print("aouioui"))
  }
  return (a+b)
}

# Réservoir
liste = c("abc", "de", "gklm")
res = liste[which.max(lapply(liste, nchar))]

matrix = matrix(ncol =T2$Nnode , nrow = T1$Nnode,
                dimnames = list(c(paste0("T1",T1$node.label),(c(paste0("T2",T2$node.label))))))

test = c(paste0("T1",T1$node.label))
T1$Nnode
comparePhylo()
T1$Nnode
list(list(paste0("T1",T1$node.label),((paste0("T2",T2$node.label)))))

