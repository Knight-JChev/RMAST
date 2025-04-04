library(ape)
library(dplyr)
library(stringr)
library(RcppHungarian)
# Fonction pour faire la phylo et la taxo.
# Retourne un vecteur avec taxo, phylo, arbre d'origine
createTaxPhy <- function(nbLeaves = 5, nb.move = 1){
  world = rtree(nbLeaves, rooted = FALSE)
  (q1 = summary(world$edge.length)[2])
  ttt = sample(1:length(world$tip.label), size = 2*nb.move, replace = F)
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
  tiplabels(world$tip.label[tip1], tip1, adj=0, font = 2, cex = 1.2, bg = "mediumpurple1")
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
  
  
  taxo[[6]] = "taxo"
  names(taxo)[6] = "name"
  
  phylo[[6]] = "phylo"
  names(phylo)[6] = "name"
  
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
treeMetrics <- function (taxo, phylo){
  nbTaxoLeaves = length(taxo$tip.label)
  leavesTaxo = taxo$tip.label # liste des noms des feuilles
  rootTaxo = taxo$node.label[1]
  edgeTaxo = edgesToDf(taxo)
  nbTaxoNodes = taxo$Nnode
  
  nbPhyloLeaves = length(phylo$tip.label)
  leavesPhylo = phylo$tip.label
  rootPhylo = phylo$node.label[1]
  edgePhylo = edgesToDf(phylo)
  nbPhyloNodes = phylo$Nnode
  
  metrics = list(nbTaxoLeaves = nbTaxoLeaves,
                 nbTaxoNodes = nbTaxoNodes,
                 leavesTaxo= leavesTaxo,
                 rootTaxo = rootTaxo,
                 edgeTaxo = edgeTaxo,
                 nbPhyloLeaves = nbPhyloLeaves,
                 nbPhyloNodes = nbPhyloNodes,
                 leavesPhylo = leavesPhylo,
                 rootPhylo = rootPhylo,
                 edgePhylo = edgePhylo)
  
  return(metrics)
}

# Création des arbres à tester 
x = createTaxPhy(nbLeaves = 50, nb.move = 5)

Taxo = x$taxo
Phylo = x$phylo
zeMatrix = matrix(nrow = Taxo$Nnode, 
                  ncol = Phylo$Nnode, 
                  dimnames = list(Taxo$node.label, Phylo$node.label))

# Mast en cas d'étages différents
maststep <- function(subRootTax, subRootPhy, trees){ #On garde les arbres d'origine
  # Initialisation des variables ####
  
    mastlist = list(c(),c(),c()) # liste de stockage
    
    whichTaxo = which(c(trees[[1]]$name, trees[[2]]$name)=="taxo")
    taxo = trees[[whichTaxo]]
    phylo = trees[[base::setdiff(1:2,whichTaxo)]]
    
    # Utiliser les subroot pour faire des sous arbres si la racine est différente
    print(paste0("SBRT 1 :", subRootTax," depuis Taxo", ";  SBRT 2 : ", subRootPhy, " depuis Phylo")) 
    
    if (subRootTax != taxo$node.label[1]){
      subTax = extract.clade(taxo, subRootTax)
    }
    else { subTax = taxo
    }
    
    if (subRootPhy != phylo$node.label[1]){
      subPhy = extract.clade(phylo, subRootPhy)
    }
    else { subPhy = phylo
    }
    
    # Calculer des métriques sur les sous arbres
    subTmetrics = treeMetrics(taxo = subTax, phylo = subPhy)
    
    # Conditions de raccourci
    ## Arrêter s'il n'y a aucune feuille en commun
    commonLeaves = (subTmetrics$leavesTaxo %in% subTmetrics$leavesPhylo)
    overlap = sum(commonLeaves, na.rm=T)
    if (overlap == 0){
      print(paste0("NOMATCH SBRT 1 :", subRootTax," depuis Taxo", ";  SBRT 2 : ", subRootPhy, " depuis Phylo")) 
      return ("")
    }
    ## Si il y a un overlap de 1, retourner la feuille en commun
    else if (overlap == 1){
      return (subTmetrics$leavesTaxo[which(commonLeaves)])
    }
    ## Si il y a un overlap de 1, retourner la feuille en commun
    else if (overlap == 2 & length(commonLeaves)>= 2){
      return (paste0(subTmetrics$leavesTaxo[which(commonLeaves)], collapse ="" ))
    }
    
    # Tableaux des paths du noeud en cours 
    currentTaxoNode = filter(subTmetrics$edgeTaxo, from == subTmetrics$rootTaxo) %>%
                      arrange(to)
    currentPhyloNode = filter(subTmetrics$edgePhylo, from == subTmetrics$rootPhylo) %>%
                      arrange(to)
  
  # Boucle Taxo sur sous-arbre de Phylo ####
    print("1 :")
    # Si on est dans un des sous arbre, faire le maststep récursivement
     for (subnodePhylo in currentPhyloNode[,2]){ # i = sous-noeuds du noeud en cours
       print(paste0("1 : Sous-noeud en cours ", subnodePhylo, " de Phylo dans ", subRootTax, " de Taxo"))
      
      if (subnodePhylo %in% subTmetrics$leavesPhylo) {  # si le sous-noeud est une feuille
        if (subnodePhylo %in% subTmetrics$leavesTaxo){ # si cette feuille appartient à l'autre sous-arbre
          mastlist[[1]] = append(mastlist[[1]], subnodePhylo)
          print(paste0("1 : Mastlist in ", (paste0(mastlist[[1]], collapse = " "))))
          }
      } else {
        # Ajoute la mastlist des enfants a celle du noeud en cours
        if (is.na(zeMatrix[subRootTax, subnodePhylo])){
          print(paste0("ENTREE mastep 1; sous noeud ", subnodePhylo, " de Phylo dans ", subRootTax, " de Taxo"))
          .GlobalEnv$zeMatrix[subRootTax, subnodePhylo] = maststep(subRootTax, subnodePhylo, list(taxo,phylo))
          print(paste0("SORTIE mastep 1; sous noeud ", subnodePhylo, " de Phylo dans ", subRootTax, " de Taxo"))
        }
        mastlist[[1]] = append(mastlist[[1]], zeMatrix[subRootTax, subnodePhylo])
        print(paste0("1 : Mastlist in ", (paste0(mastlist[[1]], collapse = " "))))
      }
     }
  
  print(paste0("1 : Mastlist end ", (paste0(mastlist[[1]], collapse = " "))))
  
  # Boucle Phylo sur sous-arbre de Taxo ####
    print("2 : ")
    # Si on est dans un des sous arbre, faire le maststep récursivement
      for (subnodeTaxo in currentTaxoNode[,2]){
        print(paste0("2 : Sous-noeud en cours ", subnodeTaxo, " de taxo dans ", subRootPhy))
        
        if (subnodeTaxo %in% subTmetrics$leavesTaxo) {  # si le sous-noeud est une feuille
          if (subnodeTaxo %in% subTmetrics$leavesPhylo){ # si cette feuille appartient à l'autre sous-arbre
            mastlist[[2]] = append(mastlist[[2]], subnodeTaxo)
            print(paste0("2 : Mastlist in ", (paste0(mastlist[[2]], collapse = " "))))
          }
        } else {
          # Ajoute la mastlist des enfants a celle du noeud en cours
          if (is.na(zeMatrix[subnodeTaxo, subRootPhy])){
            print(paste0("ENTREE mastep 2; sous noeud ", subnodeTaxo, " de Taxo dans ", subRootPhy, " de Phylo"))
            .GlobalEnv$zeMatrix[subnodeTaxo, subRootPhy] = maststep(subnodeTaxo, subRootPhy, list(taxo,phylo))
            print(paste0("SORTIE mastep 2; sous noeud ", subnodeTaxo, " de Taxo dans ", subRootPhy, " de Phylo"))
          }
          mastlist[[2]] = append(mastlist[[2]], zeMatrix[subnodeTaxo, subRootPhy])
        }
      }
  print(paste0("2 : Mastlist end ", (paste0(mastlist[[2]], collapse = " "))))
  
  # Matching des sous-arbres ####
  print("3 : ")
    # Matrice des produits cartésiens avec sous-noeuds de phylo en colonne et de taxo en ligne
    associations = matrix(nrow = nrow(currentTaxoNode), ncol = nrow(currentPhyloNode), 
                          dimnames = list(c(paste0("taxo",currentTaxoNode[,2])),
                                          c(paste0("phylo",currentPhyloNode[,2]))))
    countMat = associations # matrice compte longueur similarité
    
    for (i in 1:nrow(currentTaxoNode)){
      for (j in 1:nrow(currentPhyloNode)){
        iNode = currentTaxoNode[i,2]
        jNode = currentPhyloNode[j,2]
        print (paste0("i : ", iNode, "; j : ", jNode))
        
        # Si l'un des noeuds courant est une feuille
        if ((iNode %in% subTmetrics$leavesTaxo) && 
            (jNode %in% subTmetrics$leavesPhylo)){
          if (iNode==jNode){ # sont elles identiques ?
            associations[i,j] = iNode 
            countMat[i,j] = 1
          } else {
            associations[i,j] = NA 
            countMat[i,j] = 0 
          }
        } else if ((iNode %in% subTmetrics$leavesTaxo) || 
                 (jNode %in% subTmetrics$leavesPhylo)){
          # Sinon la feuille est-elle dans l'autre sous-arbre ?
          if ((iNode %in% subTmetrics$leavesTaxo) &&
                     (iNode %in% extract.clade(phylo, jNode)$tip.label)){
            associations[i,j] = iNode 
            countMat[i,j] = 1
            
          } else if ((jNode %in% subTmetrics$leavesPhylo) &&
                     (jNode %in% extract.clade(taxo, iNode)$tip.label)){
            associations[i,j] = jNode 
            countMat[i,j] = 1
          
          } else {
            associations[i,j] = NA 
            countMat[i,j] = 0
          }
        } else {
          print(paste0("SORTIE mastep 3; sous noeuds ", iNode, " et ", jNode))
          associations[i,j] = maststep(subRootTax = iNode, 
                                       subRootPhy = jNode, 
                                       list(taxo, phylo))
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

zeMatrix = matrix(nrow = Taxo$Nnode, 
                  ncol = Phylo$Nnode, 
                  dimnames = list(Taxo$node.label, Phylo$node.label))

resultat = maststep(subRootTax = Taxo$node.label[1], subRootPhy = Phylo$node.label[1],
     trees = list(Taxo,Phylo))

y = keep.tip(T2, c("t10","t7","t3"))

t30t6t36t44t24t15t31t4t40t20t19t32t29t43t34t25t33t35t47t12t49t42t17t14t7t21t23t38t27t41t48t2t16t50t1t11t26t3t13t46
# test zone ####
test <- function (a,b){
  if (a==b) {
    return (print("aouioui"))
  }
  return (a+b)
}

# Réservoir ####
liste = c("abc", "de", "gklm")
res = liste[which.max(lapply(liste, nchar))]

base::setdiff(1:2)
matrix = matrix(ncol =T2$Nnode , nrow = T1$Nnode,
                dimnames = list(c(paste0("T1",T1$node.label),(c(paste0("T2",T2$node.label))))))

test = c(paste0("T1",T1$node.label))
T1$Nnode
comparePhylo()
T1$Nnode
list(list(paste0("T1",T1$node.label),((paste0("T2",T2$node.label)))))

