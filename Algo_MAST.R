library(ape)
library(dplyr)
library(stringr)
library(RcppHungarian)
# Fonction pour faire la phylo et la taxo.
# Retourne un vecteur avec taxo, phylo, arbre d'origine
createTaxPhy <- function(nbLeaves = 10){
  world = rtree(nbLeaves+nb.move, rooted = FALSE)
  (q1 = summary(world$edge.length)[2])
  ttt = sample(1:length(world$tip.label), size = 2*nb.move, replace = F)
  tip1 = ttt[1:nb.move]
  tip2 = ttt[(nb.move+1):(nb.move*2)]
  
  tipsNumber = sub(".", "", world$tip.label)
  wrongTips = world$tip.label[tip2]
  
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
  tiplabels(wrongTips, tip2, adj=0, bg = "lightblue", font = 2, cex = 1.2)
  
  taxo = drop.tip(world, tip1) # Enlever une feuille
  tipToChange = match(wrongTips, taxo$tip.label)
  taxo = makeNodeLabel(taxo, method = "number", prefix = "tax")
  plot(taxo, main = "Taxo", cex = 1.2, font = 2)
  nodelabels(taxo$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(wrongTips, tipToChange, adj=0, bg = "lightblue", font = 2, cex = 1.2)
  
  phylo = drop.tip(world, tip2) # Enlever une feuille
  tipToChange = match(world$tip.label[tip1], phylo$tip.label)
  phylo = makeNodeLabel(phylo, method = "number", prefix = "phy")
  phylo$tip.label[tipToChange] = wrongTips # Changer le nom de la feuille qui est "déplacée"
  plot(phylo, main = "Phylo", cex = 1.2, font = 2)
  nodelabels(phylo$node.label, adj = c(1.3,-0.5), frame = "n", cex = 1.5, font = 2, col="red")
  tiplabels(wrongTips, tipToChange, adj=0, bg = "orchid", font = 2, cex = 1.2)
  
  # Nombre d'arrêtes entre les paires de feuilles déplacées ####
  edgeDist = c()
  notWrongTips = c()
  wrongTips = sub(".","",wrongTips)
  
  for (i in 1:nb.move) {
    edges = (length(nodepath(world, tip2[i], tip1[i])))-1
    if (edges == 2) {
      notWrongTips = append(notWrongTips, tipsNumber[tip2[i]])
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
  
  return(list(taxo = taxo, phylo = phylo, world = world, 
              notWrongTips = notWrongTips, truth = c(tipsNumber[-c(tip2,tip1)], notWrongTips),
              edgeDist = edgeDist, realWrongTips = realWrongTips, wrongTips = wrongTips))
}

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
  plot(taxo, cex = 1, main = "Taxo", font = 2)
  nodelabels(taxo$node.label, adj = c(1,-0.2), frame = "n", cex = 0.8, font = 2, col="red")
  tiplabels(taxo$tip.label[tipTo], tipTo, adj=0, font = 2, cex = 1, bg = "mediumpurple1")
  tiplabels(wrongTips, tipMoved, adj=0, bg = "lightblue", font = 2, cex = 1)
  
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
  plot(phylo, main = "Phylo", cex = 1, font = 2)
  nodelabels(phylo$node.label, adj = c(1,-0.2), frame = "n", cex = 0.8, font = 2, col="red")
  tiplabels(taxo$tip.label[tipTo], tipToAdjust, adj=0, bg = "mediumpurple1", font = 2, cex = 1)
  tiplabels(taxo$tip.label[tipMoved], tipMovedTo, adj=0, bg = "lightblue", font = 2, cex = 1)
  
  
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

# Fonction pour faire une dataframe avec le nom des feuilles et noeuds
edgesToDf <- function(tree){
  nbNodesLeaves = (length(tree$tip.label)+Nnode(tree))
  
  # Concaténer noms des feuilles et noeuds
  newLabels = c(tree$tip.label, tree$node.label)
  
  # Remplacer nombres du path par les noms des noeuds et feuilles
  edgePaths = as.data.frame(tree$edge)
  edgePaths = edgePaths %>%
    rename("from"="V1","to"="V2") %>%
    mutate(from = recode(from, !!!setNames(newLabels, as.character(1:nbNodesLeaves)))) %>%
    mutate(to = recode(to, !!!setNames(newLabels, as.character(1:nbNodesLeaves))))
  
  return(edgePaths)
}

# Métriques et repères sur l'arbre
treeMetrics <- function (taxo, phylo){

    metrics = list(nbTaxoLeaves = length(taxo$tip.label),
                 nbTaxoNodes = taxo$Nnode,
                 leavesTaxo= taxo$tip.label,
                 rootTaxo = taxo$node.label[1],
                 edgeTaxo = edgesToDf(taxo),
                 nbPhyloLeaves = length(phylo$tip.label),
                 nbPhyloNodes = phylo$Nnode,
                 leavesPhylo = phylo$tip.label,
                 rootPhylo = phylo$node.label[1],
                 edgePhylo = edgesToDf(phylo))
  
  return(metrics)
}

# Trouver le meilleur hit dans une liste
findBest <- function(hitlist) {
  # Find index of longest "vector" in the hitlist
  hits = unlist(lapply(hitlist, strsplit, split = ","), recursive = F)
  hitlength = unlist(lapply(hits, length))
  maxhits = which(hitlength == max(hitlength))
  
  # Pick max or random in case of equalities
  if (length(maxhits) > 1){
    besthits = hitlist[[sample(maxhits,1)]]
  } else {
    besthits = hitlist[[maxhits]]
  }
  return(besthits)
}

# Maximum agreement subtree
mast <- function(taxo, phylo, subRootTax, subRootPhy){ #On garde les arbres d'origine
  # Initialisation des variables ####
  
  mastlist = list(c(),c(),c()) # liste de stockage
  
  # Utiliser les subroot pour faire des sous arbres si la racine est différente
  print(paste0("SBRT 1 :", subRootTax," depuis Taxo", ";  SBRT 2 : ", subRootPhy, " depuis Phylo")) 
  print(taxo$node.label[1])
  
  if (subRootTax != taxo$node.label[1]){
    subTax = extract.clade(taxo, subRootTax)
  } else { subTax = taxo
  }
  
  if (subRootPhy != phylo$node.label[1]){
    subPhy = extract.clade(phylo, subRootPhy)
  } else { subPhy = phylo
  }
  
  # Calculer des métriques sur les sous arbres
  subTmetrics = treeMetrics(taxo = subTax, phylo = subPhy)
  
  # Conditions de raccourci ####
  ## Arrêter s'il n'y a aucune feuille en commun
  commonLeaves = (subTmetrics$leavesTaxo %in% subTmetrics$leavesPhylo)
  print(subTmetrics$leavesTaxo[which(commonLeaves)])
  overlap = sum(commonLeaves, na.rm=T)
  if (overlap == 0){
    print(paste0("NOMATCH SBRT 1 :", subRootTax," depuis Taxo", ";  SBRT 2 : ", subRootPhy, " depuis Phylo")) 
    return ("")
  
  ## Si il y a une feuille en commun; retourner la feuille
  }else if (overlap == 1){
    print(paste0("1 MATCH ", subTmetrics$leavesTaxo[which(commonLeaves)]))
    return (subTmetrics$leavesTaxo[which(commonLeaves)])
  
  ## Si il y a deux feuilles en commun; retourner les feuille
  }else if (overlap == 2 & length(commonLeaves)>= 2){
    print(paste0("2 MATCH ", paste0(subTmetrics$leavesTaxo[which(commonLeaves)], collapse = "," )))
    return (paste0(subTmetrics$leavesTaxo[which(commonLeaves)], collapse = "," ))
  }
  
  # Tableaux des paths du noeud en cours 
  currentTaxoNode = filter(subTmetrics$edgeTaxo, from == subTmetrics$rootTaxo) %>%
    arrange(to)
  currentPhyloNode = filter(subTmetrics$edgePhylo, from == subTmetrics$rootPhylo) %>%
    arrange(to)
  
  # Boucle Taxo sur sous-arbre de Phylo ####
  print("1 :")
  # Si on est dans un des sous arbre, faire le mast récursivement
  for (subnodePhylo in currentPhyloNode[,2]){ # i = sous-noeuds du noeud en cours
    if (subnodePhylo %in% subTmetrics$leavesPhylo) {  # si le sous-noeud est une feuille
      if (subnodePhylo %in% subTmetrics$leavesTaxo){ # si cette feuille appartient à l'autre sous-arbre
        mastlist[[1]] = append(mastlist[[1]], subnodePhylo)
      }
    } else {
      # Ajoute la mastlist des enfants a celle du noeud en cours
      if (is.na(doneMat[subRootTax, subnodePhylo])){
        .GlobalEnv$doneMat[subRootTax, subnodePhylo] = mast(taxo, phylo, subRootTax, subnodePhylo)
      }
      mastlist[[1]] = append(mastlist[[1]], doneMat[subRootTax, subnodePhylo])
    }
  }
  
  # Boucle Phylo sur sous-arbre de Taxo ####
  print("2 : ")
  # Si on est dans un des sous arbre, faire le mast récursivement
  for (subnodeTaxo in currentTaxoNode[,2]){
    if (subnodeTaxo %in% subTmetrics$leavesTaxo) {  # si le sous-noeud est une feuille
      if (subnodeTaxo %in% subTmetrics$leavesPhylo){ # si cette feuille appartient à l'autre sous-arbre
        mastlist[[2]] = append(mastlist[[2]], subnodeTaxo)
      }
    } else {
      # Ajoute la mastlist des enfants a celle du noeud en cours
      if (is.na(doneMat[subnodeTaxo, subRootPhy])){
        .GlobalEnv$doneMat[subnodeTaxo, subRootPhy] = mast(taxo, phylo, subnodeTaxo, subRootPhy)
      }
      mastlist[[2]] = append(mastlist[[2]], doneMat[subnodeTaxo, subRootPhy])
    }
  }

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

      # Si les deux noeuds courant sont des feuilles
      if ((iNode %in% subTmetrics$leavesTaxo) && 
          (jNode %in% subTmetrics$leavesPhylo)){
        if (iNode==jNode){ # sont elles identiques ?
          associations[i,j] = iNode 
          countMat[i,j] = 1
        } else {
          associations[i,j] = NA 
          countMat[i,j] = 0 
        }
      # Sinon la feuille est-elle dans l'autre sous-arbre ?
      } else if ((iNode %in% subTmetrics$leavesTaxo) || 
                 (jNode %in% subTmetrics$leavesPhylo)){
        
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
        associations[i,j] = mast(taxo, phylo,
                                 subRootTax = iNode, 
                                 subRootPhy = jNode)
        countMat[i,j] = length(strsplit(associations[i,j], split = ","))
      }
    }  
  }
  
  # Inverser la countMat pour résoudre maximisation avec Algo Hongrois 
  countMat = abs(countMat - max(countMat))
  #print("countMatrix :")
  #print(countMat)
  bestmatches = HungarianSolver(countMat)$pairs
  
  # Retire les matchs sans optimum dans une matrice rectangle
  if (!(identical(which(bestmatches[,2]==0), integer(0)))) # Vérifie que le which n'es pas vide
    bestmatches = bestmatches[-which(bestmatches[,2]==0),] 
  
  # Concaténation du meilleur groupe d'associations
  for (i in 1:nrow(bestmatches)){
    mastlist[[3]] = append(mastlist[[3]], associations[bestmatches[i,1], bestmatches[i,2]])
  }
  mastlist[[3]] = gsub("NA,*|,*NA", "", paste0(mastlist[[3]], collapse = ",")) 

  # Choisir le max de la mastlist, aléatoire si plusieurs max ####
  besthits = rep(1,length(mastlist))
  
  for (i in 1:length(mastlist)){
    besthits[i] = findBest(mastlist[[i]])
  }
  res = findBest(besthits)
  
  print(paste0("--FIN-- SBRT 1 :", subRootTax," depuis Taxo", ";  SBRT 2 : ", subRootPhy, " depuis Phylo"))
  print(paste0("besthits :", besthits))
  print(besthits)
  print(paste0("meilleur = ", res))
  return (res) # Meilleur aléatoire
}

# Benchmark sur des arbres aléatoires
benchmark <- function(nbRepeats, nbLeaves, nb.move, phyType = "Normal", dist){
  # Dossiers et fichiers de stockage
  time = format(Sys.time(), "%H%M%S")
  if (phyType == "Alt"){
    dirNameFig = paste0("Benchmark_data/", dist,"dist_",nb.move,"moves_",nbRepeats,"rep_Fig_",time)
    fileName = paste0(dist,"dist_",nb.move,"moves_",nbLeaves,"leaves_",nbRepeats,"rep_",time)
  } else {  
    dirNameFig = paste0("Benchmark_data/",  nb.move,"moves_", nbRepeats,"rep_Fig_",time)
    fileName = paste0(nb.move,"moves_",nbLeaves,"leaves_",nbRepeats,"rep_",time)
  }
  
  ifelse(!dir.exists(file.path(dirNameFig)), dir.create(file.path(dirNameFig)), "Figure Directory already exists")

  # Variable de stockage des résultats
  metrics = data.frame() 
  
  
  # Boucle pour atteindre le nombre de répétitions
  print(paste0("Tree 1 / ",nbRepeats, " ongoing"))
  for (rep in 1:nbRepeats){
    
    # Création + sauvegarde des arbres aléatoires avec feuilles déplacées
    png(filename = paste0(dirNameFig,"/",fileName,"_",rep), width = 1920, height = 1080)
    if (phyType == "Alt"){
      random = createTaxPhyAlt(nbLeaves = nbLeaves, nb.move = nb.move, dist = dist)
    } else  random = createTaxPhy(nbLeaves = nbLeaves, nb.move = nb.move)
    dev.off()
    
    Taxo = random$taxo
    Phylo = random$phylo

    # Matrice d'enregistrement des résultats de mast en cours
    .GlobalEnv$doneMat = matrix(nrow = Taxo$Nnode, 
                      ncol = Phylo$Nnode, 
                      dimnames = list(Taxo$node.label, Phylo$node.label))
    
    # Résultat du mast sous forme de vecteur avec numéro des feuilles
    resultat = strsplit(mast(taxo = Taxo, phylo = Phylo,
                             subRootTax = Taxo$node.label[1], 
                             subRootPhy = Phylo$node.label[1]),"t")[[1]][-1]
    
    # Calcul des métriques : positif = feuille déplacée non retenue
    TPtips = setdiff(random$realWrongTips, resultat) # Feuilles dégagées censées l'être
    FPtips = setdiff(random$truth, resultat) # Feuilles dégagées censées être retenues
    TNtips = resultat[which(resultat %in% random$truth)] # Feuilles retenues censées l'être
    FNtips = resultat[which(resultat %in% random$realWrongTips)] # Feuilles retenues censées être dégagées
    
    TP = length(TPtips)
    FP = length(FPtips)
    TN = length(TNtips)
    FN = length(FNtips)

    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    accuracy = (TP + TN) / length(c(random$truth, random$realWrongTips))
    
    which((setdiff(random$realWrongTips, resultat) %in% random$realWrongTips))

    # Remplissage d'une table pour l'arbre en cours
    df = data.frame(id = rep,
                    nbRepeats = nbRepeats,
                    nbLeaves = nbLeaves,
                    nb.move = as.factor(nb.move),
                    obsTreeSize = length(resultat),
                    expErrorSize = length(random$realWrongTips),
                    TP = TP,
                    FP = FP,
                    TN = TN,
                    FN = FN,
                    Accuracy = accuracy,
                    Precision = precision,
                    Recall = recall)
    df$resultat = list(resultat)
    df$truth = list(random$truth)
    df$wrongTips = list(random$wrongTips)
    df$realWrongTips = list(random$realWrongTips)
    df$TPtips = list(TPtips)
    df$TNtips = list(TNtips)
    df$FPtips = list(FPtips)
    df$FNtips = list(FNtips) 
    df$edgeDist = list(random$edgeDist)
    df$perfect = as.factor(ifelse(df$Accuracy==1, 1,0))
    
    if (phyType == "Alt"){
      df$dist = as.factor(dist)
      df = df %>% relocate(dist, .after = nb.move)
    }
    
    metrics = rbind(metrics, df)
    print(df[1, 1:12])
    
    print(paste0("Tree  ",rep," / ",nbRepeats, " done"))
  }
  
  # Sauvegarde des données
  saveRDS(metrics, file = paste0("Benchmark_data/",fileName,".rds"))
  
  return (metrics)
}

# Parallélisation du benchmark
bench_par <- function(arglist){
  
  Ncpus <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(cl)
  
  res <- foreach::foreach(i=1:length(arglist), 
                   .export =c("benchmark", "createTaxPhyAlt", "mast", 
                              "treeMetrics", "edgesToDf"),
                   .packages = c("ape","RcppHungarian", 
                                 "stringr", "dplyr")) %dopar% {
                                   return(do.call(what = benchmark, arglist[[i]]))
                                 }
  
  parallel::stopCluster(cl)
  return(res)
}

# Matrice d'enregistrement des résultats de mast en cours
arbreTax = read.tree(text = "((7990KR701906,7992KJ473717,11894MT795184,18920KT375565,22559OR546136,26308MK978155,4556MH248251,3295LC549804)Nemipterus,4439AY484975)Eupercaria;0")
arbrePhy = read.tree(text = "((26308MK978155:0.06042889,((11894MT795184:0.00958929,18920KT375565:0.00922456)0.981669:0.01934093,(7990KR701906:0.00000001,7992KJ473717:0.00000001)-1.000000:0.02026485)0.334856:0.00186894)0.334044:0.00724549,4556MH248251:0.00993826,(22559OR546136:0.08719328,(4439AY484975:0.00578967,3295LC549804:0.05351206)0.998239:0.04709670)0.662670:0.02002656);")

arbrePhy = makeNodeLabel(arbrePhy, method ="number")

.GlobalEnv$doneMat = matrix(nrow = arbreTax$Nnode, 
                              ncol = arbrePhy$Nnode, 
                              dimnames = list(arbreTax$node.label, arbrePhy$node.label))

mast(taxo = arbreTax, phylo = arbrePhy, subRootTax = arbreTax$node.label[1], subRootPhy = arbrePhy$node.label[1])

# Setup parallélistation ----
  ##toRun = list(list(200, 100, 10, "Alt", 4),
  ##             list(200, 100, 15, "Alt", 4))
  ##x = bench_par(toRun)

# Lancement manuel d'un benchmark ####
  ##x = benchmark(nbRepeats = 30, nbLeaves = 15, nb.move = 1, dist = 3, phyType = "Alt")

# Lancement manuel d'une instance ####
  random = createTaxPhyAlt(nbLeaves = 6, nb.move = 1, dist = 3)
  Taxo = random$taxo
  Phylo = random$phylo
  doneMat = matrix(nrow = Taxo$Nnode, ncol = Phylo$Nnode, dimnames = list(Taxo$node.label, Phylo$node.label))
  resultat = mast(taxo = Taxo, phylo = Phylo ,subRootTax = Taxo$node.label[1], subRootPhy = Phylo$node.label[1])
  gsub("NA,*|,*NA", "","1,2,NA")
  
       