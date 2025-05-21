library(dplyr)
library(ggplot2)
library(stringr)

resdir ="/home/knight/Bureau/Rouen_M1/Stage_LECA/KAST/Results/"

# Vecteur de seuils de similarités
comptes = list.files(resdir, pattern = "0\\.[0-9]*\\.tsv$")
kast_tables = list.files(resdir, pattern="max\\.tsv")

# Initialise Data frames
all_clusters = data.frame()
all_kast = data.frame()

# Lecture des données ####
# Build the table with cluster distribution
for (file in comptes) {
  clusters = as.data.frame(t(read.table(paste0(resdir,file), 
    sep="\t", header = F, colClasses=c("integer"))))
  
  barcode = str_extract(file, "DB\\.(.+?)\\.", group = 1)
  thr = str_extract(file, "(0\\.[0-9]+)\\.tsv", group = 1)
    
  clusters = clusters %>%
    mutate(similarity = thr, barcode = barcode) # Add columns for later grouping
  
  all_clusters = rbind(all_clusters, clusters)
}

# Build the table with kast results
for (file in kast_tables){
  kast = read.table(paste0(resdir, file), sep="\t",
    header = T, row.names=1)
   
  barcode = str_extract(file, "table_(.+?)\\.", group = 1)
  thr = str_extract(file, "_(0\\.[0-9]+)_", group = 1)
  
  kast = kast %>%
    select(!c(taxo,phylo)) %>% # Remove trees
    mutate(similarity = thr, barcode = barcode) # Add columns for later grouping
  
  all_kast = rbind(all_kast, kast)
}

# Conversion to factors
all_clusters$similarity = as.factor(all_clusters$similarity)
all_clusters$barcode = as.factor(all_clusters$barcode)

all_kast$similarity = as.factor(all_kast$similarity)
all_kast$barcode = as.factor(all_kast$barcode)

# Plots ####
# Dot plot
all_clusters %>% dplyr::count(V1, barcode, similarity) %>%
  ggplot() +
  facet_grid(barcode~similarity) +
  geom_point(aes(x = V1, y = n, colour = similarity)) +
  geom_vline(aes(xintercept = 20), colour = "red") +
  geom_vline(aes(xintercept = 150), colour = "red") +
  labs(x="Cluster size",
       y = "Nombre d'observation") +
  scale_x_log10()+
  scale_y_log10()


# Density plot à réfléchir  
all_clusters %>% dplyr::count(V1,barcode, similarity) %>%
  ggplot() + 
  facet_grid(barcode~.) +
  geom_density(aes(x = V1, colour = similarity))+
  labs(x = "Cluster size") +
  scale_x_log10()

# Manipulation de tableaux ####
# Nombre de clusters de taille < 20
all_clusters %>%
  group_by(similarity, barcode) %>%
  count(V1, barcode, similarity, name = "count") %>%
  filter(V1<20) %>%
  tally(count)

# Nombre de clusters analysés
all_kast %>%
  group_by(similarity, barcode) %>%
  count(similarity, barcode)
