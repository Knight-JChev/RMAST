library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggprism)
library(ggpubr)

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
all_clusters = tibble(all_clusters) %>% filter(!is.na(V1))

all_kast$similarity = as.factor(all_kast$similarity)
all_kast$barcode = as.factor(all_kast$barcode)
all_kast = tibble(all_kast)

# Manipulation tableaux clusters ####
# Clusters > 20 seq
tmp = all_clusters %>%
  mutate(size = if_else(V1>=20, "Large cluster", "Small cluster")) %>%
  group_by(similarity, barcode, size) %>%
  summarise(nb_seq = sum(V1)) # Nb de séquences

clusters_sizes = all_clusters %>%
  mutate(size = if_else(V1>=20, "Large cluster", "Small cluster")) %>%
  group_by(similarity, barcode, size) %>%
  count(similarity, barcode, size, name="cluster_count") %>% # Nb de clusters
  left_join(tmp) # Join les tables

# Dot plot
all_clusters %>% dplyr::count(V1, barcode, similarity) %>%
  ggplot() +
  facet_grid(barcode~similarity) +
  geom_point(aes(x = V1, y = n, colour = barcode)) +
  geom_vline(aes(xintercept = 20), colour = "black", linetype = 2) +
  labs(x="Taille du cluster",
       y = "Nombre de clusters") +
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 9),
        strip.text = element_text(face = 2)) +
  coord_cartesian(ylim=c(1,1e+5))+
  scale_color_brewer(palette = "Set2")+
  scale_x_log10()+
  scale_y_log10()

# Autre plot pour voir les évolutions
ggplot(clusters_sizes)+
  geom_line(aes(x = as.numeric(similarity), y = nb_seq, colour = barcode),
            linewidth = 1.5)+
  facet_grid(size~barcode)+
  labs(x = "Similarité",
       y = "Nombre de séquences")+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 9),
        strip.text.y = element_text(size = 10),
        strip.text = element_text(face = 2)) +
  scale_x_continuous(breaks = c(0.85,0.9,0.95,0.99))+
  scale_color_brewer(palette = "Set2")+
  scale_y_log10()


ggplot(clusters_sizes)+
  geom_line(aes(x = as.numeric(similarity), y = cluster_count, colour = barcode),
            linewidth = 1.5)+
  facet_grid(size~barcode)+
  labs(y = "Nombre de clusters",
       x = "Similarité") +
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 9),
        strip.text.y = element_text(size = 10),
        strip.text = element_text(face = 2)) +
  scale_x_continuous(breaks = c(0.85,0.9,0.95,0.99))+
  scale_color_brewer(palette = "Set2")+
  scale_y_log10()

# KAST et MAST ####
# Kast et mast normalisés
all_kast = all_kast %>% 
  mutate(real_tree_size = if_else(tree_size>150, true =150, false = tree_size)) %>%
  mutate(norm_mast = mast/real_tree_size) %>%
  mutate(norm_kast = kast/real_tree_size)

# Nombre de clusters analysés
counts = all_kast %>% 
  group_by(similarity, barcode) %>%
  tally(name = "count") %>%
  mutate(count = paste0("n = ", count))

# Garde le max du mast normalisé pour annotation figure
pos_mast = all_kast %>% 
  group_by(similarity, barcode) %>%
  filter(norm_mast == max(norm_mast)) %>%
  distinct(norm_mast) %>%
  mutate(pos_mast = max(norm_mast)+0.05) %>%
  select(similarity, barcode, pos_mast)

# Garde le max du kast normalisé pour annotation figure + join les tableaux
kast_mast_counts = all_kast %>% 
  group_by(similarity, barcode) %>%
  filter(norm_kast == max(norm_kast)) %>%
  distinct(norm_kast) %>%
  mutate(pos_kast = max(norm_kast)+0.03) %>%
  select(similarity, barcode, pos_kast) %>%
  left_join(pos_mast) %>%
  left_join(counts)

# Violin plots
## Mast
MASTplot = ggplot(all_kast) +
  facet_grid(.~similarity) + 
  geom_violin(aes(x = barcode, y = norm_mast, fill = barcode),
              draw_quantiles = c(0.25,0.5,0.75), show.legend = F) +
  geom_text(data = kast_mast_counts, aes(x = barcode, y = -0.15, label = count),
            size=4)+
  labs(x = "Barcode",
       y = "MAST normalisé")+
  theme_bw()+
  coord_cartesian(ylim=c(0,0.8),
                  clip = "off")+ # Don't crop text out of plot
  theme(axis.title.x.bottom = element_blank(),
        axis.ticks.length.x.bottom = unit(.2,"cm"),
        axis.text.x = element_blank(),
        plot.margin = unit(c(1,1,2,1), "lines"),
        axis.text = element_text(size = 12 ),
        axis.title = element_text(size = 14)) + #Widens margins
  scale_fill_brewer(palette="Set2")
MASTplot

KASTplot= ggplot(all_kast) +
  facet_grid(.~similarity) + 
  geom_violin(aes(x = barcode, y = norm_kast, fill = barcode),
              draw_quantiles = c(0.25,0.5,0.75), show.legend = F) +
   labs(x = "Barcode",
       y = "KAST normalisé") +
  theme_bw()+
  coord_cartesian(ylim=c(0,0.8), clip="off")+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size = 12 ),
        axis.title = element_text(size = 14))+
  scale_x_discrete(labels = c("Inse01","Sper01","Vert01"))+
  scale_fill_brewer(palette="Set2")+
   annotation_ticks(sides="t", type = "major", outside = T,
                    linewidth = 0.4)
  KASTplot 

ggarrange(MASTplot,NULL, KASTplot, ncol = 1, nrow = 3,
          heights = c(1,-0.05,1))


## Tests ####
