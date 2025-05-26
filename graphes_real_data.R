library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

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

all_kast %>%
  ggplot() +
  facet_grid(.~barcode) +
  geom_boxplot(aes(y = mast, x = similarity ))

# Manipulation de tableaux ####
# Nombre de clusters de taille < 20
all_clusters %>%
  group_by(similarity, barcode) %>%
  count(V1, barcode, similarity, name = "count") %>%
  filter(V1>20) %>%
  tally(count)

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
  scale_fill_brewer(palette="Spectral")
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
  scale_fill_brewer(palette="Spectral")+
   annotation_ticks(sides="t", type = "major", outside = T,
                    linewidth = 0.4)
  KASTplot 

ggarrange(MASTplot,NULL, KASTplot, ncol = 1, nrow = 3,
          heights = c(1,-0.05,1))


## Tests ####
# Generate data
df <- data.frame(y=c("cat1","cat2","cat3"),
                 x=c(12,10,14),
                 n=c(5,15,20))

# Create the plot
ggplot(df,aes(x=x,y=y,label=n)) +
  geom_point()+
  geom_text(y = I(4), # Set text's position to the right end of the plot
            hjust = 0,
            size = 8) +
  coord_cartesian(xlim = c(10, 14), # This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  theme(plot.margin = unit(c(1,3,1,1), "lines")) # This widens the right margin
