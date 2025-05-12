library(dplyr)
library(ggplot2)

resdir ="/home/knight/Bureau/Rouen_M1/Stage_LECA/Vert01/Results"

sim_thr =c(0.75, 0.8, 0.85, 0.90, 0.95)
thr = 0.95
all_clusters = data.frame()
all_kast = data.frame()

for (thr in sim_thr) {
  clusters = as.data.frame(t(read.table(
    paste0(resdir,"/DB.Vert01.fasta_comptes_cluster_", thr ,".tsv"), 
    sep="\t", header = F, colClasses=c("integer"))))
  
  clusters = clusters %>%
    mutate(similarity = thr)
  
  all_clusters = rbind(all_clusters, clusters)
  
  kast = read.table(
    paste0(resdir,"/kast_table_", thr,"_similarity.tsv"), sep="\t",
    header = T, row.names=1)
  kast = kast %>%
    mutate(similarity = thr)
  
  all_kast = rbind(all_kast, kast)
  
}
all_clusters$similarity = as.factor(all_clusters$similarity)
all_kast$similarity = as.factor(all_kast$similarity)


cluster_test = as.data.frame(all_clusters[10<all_clusters[,1],])

all_clusters %>% dplyr::count(V1, similarity) %>%
  ggplot() + geom_point(aes(x = log(V1), y = log(n))) +
  geom_vline(aes(xintercept = log(11)), colour = "red") +
  geom_vline(aes(xintercept = log(150)), colour = "red")

all_clusters %>% dplyr::count(V1, similarity) %>%
  ggplot() + geom_point(aes(x = V1, y = n, colour = similarity)) +
  geom_vline(aes(xintercept = 11), colour = "lightgreen") +
  geom_vline(aes(xintercept = 150), colour = "lightgreen")+
  labs(x="Cluster size",
       y = "Nombre d'observation") +
  scale_x_log10()+
  scale_y_log10()
  
all_clusters %>% dplyr::count(V1, similarity) %>%
  ggplot() + geom_density(aes(x = V1, colour = similarity))+
  labs(x = "Cluster size") +
  scale_x_log10()


class(as.data.frame(all_clusters))

