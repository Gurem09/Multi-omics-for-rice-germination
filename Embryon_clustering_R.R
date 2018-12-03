setwd("C:/Users/guill/Documents/APT/3A/Cours/Projet_fil_rouge_germination/Données/CSV-20181119T140830Z-001/CSV")

#Ressources : https://uc-r.github.io/kmeans_clustering

library(dplyr)
library(tidyr) 

library(tidyverse)
library(cluster)
library(factoextra)

############# Analyse base transcriptomes
Germ_E_T = read.csv('Germination_E_T.csv', sep=';', stringsAsFactors = FALSE)
Germ_E_T_pb_gn =
  Germ_E_T %>%
  select(probe,MSU_id)
#préparation des bases
Germ_E_T_R1 =
  Germ_E_T %>%
  select( E0_R1,E4_R1,E8_R1,E12_R1,E16_R1,E24_R1)
row.names(Germ_E_T_R1) = Germ_E_T[,1]
Germ_E_T_R2 =
  Germ_E_T %>%
  select( E0_R2,E4_R2,E8_R2,E12_R2,E16_R2,E24_R2)
row.names(Germ_E_T_R2) = Germ_E_T[,1]
Germ_E_T_R3 =
  Germ_E_T %>%
  select( E0_R3,E4_R3,E8_R3,E12_R3,E16_R3,E24_R3)
row.names(Germ_E_T_R3) = Germ_E_T[,1]

#Clustering k_means
center_number = 4
k1 = kmeans(Germ_E_T_R1, centers = center_number)
str(k1)
k2 = kmeans(Germ_E_T_R2, centers = center_number)
str(k2)
k3 = kmeans(Germ_E_T_R3, centers = center_number)
str(k3)

#Plot
p1 <- fviz_cluster(k1, geom = "point", data = Germ_E_T_R1) + ggtitle("R1")
p2 <- fviz_cluster(k2, geom = "point",  data = Germ_E_T_R2) + ggtitle("R2")
p3 <- fviz_cluster(k3, geom = "point",  data = Germ_E_T_R3) + ggtitle("R3")

library(gridExtra)
grid.arrange(p1, p2, p3, nrow = 2)

#  Affectation des clusters
library(tibble)
R1_cluster=as.data.frame(k1$cluster)
R1_cluster=rownames_to_column(R1_cluster, var="probe")
R2_cluster=as.data.frame(k2$cluster)
R2_cluster=rownames_to_column(R2_cluster, var="probe")
R3_cluster=as.data.frame(k3$cluster)
R3_cluster=rownames_to_column(R3_cluster, var="probe")
R1vsR2vsR3_cluster =
  R1_cluster %>%
  full_join(R2_cluster, by="probe")%>%
  full_join(R3_cluster, by="probe")%>%
  full_join(Germ_E_T_pb_gn, by="probe")%>%
  mutate(k1vsk2vsk3 = case_when(
    k1$cluster==k2$cluster & k2$cluster==k3$cluster ~ 1,
    TRUE ~ 0
  ))
  
############# Analyse base métabolomes :
Germ_E_M = read.csv('Germination_E_M.csv', sep=';', stringsAsFactors = FALSE)
cmd_to_elim = c("Fructose","Raffinose","Alanine","Sucrose","Glucose","Pyroglutamate", "Phosphate")
Germ_E_M_1 =
  Germ_E_M %>%
  filter(!Compound %in% cmd_to_elim)%>%
  group_by(Class,Compound)%>%
  summarize_all(funs(mean))%>%
  ungroup()%>%
  select( Compound,E0.1,E4.1,E8.1,E12.1,E16.1,E24.1)
#row.names(Germ_E_M_1) = Germ_E_M_1[,1]
Germ_E_M_1[,2:7]=scale(Germ_E_M_1[,2:7])
center_number = 11
km = kmeans(Germ_E_M_1[,2:7], centers = center_number)
str(km)
fviz_cluster(km, geom = "point", data = Germ_E_M_1[,2:7])

m1_cluster=as.data.frame(km$cluster)
Germ_E_M_2 =
  Germ_E_M_1 %>%
  bind_cols(m1_cluster)%>%
  select("Compound", "km$cluster")