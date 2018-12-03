setwd("C:/Users/Raphaëlle/Desktop/Cours/Cours APT 3A/Projet_Fil_Rouge")
install.packages('pillar')
install.packages('dplyr')
install.packages('tidyr')

library(dplyr)
library(tidyr)

Germ_A_T = read.csv('Germination_A_T.csv', sep=';', stringsAsFactors = FALSE)
Germ_A_M = read.csv('Germination_A_M.csv', sep=';', stringsAsFactors = FALSE)