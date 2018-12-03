setwd("C:/Users/Raphaëlle/Desktop/Cours/Cours APT 3A/Projet_Fil_Rouge")
install.packages('pillar')
install.packages('dplyr')
install.packages('tidyr')

library(dplyr)
library(tidyr)

Germ_E_T = read.csv('Germination_E_T.csv', sep=';', stringsAsFactors = FALSE)
Germ_E_M = read.csv('Germination_E_M.csv', sep=';', stringsAsFactors = FALSE)