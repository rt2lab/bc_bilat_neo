source('/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/all_Rmd_cliniques/src/setup_BC_BILAT_NEO_clinique.R', local  = TRUE)

library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)

load(file="~/RT2Lab/BC_BILAT_NEO/clinique/data/processed/unilat_and_bilat.RData")
nrow(unilat_and_bilat) # 18190

# And only for unique patient
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/unilat_and_bilat_unique_pat.Rdata")
nrow(unilat_and_bilat_unique_pat) # 17575