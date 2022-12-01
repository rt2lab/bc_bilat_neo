################################################## 
#                       SETUP
################################################## 

# Packages
library(ggpubr)

# Data and masterfiles
load("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/data/raw/CD3_manual_and_digital_de_ident.RData")
head(CD3_manual_and_digital_de_ident)

p_str_tils_versus_CD3_pos <- ggscatter(CD3_manual_and_digital, x = "CD3_pos", y = "str_til_perc",
                                       add = "reg.line",  shape=21,size = 4,fill = "#619CFF",
                                       add.params = list(color = "#619CFF", fill = "#619CFF"), conf.int = TRUE) +
                                        stat_cor(method = "pearson") + coord_cartesian(ylim = c(0,75))+
                                        theme_bw() + xlab("Number of CD3 + cells  ") + ylab("Str TIL levels (%), manual assessment")

p_str_tils_versus_CD3_pos_tissue <- ggscatter(CD3_manual_and_digital, x = "% CD3 Positive Tissue", y = "str_til_perc",
                                              add = "reg.line", shape=21, size = 4,fill = "#619CFF",
                                              add.params = list(color = "#619CFF", fill = "#619CFF"), conf.int = TRUE) + 
                                              stat_cor(method = "pearson") + coord_cartesian(ylim = c(0,75))+
                                              theme_bw() + xlab("% CD3 Positive Tissue") + ylab("Str TIL levels (%), manual assessment")

p_it_tils_versus_CD3_pos <- ggscatter(CD3_manual_and_digital, x = "CD3_pos", y = "it_til_perc",
                                      add = "reg.line",shape=21,size = 4,fill = "#CFE1FF",
                                      add.params = list(color = "#CFE1FF", fill = "#CFE1FF"), conf.int = TRUE) + 
                                      stat_cor(method = "pearson") + coord_cartesian(ylim = c(0,65))+
                                      theme_bw() + xlab("Number of CD3 + cells  ") + ylab("IT TIL levels (%), manual assessment")

p_it_tils_versus_CD3_pos_tissue <- ggscatter(CD3_manual_and_digital, x = "% CD3 Positive Tissue", y = "it_til_perc",
                                             add = "reg.line", shape=21, size = 4,fill = "#CFE1FF",
                                             add.params = list(color = "#CFE1FF", fill = "#CFE1FF"), conf.int = TRUE) + 
                                            stat_cor(method = "pearson") + coord_cartesian(ylim = c(0,65))+
                                            theme_bw() + xlab("% CD3 Positive Tissue") + ylab("IT TIL levels (%), manual assessment")

p_compil_tils_cd3 <- plot_grid(p_str_tils_versus_CD3_pos,p_str_tils_versus_CD3_pos_tissue,
                               p_it_tils_versus_CD3_pos,p_it_tils_versus_CD3_pos_tissue,
                               nrow=2, ncol = 2,
                               labels = c("D","E","F","G"))

save_plot(p_compil_tils_cd3 , file = "~/RT2Lab/bc_bilat_neo_git/figures/FigS9_p_compil_tils_cd3.pdf",
          base_height = 8, base_width = 8)
