# Qualitative immune infiltration analysis with deconvolution


# Packages
library(tidyverse)
library(ggpubr)
library(corrplot)
library(gdata)

# Data and masterfiles
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/processed/results_CIBERSORT_BC_BILAT_NEO_LM22_df_3.RData")
head(results_CIBERSORT_BC_BILAT_NEO_LM22_df_3)

load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/submission_R1_Nature/R2_nature/panels_tils_total.RData")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/submission_R1_Nature/R2_nature/panels_macrophages_total.RData")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/immunofluo/data/processed/panels_macrophages.RData")

head(panels_tils_total)

# Functions and colors
colsubtype_cvd_pCR3 <- c(  "PT"     = "#49AAE3",  "RD"     = "#57EBC3")

# We applied on the 20 samples of the cohort the CIBERSORT algorithm using the “absolute” mode to deconvolute RNAseq expression profiles
# into 22 subsets of immune subpopulations. 

# Cf how to preprocess it 
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/processed/CIBERSORT_BC_BILAT_NEO_LM22_3.RData")
CIBERSORT_BC_BILAT_NEO_LM22_3$all   <- "all"
CIBERSORT_BC_BILAT_NEO_LM22_3$value <- CIBERSORT_BC_BILAT_NEO_LM22_3$value*100
head(CIBERSORT_BC_BILAT_NEO_LM22_3)
nrow(CIBERSORT_BC_BILAT_NEO_LM22_3) # 440

p_boxplot_global_LM22_wp_abs <- CIBERSORT_BC_BILAT_NEO_LM22_3 %>% 
                                ggplot( aes(x=reorder(celltype,value), y=value, fill = celltype)) +
                                geom_boxplot()  + 
                                # geom_jitter()  + 
                                xlab('') + ylab('Absolute abundance score') + 
                                coord_flip(ylim = c(0,110))+
                                theme_bw() + 
                                theme(axis.text.x = element_text( hjust=1),
                                      legend.position = "bottom",
                                      legend.text=element_text(size=12),legend.title=element_text(size=12))
p_boxplot_global_LM22_wp_abs

p_boxplot_global_LM22_rel_PT_RD_abs    <- CIBERSORT_BC_BILAT_NEO_LM22_3 %>% ggplot( aes(x=reorder(celltype,value), y=value, fill = tumortype)) +
                                          geom_boxplot()  + 
                                          geom_jitter(size=0.5)  + 
                                          xlab('') + ylab('Absolute abundance score') + 
                                          coord_flip(ylim = c(0,110))+
                                          theme_bw() + 
                                          stat_compare_means( label = "p.format",label.y = 100) +
                                          scale_fill_manual(name = "", values = colsubtype_cvd_pCR3, drop = TRUE)  +
                                          theme(axis.text.x = element_text(hjust=1),
                                                legend.position = "bottom",
                                                legend.text=element_text(size=12),legend.title=element_text(size=12))
p_boxplot_global_LM22_rel_PT_RD_abs

p_boxplot_global_LM22_wp_PT_RD_abs <- plot_grid(p_boxplot_global_LM22_wp_abs + theme(legend.position = "none") , 
                                                p_boxplot_global_LM22_rel_PT_RD_abs  + 
                                                  theme(legend.position = "right", axis.ticks.y   = element_blank(),
                                                        axis.text.y = element_blank()), 
                                                labels = "AUTO", rel_widths = c(6,5))     
p_boxplot_global_LM22_wp_PT_RD_abs

save_plot(file= "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/FigS8_p_boxplot_global_LM22_wp_PT_RD_abs.pdf", 
          p_boxplot_global_LM22_wp_PT_RD_abs, nrow=1, ncol=2,  base_height=6.5, base_width=6)							 

# The top 3 most abundant immune subpopulations were M2 macrophages, 
# CD4 memory resting T cells and M1 macrophages (FigS8A). 
# CD4 memory T cells and M2 macrophages were increased in RD compared with PT (FigS8B). 


# To further characterize the microenvironment of these samples using an orthogonal experimental approach, we performed immunofluorescence stainings 
# using 2 antibodies panels (Panel1: CD8 CD45RO CD20 CD4 and FoxP3; Panel 1: CD4 – CD8 – CD45RO – FOXP3 – CD20 – pan-cytokeratin; FigS8A-F); 
# Panel 2: CD68 – CD163 – CD138 – CKIT – pan-cytokeratin; FigS8G-H) to assess the concordance between the immune subpopulations inferred by gene expression
# and the number of cells stained on pathologic slides. 

# Cytotoxic CD8+
p_Somme_de_CD8pos_Cytotox_T_cells <- ggscatter(panels_tils_total, x = "Somme_de_CD8pos", y = "Cytotox_T_cells",
                                               add = "reg.line", shape=21, size = 2.5, fill = "lightgray",
                                               add.params = list(color = "#BE9A32", fill = "#BE9A32"), 
                                               conf.int = TRUE) + stat_cor(method = "pearson") + theme_bw() + 
                                               xlab("Number of CD8+ cells") + ylab("Cytotoxic T cells \n Absolute abundance score")

# B cells : CD 20+ 
p_Somme_de_CD20pos_B_cells_naive <- ggscatter(panels_tils_total, x = "Somme_de_CD20pos", y = "B_cells_naive",
                                              add = "reg.line",shape=21, size = 2.5,fill = "lightgray",
                                              add.params = list(color = "#F9766D", fill = "#F9766D"), conf.int = TRUE) +
                                              stat_cor(method = "pearson") + theme_bw() + 
                                              xlab("Number of CD20+ cells") + ylab("B cells naive \n Absolute abundance score")

# T reg  CD4+ / FOXP3+ 
p_Somme_de_CD4pos_FOXP3pos_T_reg <- ggscatter(panels_tils_total, x = "Somme_de_CD4pos_FOXP3pos", y = "T_reg",
                                              add = "reg.line", shape=21, size = 2.5,fill = "lightgray",
                                              add.params = list(color = "#53B867", fill = "#53B867"), conf.int = TRUE) + 
                                              stat_cor(method = "pearson") +   theme_bw() +
                                              theme_bw() + xlab("Number of CD4+/FOXP3+ cells") + ylab("T Regs \n Absolute abundance score")
#  Effector memory CD45RO+
p_Somme_de_CD45R0pos_T_cells_CD4_memory_activated <- ggscatter(panels_tils_total, x = "Somme_de_CD45R0pos", y = "T_cells_CD4_memory_activated",
                                             add = "reg.line",shape=21, size = 2.5,fill = "lightgray",
                                             add.params = list(color = "#76AE33", fill = "#76AE33"),conf.int = TRUE) + 
                                              stat_cor(method = "pearson") + theme_bw() + 
                                              xlab("Number of CD45RO+ cells") + ylab("T Cells memory activated \n Absolute abundance score")


# CD68+/CD163+ and Macrophages M2
p_Somme_de_CD68pos_CD163pos_MacroM2 <- ggscatter(panels_macrophages_total, x = 'Somme_de_cd68pos/cd163pos', y = "Macrophages_M2",
                                                 add = "reg.line",  shape=21, size = 2.5,fill = "lightgray",
                                                 add.params = list(color = "#7B95FF", fill = "#7B95FF"), conf.int = TRUE) + 
                                                stat_cor(method = "pearson") +   theme_bw() + 
                                                xlab("Number of CD68+/CD163+ cells") + ylab("Macrophages M2 \n Absolute abundance score")

# CD138+ and PLasma cells M2
p_Somme_de_CD138pos_plasma_cells <- ggscatter(panels_macrophages_total, x = 'Somme_de_cd138pos', y = "Plasma_cells",
                                              add = "reg.line", shape=21, size = 2.5,fill = "lightgray",
                                              add.params = list(color = "#DB8E01", fill = "#DB8E01"),conf.int = TRUE) +
                                              stat_cor(method = "pearson") + theme_bw() + 
                                              xlab("Number of CD138+ cells") + ylab("Plasma cells \n Absolute abundance score")

# CD68+/CD163+ and Macrophages M2
p_Somme_de_CD68pos_MacroM0_M1_M2 <- ggscatter(panels_macrophages_total, x = 'Somme_de_cd68pos', y = "Macrophages_M0_M1_M2",
                                              add = "reg.line", shape=21, size = 2.5,fill = "lightgray",
                                              add.params = list(color = "#BE9A32", fill = "#BE9A32"), conf.int = TRUE) + 
                                              stat_cor(method = "pearson") +   theme_bw() + 
                                              xlab("Number of CD68+ cells") + ylab("Macrophages M0, M1 and M2 \n Absolute abundance score")

# CKIT+/CK- cells  ( Only on CK- for mast cells)
p_Somme_de_ckitpos_mast_cells <- ggscatter(panels_macrophages %>% filter(Tissue_Category == "CK-"),
                                           x = 'Somme_de_ck-/ckitpos', y = "Mast_cells_resting",
                                           add = "reg.line",  shape=21, size = 2.5, fill = "lightgray",
                                           add.params = list(color = "#EF67EB", fill = "#EF67EB"), conf.int = TRUE) + 
                                            stat_cor(method = "pearson") + theme_bw() + 
                                            xlab("Number of CKIT+/CK- cells") + ylab("Mast cells resting \n Absolute abundance score")
p_Somme_de_ckitpos_mast_cells

p_correl_cibersort_immunoflu_compil_8panels <- plot_grid(     p_Somme_de_CD8pos_Cytotox_T_cells ,
                                                              p_Somme_de_CD4pos_FOXP3pos_T_reg,
                                                              p_Somme_de_CD45R0pos_T_cells_CD4_memory_activated,
                                                              p_Somme_de_CD20pos_B_cells_naive,
                                                              p_Somme_de_CD68pos_CD163pos_MacroM2 ,
                                                              p_Somme_de_CD68pos_MacroM0_M1_M2,
                                                              p_Somme_de_CD138pos_plasma_cells,
                                                              p_Somme_de_ckitpos_mast_cells + coord_cartesian(xlim = c(0,410)),
                                                              labels = c("I","J","K","L","M","N","O","P") ,
                                                              ncol = 2)
p_correl_cibersort_immunoflu_compil_8panels

save_plot(p_correl_cibersort_immunoflu_compil_8panels, 
          file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig7_p_correl_cibersort_immunoflu_compil_8panels.pdf", 
          base_height = 13, base_width = 9)

# The correlation coefficient between both metrics was statistically significant 
# regarding cytotoxic T cells (CD8+ cells), T regs (CD4+/FOXP3+ cells), Mast cells resting (CKIT+/CK- cells) (ExtendedFig7I-J-Q), 
# was nearly significant for T cells memory activated (CD45RO+ cells) and B cells naïve (CD20+) (ExtendedFig7K-L);
# and was not significant in macrophages M2 (CD68+/CD163+ cells) or macrophages of any type (CD68+ cells) and plasma cells (CD138+ cells)(ExtendedFig7M-N-O). 


# We average Patient6 for each
mat_patient_6_L_R <- CIBERSORT_BC_BILAT_NEO_LM22_3 %>% 
  filter(tumortype == "PT" ) %>% 
  select(samplename,patient,patient_sideAB,tumortype,celltype,side,value) %>% filter(patient == "Patient6") %>% group_by(patient,celltype,side) %>% summarise(value = mean(value))
head(mat_patient_6_L_R)

mat_patient_no6_L_R <- CIBERSORT_BC_BILAT_NEO_LM22_3 %>% 
  filter(tumortype == "PT" ) %>% 
  select(samplename,patient,patient_sideAB,tumortype,celltype,side,value) %>% filter(patient != "Patient6") %>% select(patient, celltype, side, value)
head(mat_patient_no6_L_R)

mat_all_patients_cibersort <- rbind(mat_patient_6_L_R,mat_patient_no6_L_R)

p_L_R_lines_cibersort  <-   ggplot(mat_all_patients_cibersort  ,
                                   aes	(x=celltype,y=value, color =celltype, shape = side)) +
                            theme_bw()+
                            theme(axis.ticks.x = element_blank(), 
                                  legend.position = "bottom")+
                            facet_wrap( . ~ patient)   + geom_point(size = 2) +  
                            coord_flip()+
                            geom_line(aes	(x=celltype,y=value,group=celltype, colour =celltype),
                                      size = 1,
                                      position = position_dodge(width=0.1))	+
                            scale_shape_manual(values = c(8,17), name = " ", labels =  c("Left","Right")) +
                            guides(colour = FALSE) +
                            ylab("Absolute abundance score")+ xlab("") + 
                            ggtitle("Difference between left and right tumors regarding immune subsets levels")
p_L_R_lines_cibersort

# Create the plot with arrows for PT and RD
collevelside                  <- c("R"="red3","L"="limegreen")
coltumortype                  <- c("PT"="#dd9d03", "RD"="#f6ce48")
coltumortype_shape            <- c(coltumortype,collevelside)
coltumortype_shape["PT"]      <- c(15)
coltumortype_shape["RD"]      <- c(4)
coltumortype_shape["R"]       <- c(8)
coltumortype_shape["L"]       <- c(11)

p_PT_RD_lines_horiz_cibersort  <- ggplot(CIBERSORT_BC_BILAT_NEO_LM22_3 %>%
                                        filter(tumortype2 %in% c("PT_with_RD","RD"))  %>% 
                                        select(samplename,patient,patient_sideAB,tumortype,celltype,side,value)  ,
                                      aes	(x=celltype,y=value, color =celltype, shape = tumortype)) +
                                      theme_bw()+
                                      theme(axis.ticks.x = element_blank(), 
                                            legend.position = "bottom")+
                                      facet_wrap( . ~ patient_sideAB)   + geom_point(size = 2) +  
                                      coord_flip()+
                                      geom_line(aes	(x=celltype,y=value,group=celltype, colour =celltype),
                                                size = 1,
                                                arrow = arrow(length=unit(0.2,"cm") ,type = "closed"  , angle = 35),
                                                position = position_dodge(width=0.1))	+
                                      scale_shape_manual(values = c(1,5), name = " " , labels = c("Primary tumor (PT)","Residual disease (RD)")) +
                                      guides(colour = FALSE) + 
                                      ylab("Absolute abundance score") + xlab("") + 
                                      ggtitle("Difference between PT and RD regarding immune subsets levels")
p_PT_RD_lines_horiz_cibersort

p_lines_horiz_compil_L_R_PT_RD     <-   plot_grid(p_L_R_lines_cibersort ,
                                                  p_PT_RD_lines_horiz_cibersort + xlab(""), 
                                                  labels = c("AUTO"),
                                                  ncol = 2, nrow = 1, 
                                                  align="hv")

save_plot(file="/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig8.pdf", 
          p_lines_horiz_compil_L_R_PT_RD, 
          base_height=9, base_width=17)							 

# At the patient level, the immune composition of the paired contralateral tumors were different regarding several immune subsets 
# (Macrophages M0, M1, M2, T cells CD4 memory activated and resting) (Extended Fig8A), while the variation of the immune composition between PT and RD 
# mostly concerned increasing levels of macrophages M2, M0 and T cells CD4 memory resting (Extended Fig8B). 


# See previous code
# preprocess_correlation_matrix_cibersort_LM22.r
mat_wide          <- results_CIBERSORT_BC_BILAT_NEO_LM22_df_3 %>%  rename( sample = rowname) %>% 
                                                            select(-Correlation, - RMSE, -sig.score, - Pvalue)

levelsample_corr <- c(  "PT1_R","PT1_L","PT2_R","PT2_L","PT3_R","PT3_L","PT4_R","PT4_L","PT5_R","PT5_L",
                        "PT6A_R","PT6B_R","PT6A_L","PT6B_L",
                        "RD1_L","RD3_R","RD3_L","RD4_R","RD5_R","RD6B_R")
mat_wide$sample  <- factor(mat_wide$sample, levels = levelsample_corr)
mat_wide_tmp     <- mat_wide[match(mat_wide$sample,levelsample_corr),]
mat_wide         <- mat_wide_tmp
rep.row<-function(x,n){
                          matrix(rep(x,each=n),nrow=n)
                        }
levelsample      <- mat_wide$sample

mat_sum_squared           <- matrix(NA, ncol = 20, nrow = 20) %>% as.data.frame()
colnames(mat_sum_squared) <- levelsample
rownames(mat_sum_squared) <- levelsample

for (i in 1 : length(levelsample)){
  # i=1
  print(i)
  tmp_pat     <- levelsample[i]; print(tmp_pat)
  tmp_mat_pat <- mat_wide %>% filter(sample == tmp_pat) %>% select(-sample) %>% as.matrix()
  tmp_mat_pat_unnamed <- as.vector(unname(tmp_mat_pat))
  tmp_mat_tmp_pat <- matrix(rep(tmp_mat_pat_unnamed,20), nrow = 20,byrow = TRUE) %>% as.data.frame()
  M1 <- mat_wide %>% select(-sample)
  M2 <- tmp_mat_tmp_pat
  tmp_mat_pat_soustract                      <-  M1 - M2
  tmp_mat_pat_soustract_squared              <-  tmp_mat_pat_soustract*tmp_mat_pat_soustract 
  sum_all_difference_squared_tmp_sample      <- rowSums(tmp_mat_pat_soustract_squared)
  mat_sum_squared[i,]                        <- sum_all_difference_squared_tmp_sample
}
mat_sum_squared <- mat_sum_squared %>% as.matrix()
head(mat_sum_squared)
# mat_sum_squared  <- mat_sum_squared*100

# Create correlation matrix 

colxx <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig9_correlogram_cibersort.pdf", height = 10, width = 15)
corrplot(mat_sum_squared, 
         type="upper", 
         # method="color",
         is.corr = FALSE,addCoef.col = TRUE,
         # title = "toto",
         cl.lim=c(0,1.5),diag = F,
         # col=colxx(200),
         # col=colorRampPalette(c("white","lightpink", "red","brown"))(100),
         cl.pos = "n", number.cex=0.75,tl.col="black") # ,) #
dev.off()

# We compared the predicted immune contexture patterns between samples of the cohort using a dissimilarity index (Idissimilarity). 
# The higher the dissimilarity index, the more the composition of the immune infiltration differs. 

# Overall matrix
all        <- mat_sum_squared %>% upperTriangle()
all %>% mean() # 0.337

# Pairs left and right (individual)
pairs_L_R  <- c(  mat_sum_squared["PT1_L","PT1_R"],
                  mat_sum_squared["PT2_L","PT2_R"],
                  mat_sum_squared["PT3_L","PT3_R"],
                  mat_sum_squared["PT4_L","PT4_R"],
                  mat_sum_squared["PT5_L","PT5_R"]) 
pairs_L_R %>% mean() # 0.2228

# Pairs PT and RD (individual)
pairs_PT_RD <- c(   mat_sum_squared["PT1_L","RD1_L"],
                    mat_sum_squared["PT3_R","RD3_R"],
                    mat_sum_squared["PT3_L","RD3_L"],
                    mat_sum_squared["PT4_R","RD4_R"],
                    mat_sum_squared["PT5_R","RD5_R"],
                    mat_sum_squared["PT6B_R","RD6B_R"]) 
pairs_PT_RD %>% mean() # 0.29

# Neither the mean dissimilarity indices
# (Extended Fig9) of paired left and right tumor (green bordered squares, mean Idissimilarity= 0.22) nor paired PT and RD 
# (yellow bordered squared, mean Idissimilarity= 0.29) were statistically different from the rest of the samples. 

# Matrix PT versus PT
PT_vs_PT <- mat_sum_squared[grep("PT", rownames(mat_sum_squared)),grep("PT",colnames(mat_sum_squared))]  %>% upperTriangle()
PT_vs_PT %>% mean() # 0.24

# Matrix RD versus RD
RD_vs_RD <- mat_sum_squared[grep("RD", rownames(mat_sum_squared)),grep("RD",colnames(mat_sum_squared))]  %>% upperTriangle()
RD_vs_RD %>% mean() # 0.37
wilcox.test(PT_vs_PT,RD_vs_RD)   # <0.01

# Matrix PT versus RD
PT_vs_RD  <- mat_sum_squared[grep("PT", rownames(mat_sum_squared)),grep("RD",colnames(mat_sum_squared))]  %>% upperTriangle()
PT_vs_RD %>% mean() # 0.49

# At the cohort level, the dissimilarity was lower among the PT samples (blue area) than among the RD samples (orange area) 
# (mean Idissimilarity: 0.24 versus 0.37, p=0.009), while the greatest difference was seen between PTs samples compared with RDs samples
# (yellow area, mean Idissimilarity = 0.49). 
# These results suggest that the composition of the immune microenvironment is strongly associated with the pretreated or non-pretreated character of the sample 
# (PT or RD), in line with the changes in the immune contexture induced by NAC treatment. 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  