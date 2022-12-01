# Transcriptomic alterations

# Tumor clustering and principal component analysis (PCA)


# Load data and setup

library(Hmisc)

# clinical data
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/data_clin_min.RData")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/data_clin_min_new.RData")

# RNASeq data

# Raw
tpm <- read.delim("~/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer1/new_pca_with_patient5/featurecounts_tpm2.txt", comment.char="#")
colnames(tpm)[which(colnames(tpm) == "RD6_R")] <- "RD6B_R"
# save(tpm, file = "/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/tpm.Rdata")
head(tpm)


# mv genes 
# already preprocessed
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/mvgenes_counts.RData")
colnames(mvgenes_counts)[colnames(mvgenes_counts)=="RD6_R"]  <- "RD6B_R"
mvgenes_counts <- mvgenes_counts[,match(data_clin_min$samplename,colnames(mvgenes_counts))] # reorder genes  # 2846 20

# Functions and colors
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/trim_BS.r', local=TRUE)
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/fct.forclustering_CL_JA.R')					
my_palette 		<- colorRampPalette(c("green", "black", "red"))(n = 299)

colsubtype_cvd     <- c("luminal" = "#0072B2",
                        "TNBC"    = "#D55E00", 
                        "HER2+"   = "#009E73") #		[1] "#F8766D" # rouge "#00BA38" vert "#619CFF" bleu

colsubtype_cvd_pCR <- c("luminal" = "#0072B2",
                        "TNBC"    = "#D55E00", 
                        "HER2+"   = "#009E73",
                        "pCR"     = "#6ABBEB",
                        "No pCR"  = "#044970",
                        "PT"     = "#49AAE3",
                        "RD"     = "#57EBC3",
                        "different pCR status"  = "#F5AE78") 

collevelside2   <- c("R"="#D6F8D6","L"="#9DC3C2")
# collevelBRCA_2   <- c("BRCAmut"="#0B3954","BRCAwt"="#F4FFFD")
collevelBRCA_3   <- c("BRCAmut"="#0B3954","BRCAwt"="#ADD7F6")
collevelBRCA_4   <- c("BRCA2 mut"="#0B3954","BRCA1 mut"="#D64550","wt"="#ADD7F6")

levelsample <- c( "PT1_R","PT1_L","RD1_L",
                  "PT2_R","PT2_L",
                  "PT3_R","RD3_R","PT3_L","RD3_L",
                  "PT4_R","RD4_R","PT4_L",
                  "PT5_R","RD5_R","PT5_L",
                  "PT6A_R","PT6B_R","RD6B_R","PT6A_L","PT6B_L")

load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/data/raw/ashp_colors.RData")

# We performed unsupervised hierarchical clustering based on transcriptomic profile of the most variable genes. 

dataf2      <- mvgenes_counts 
head(mvgenes_counts)
cl.col 			<- clustering(dataf2,metric='pearson',method='ward')		# Samples
cl.row 			<- clustering(t(dataf2),metric='pearson',method='ward')		# genes
dataf 			<- trim.heatmap(dataf2,0.99)

# Cut into 4  clusters
cluster_genes_4 <- cutree(cl.row,4)

# And see enrichment on GSEA website

cluster4_genes_new4 <- rownames(dataf)[which(cluster_genes_4==1)] # Cluster 1 # 1107  ,
# HALLMARK_KRAS_SIGNALING_DN,HALLMARK_KRAS_SIGNALING_UP,HALLMARK_ESTROGEN_RESPONSE_LATE,HALLMARK_INFLAMMATORY_RESPONSE,
# HALLMARK_COAGULATION,HALLMARK_COMPLEMENT,HALLMARK_HYPOXIA,HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,HALLMARK_MYOGENESIS,HALLMARK_XENOBIOTIC_METABOLISM

cluster4_genes_new2 <- rownames(dataf)[which(cluster_genes_4==2)] # Cluster 2  617# 
#HALLMARK_TNFA_SIGNALING_VIA_NFKB# HALLMARK_MYOGENESIS# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION# HALLMARK_ADIPOGENESIS# 
# HALLMARK_HYPOXIA# HALLMARK_APICAL_JUNCTION# HALLMARK_XENOBIOTIC_METABOLISM# HALLMARK_FATTY_ACID_METABOLISM# HALLMARK_KRAS_SIGNALING_UP# HALLMARK_UV_RESPONSE_UP

cluster4_genes_new1 <- rownames(dataf)[which(cluster_genes_4==3)] # Cluster 3 734 
# HALLMARK_ESTROGEN_RESPONSE_EARLY,# HALLMARK_ESTROGEN_RESPONSE_LATE,# HALLMARK_KRAS_SIGNALING_DN,# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
# HALLMARK_MYOGENESIS,# HALLMARK_PANCREAS_BETA_CELLS,# HALLMARK_SPERMATOGENESIS,# HALLMARK_COAGULATION,# HALLMARK_UV_RESPONSE_DN,# HALLMARK_ANDROGEN_RESPONSE

cluster4_genes_new3 <- rownames(dataf)[which(cluster_genes_4==4)] # Cluster 4  390,
# HALLMARK_G2M_CHECKPOINT,HALLMARK_E2F_TARGETS,HALLMARK_MITOTIC_SPINDLE,HALLMARK_SPERMATOGENESIS,
# HALLMARK_MTORC1_SIGNALING,HALLMARK_ESTROGEN_RESPONSE_LATE,HALLMARK_ESTROGEN_RESPONSE_EARLY,HALLMARK_GLYCOLYSIS,HALLMARK_CHOLESTEROL_HOMEOSTASIS

# write.table(cluster4_genes_new1, file = "/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/processed/cluster4_genes_new1.txt",row.names = FALSE, col.names = FALSE)
# write.table(cluster4_genes_new2, file = "/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/processed/cluster4_genes_new2.txt",row.names = FALSE, col.names = FALSE)
# write.table(cluster4_genes_new3, file = "/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/processed/cluster4_genes_new3.txt",row.names = FALSE, col.names = FALSE)
# write.table(cluster4_genes_new4, file = "/Users/ahamypet/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/data/processed/cluster4_genes_new4.txt",row.names = FALSE, col.names = FALSE)


# Cluster New 1 (light blue) :  734 genes ESR1 - PGR - AR - TFF1 - GATA3 etc...
# Split les TNBC et les luminales
# high expression of most of the genes.
# To name cluster 1 : Estrogen responsive genes
# Highly relevant +++++

# Cluster new 2 (dark blue) : FOS JUN DUSP1  
# Fortement exprimé dans certaines luminales pour et une TNBC pour la quasi totalité des gènes
# Assez homogène pour tout le cluster
# EMT - TNF alpha - KRAS signaling
# HALLMARK_TNFA_SIGNALING_VIA_NFKB# HALLMARK_MYOGENESIS# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION# HALLMARK_ADIPOGENESIS# 
# HALLMARK_HYPOXIA# HALLMARK_APICAL_JUNCTION# HALLMARK_XENOBIOTIC_METABOLISM# HALLMARK_FATTY_ACID_METABOLISM# HALLMARK_KRAS_SIGNALING_UP# HALLMARK_UV_RESPONSE_UP
# ??? ASHP does not really understand why : because luminal B 
# Relevance? check with Ivan or Florence or Thierry ?

# Cluster new 3 (yellow) :  390 genes AURKA CCNE1 BIRC5 KIF23 PRC1 MYBL1 MYBL2
# Nothing very clear regarding samples 
# Split bien 2 sous groupes de luminales et 6  des TNBC (exclut 2 samples)
# HALLMARK_G2M_CHECKPOINT,HALLMARK_E2F_TARGETS,HALLMARK_MITOTIC_SPINDLE,HALLMARK_SPERMATOGENESIS,
# up in TNBC and luminal B

# Cluster new 4 (orange) : IDO1 , ROS1 
# Un peu a boire et a manger ; pas de cluster clair
# the biggest cluster 1107
# HALLMARK_KRAS_SIGNALING_DN,HALLMARK_KRAS_SIGNALING_UP,HALLMARK_ESTROGEN_RESPONSE_LATE,HALLMARK_INFLAMMATORY_RESPONSE,
# HALLMARK_COAGULATION,HALLMARK_COMPLEMENT,HALLMARK_HYPOXIA,HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,HALLMARK_MYOGENESIS,HALLMARK_XENOBIOTIC_METABOLISM

rlab <- as.data.frame(rownames(dataf2))
colnames(rlab) <- "gene"
rownames(rlab) <- rownames(dataf2) 
head(rlab)
rlab[cluster4_genes_new1,"gene"] <- '#219EBC'
rlab[cluster4_genes_new2,"gene"] <- '#023047'
rlab[cluster4_genes_new3,"gene"] <- '#FFB703'
rlab[cluster4_genes_new4,"gene"] <- '#FB8500'
head(rlab)

allcolcluster        <- c("#219EBC","#023047","#FFB703","#FB8500")
names(allcolcluster) <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
color_genes <- factor(rlab[,"gene"],     levels=((allcolcluster)))

rlab <- as.matrix(rlab) 
head(rlab)

data_clin_min[]         <- lapply(data_clin_min, as.factor)
data_clin_min$samplename

labSamples 	<- data_clin_min
clab 		    <- labSamples
clab

clab[,"patient"]      <- allcollevelpatient[match(clab[,"patient"],names(allcollevelpatient))]
clab[,"subtype"]      <- colsubtype_cvd [match(clab[,"subtype"]     ,names(colsubtype_cvd)) ]
clab[,"tumortype"]    <- colsubtype_cvd_pCR    [match(clab[,"tumortype"]   ,names(colsubtype_cvd_pCR   )) ]
clab[,"tumortype2"]   <- coltumortype2   [match(clab[,"tumortype2"] , names(coltumortype2  ))] # 
clab[,"side"]         <- collevelside2    [match(clab[,"side"]       , names(collevelside2   ))]
clab[,"BRCA_mut"]     <- collevelBRCA_3    [match(clab[,"BRCA_mut"]   , names(collevelBRCA_3   ))]
clab[,"BRCA1_2_mut"]  <- collevelBRCA_4    [match(clab[,"BRCA1_2_mut"]   , names(collevelBRCA_4   ))]


legendClab1 	<- c(  "subtype",names(collevelsubtype),"",  "BRCA status",names(collevelBRCA_4),"",
                    "side",names(collevelside2),"",  "sample type",names(coltumortype),"",
                    "patient",names(allcollevelpatient),"","","","",  "Gene cluster",names(allcolcluster))

color_patients <- factor(clab[,"patient"],     levels=((allcollevelpatient)))
color_brca12   <- factor(clab[,"BRCA1_2_mut"], levels=((collevelBRCA_4)))
color_side     <- factor(clab[,"side"],        levels=((collevelside2)))

legendClab2 	<- c(  "white", rev(levels(as.factor(clab[,"subtype"])    ))   ,"white", 
                    "white", levels(color_brca12)   ,      "white",  "white", levels(color_side)               ,"white", 
                    "white", levels(as.factor(clab[,"tumortype"]  ))        ,"white", 
                    "white", levels(color_patients)                                , "white",  "white",  "white",  "white",   "white", levels(color_genes)
                  )

colnames(rlab)[colnames(rlab) == "gene" ] <- "Gene \n cluster"

clab 			         <- as.matrix(clab[,c("patient","tumortype","side","BRCA1_2_mut","subtype")]) 
colnames(clab)[colnames(clab) == "BRCA1_2_mut" ]  <- "BRCA status"
colnames(clab)[colnames(clab) == "tumortype" ] <- "sample type"

pdf("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/Figures/Heatmap_RNASeq_BC_bilat_neo.pdf", height = 8)
                            main_title = ("")
                            par(cex.main=0.9)
                            heatmap.3(dataf, 
                                      Rowv = as.dendrogram(cl.row), 
                                      Colv = as.dendrogram(cl.col)			,	# ,
                                      ColSideColors	= clab,
                                      RowSideColors	= t(rlab),
                                      na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,11),
                                      symbreaks=FALSE, symkey=FALSE, 
                                      key=TRUE, KeyValueName="Gene Expression", keysize = 1.1, 
                                      density.info="none", trace="none", main=main_title, 
                                      labRow =  FALSE,
                                      ColSideColorsSize=6,
                                      RowSideColorsSize=1.3,
                                      col=my_palette)
                            legend("topright",legend	= c( legendClab1),
                                   fill	= c(  	legendClab2),
                                   border=FALSE, bty="n", y.intersp = 0.8, cex=0.6)
dev.off()

# The clustering first split samples into a group of luminal tumors and a group of TNBCs 
# (except for the ER weakly positive (+25%) luminal tumor PT5_L clustering into the TNBC group) (Fig6A), 
# while gene clustering split the 2846 genes into four main clusters. 


# The group of luminal tumors was globally enriched in genes from cluster 1 (early and late response to estrogens); 
# and a subset of luminal tumors were also enriched in genes from the cluster 2 (TNF signaling, myogenesis, epithelial mesenchymal transition). 
# The majority of TNBC samples were enriched in genes from cluster 3 (G2M checkpoints, E2F targets, cellular cycle). 
# Genes from cluster 4 showed no clear enrichment in specific pathways and were expressed in complex patterns by both luminal 
# and TNBC tumor subsets. Within each subgroup, the PT samples consistently clustered with their related RD rather than the tumor
# from the contralateral side. 

# Perform PCA analysis

# PCA data table
pca_data     <- tpm %>% select(-c(Length)) %>% column_to_rownames("Geneid") ; head(pca_data)
# Order samples 
pca_data     <- pca_data[,levelsample] 
pca_data %>%  rowMeans(.) %>% order(., decreasing=TRUE) -> select
sample_names <- names(pca_data)

# ashp_colors
# samples_recap_DNA_RNA_color <- read_excel("~/RT2Lab/bc_bilat/NGS/WES/preprocessing_WES_data/data/raw/Dataset_KDI/samples_recap_DNA_RNA_color_v3.xlsx")
# samples_recap_DNA_RNA_color <- samples_recap_DNA_RNA_color %>% filter(!is.na(num_histo))
# head(samples_recap_DNA_RNA_color)
# 
# ashp_colors_tmp    <- samples_recap_DNA_RNA_color %>% dplyr :: select(sample_name, color) %>% unique()
# ashp_colors        <- ashp_colors_tmp$color
# names(ashp_colors) <- ashp_colors_tmp$sample_name

pca_plots <- function(nb_genes=3000){
  #   nb_genes = 3000
  top.pca_data = pca_data[select[1:nb_genes],] %>% t()
  
  pca_prcomp=prcomp(top.pca_data) 
  names(pca_prcomp)
  
  ## Calculating the variance covered by PCs
  pca_data_perc=round(100*pca_prcomp$sdev^2/sum(pca_prcomp$sdev^2),1)
  
  df_pca_data = data.frame(PC1 = pca_prcomp$x[,1], PC2 = pca_prcomp$x[,2], 
                           sample = factor(sample_names, 
                                           levels = sample_names))
  
  df_pca_data <- df_pca_data %>% left_join(.,data_clin_min_new %>% rename(sample = samplename ) %>% select(sample, pCR_status) )
  
  p1 <- ggplot(df_pca_data,
               aes(PC1, PC2, color = sample, label = sample, shape = pCR_status))+
    geom_point(size=4)+
    scale_color_manual(name="Sample Name",
                       labels = names(ashp_colors),
                       values = ashp_colors) +
    labs(#title = "BILAT project -- PCA plot -- PC1 vs PC2",
      subtitle = paste0("top ", nb_genes, " variable genes"),
      x = paste("PC1, ", pca_data_perc[1], "%", sep=""), 
      y = paste("PC2, ", pca_data_perc[2], "%", sep=""))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line.x = element_line(colour = "grey30", size = 0.5),
          axis.line.y = element_line(colour = "grey30", size = 0.5),
          plot.margin = margin(10,10,10,10),
          legend.background = element_rect(fill = "white"))+
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5), #title
          plot.subtitle = element_text(size = 16, hjust = 0.5), #subtitle
          axis.title.x=element_text(size=18),  # X axis title
          axis.title.y=element_text(size=18),  # Y axis title
          axis.text.x=element_text(size=14, angle = 0, vjust=.5),  # X axis text
          axis.text.y=element_text(size=14), # Y axis text
          legend.title = element_text(size =18, face="bold"),
          legend.text = element_text(size =14),
          legend.position = "none")+
    geom_hline(yintercept = 0, color = "grey60", size = 0.2)+
    geom_vline(xintercept = 0, color = "grey60", size = 0.2)+
    ggrepel::geom_text_repel(size = 6, point.padding = 0.01)
  p1
  
  # df_pca_data = data.frame(PC2 = pca_prcomp$x[,1], PC3 = pca_prcomp$x[,3], # Previous line from mathias ; error ??
  df_pca_data = data.frame(PC2 = pca_prcomp$x[,2], PC3 = pca_prcomp$x[,3], 
                           sample = factor(sample_names, 
                                           levels = sample_names))
  
  df_pca_data <- df_pca_data %>% left_join(.,data_clin_min_new %>% rename(sample = samplename ) %>% select(sample, pCR_status) )
  
  p2 <- ggplot(df_pca_data,
               aes(PC2, PC3, color = sample, label = sample, shape = pCR_status))+
    geom_point(size=4)+
    scale_color_manual(name="Sample Name",
                       labels = names(ashp_colors),
                       values=ashp_colors) +
    labs(#title = "BILAT project -- PCA plot -- PC2 vs PC3",
      subtitle = paste0("top ", nb_genes, " variable genes"),
      x = paste("PC2, ", pca_data_perc[2], "%", sep=""), 
      y = paste("PC3, ", pca_data_perc[3], "%", sep=""))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line.x = element_line(colour = "grey30", size = 0.5),
          axis.line.y = element_line(colour = "grey30", size = 0.5),
          plot.margin = margin(10,10,10,10),
          legend.background = element_rect(fill = "white"))+
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5), #title
          plot.subtitle = element_text(size = 16, hjust = 0.5), #subtitle
          axis.title.x=element_text(size=18),  # X axis title
          axis.title.y=element_text(size=18),  # Y axis title
          axis.text.x=element_text(size=14, angle = 0, vjust=.5),  # X axis text
          axis.text.y=element_text(size=14), # Y axis text
          legend.title = element_text(size =18, face="bold"),
          legend.text = element_text(size =14),
          legend.position = "bottom")+
    scale_shape_manual(name = " ",
                       values = c(16,17),
                       labels= c("Pair PT/RD","Tumor reaching pCR")) +
    guides(color = FALSE) +
    guides(shape = guide_legend(nrow=2)) +
    geom_hline(yintercept = 0, color = "grey60", size = 0.2)+
    geom_vline(xintercept = 0, color = "grey60", size = 0.2)+
    ## adding label text, point.padding - protected distance around points
    ggrepel::geom_text_repel(size = 6, point.padding = 0.01)
  p2
  
  df_pca_data = data.frame(PC1 = pca_prcomp$x[,1], PC3 = pca_prcomp$x[,3], 
                           sample = factor(sample_names, 
                                           levels = sample_names))
  
  df_pca_data <- df_pca_data %>% left_join(.,data_clin_min_new %>% rename(sample = samplename ) %>% select(sample, pCR_status) )
  
  p3 <- ggplot(df_pca_data,
               aes(PC1, PC3, color = sample, label = sample, shape = pCR_status))+
    geom_point(size=4)+
    scale_color_manual(name="Sample Name",
                       labels = names(ashp_colors),
                       values=ashp_colors) +
    labs(#title = "BILAT project -- PCA plot -- PC1 vs PC3",
      subtitle = paste0("top ", nb_genes, " variable genes"),
      x = paste("PC1, ", pca_data_perc[1], "%", sep=""), 
      y = paste("PC3, ", pca_data_perc[3], "%", sep=""))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line.x = element_line(colour = "grey30", size = 0.5),
          axis.line.y = element_line(colour = "grey30", size = 0.5),
          plot.margin = margin(10,10,10,10),
          legend.background = element_rect(fill = "white"))+
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5), #title
          plot.subtitle = element_text(size = 16, hjust = 0.5), #subtitle
          axis.title.x=element_text(size=18),  # X axis title
          axis.title.y=element_text(size=18),  # Y axis title
          axis.text.x=element_text(size=14, angle = 0, vjust=.5),  # X axis text
          axis.text.y=element_text(size=14), # Y axis text
          legend.title = element_text(size =18, face="bold"),
          legend.text = element_text(size =14),
          legend.position = "none")+
    # scale_shape_manual(name=" ") +
    geom_hline(yintercept = 0, color = "grey60", size = 0.2)+
    geom_vline(xintercept = 0, color = "grey60", size = 0.2)+
    ggrepel::geom_text_repel(size = 6, point.padding = 0.01)
  p3
  
  fig <- ggpubr::ggarrange(p1,p3,p2, nrow = 3, ncol = 1)
}

save_plot(fig, file = "~/RT2Lab/BC_BILAT_NEO/codes_git/Figures/pca_3000_genes.pdf", 
          base_width = 9, base_height = 20)

# Similar results were seen after principal component analysis (PCA) using the 3000 most variable genes (Fig6B). 
# This suggests that PT and RD are closer from a transcriptomic point of view than are left and right tumors from the same patient.  
