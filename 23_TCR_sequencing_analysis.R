# T-cell receptor (TCR) sequencing analysis


################################################## 
#                       SETUP
################################################## 
rm(list = ls())

# Packages
library(immunarch)
library(ggplot2)
library(dplyr)
library(cowplot)

# Data and masterfiles
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer2/newmicr/bilat_immdata2.Rdata")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer2/newmicr/df_all_data.Rdata")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/data_clin_min_new.RData")
head(data_clin_min_new)

# Functions and colors
colsubtype_cvd     <- c("luminal" = "#0072B2","TNBC"    = "#D55E00","HER2+"   = "#009E73") # [1] "#F8766D" # rouge "#00BA38" vert "#619CFF" bleu

colsubtype_cvd_pCR <- c("luminal" = "#0072B2",
                        "TNBC"    = "#D55E00",
                        "HER2+"   = "#009E73",
                        "pCR"     = "#6ABBEB",
                        "No pCR"  = "#044970",
                        "different pCR status"  = "#F5AE78")

new_sample_order = c('PT1_L', 'RD1_L', 'PT1_R', 'PT2_L', 'PT2_R','PT3_L', 'RD3_L', 'PT3_R','RD3_R',
                     'PT4_L','PT4_R', 'RD4_R','PT5_L','PT5_R', 'RD5_R',
                     'PT6A_L','PT6B_L', 'PT6A_R','PT6B_R', 'RD6_R' )
new_sample_order

immdata <- list()

immdata$data <- immdata_data
immdata$data <- immdata_data[new_sample_order]
names(immdata$data)[names(immdata$data) == "RD6_R"] <- "RD6B_R" 
immdata$meta <- immdata_meta[match(new_sample_order,immdata_meta$Sample),]
immdata$meta[which(immdata$meta$Sample == "RD6_R"),"Sample"] <- "RD6B_R" 
immdata$meta$Patient <- immdata$meta$Sample
immdata$meta$Patient <- gsub("PT","",immdata$meta$Patient)
immdata$meta$Patient <- gsub("RD","",immdata$meta$Patient)
immdata$meta$Patient <- gsub("_L","",immdata$meta$Patient)
immdata$meta$Patient <- gsub("_R","",immdata$meta$Patient)
immdata$meta$Patient <- gsub("_A","",immdata$meta$Patient)
immdata$meta$Patient <- gsub("_B","",immdata$meta$Patient)
immdata$meta$Patient <- paste0("Patient ",immdata$meta$Patient)

.col = immunarch:::process_col_argument("aa")
.col

.data = immdata$data
for (i in 1:length(.data)) {
  if (immunarch:::has_class(.data[[i]], "data.table")) {
    .data[[i]] <- .data[[i]] %>% lazy_dt() %>% select(.col) %>% 
      collect(n = Inf)
  }
  else {
    .data[[i]] <- .data[[i]] %>% select(.col) %>% 
      collect(n = Inf)
  }
}


.data %>% unlist() %>% as.data.frame() %>% 
  rownames_to_column(var="sample_col") %>% 
  magrittr::set_colnames(c("sample_col", "seq"))       %>% 
  mutate(sample=stringr::str_split(sample_col , "\\.") %>% map_chr(., 1)) %>% 
  dplyr::select(sample, seq ) %>% unique() %>% dplyr::select(seq) %>% group_by(seq) %>% summarise(count=n())  %>%
  arrange(desc(count)) -> seq_counts

sample_clonotypes_p2 <- seq_counts  %>%  dplyr::select(count) %>% group_by(count) %>%
                  ggplot(aes(x=count)) +
                  geom_histogram(binwidth = 0.5, bins=30, color="darkblue", fill="lightblue") + scale_y_log10() +
                  theme_bw() +
                  scale_x_continuous(breaks=seq(1,20) ) +
                  xlab("Occurrence (#samples)")+
                  ylab("Number of Clonotypes (log10 scale)") +
                  ggtitle("Clonotype (CDR3 AA) Occurrence - Number of samples -- log10 scale") 

# The large majority of clonotypes retrieved were unique to a sample (Extended Fig10A) but some sequences were found in multiple samples. 

# Cf github
# cf https://github.com/antigenomics/vdjdb-db

# score 	description
# 0 	Low confidence/no information - a critical aspect of sequencing/specificity validation is missing
# 1 	Moderate confidence - no verification / poor TCR sequence confidence
# 2 	High confidence - has some specificity verification, good TCR sequence confidence
# 3 	Very high confidence - has extensive verification or structural data

# Cf former script_cdj_PG_BC_BILAT_NEO

vdjdb  <- read.csv("~/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/MiXCR/vdjdb/vdjdb-2018-06-04 2/vdjdb.txt", sep = "\t")
vdjdb2 <- read.csv("~/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/MiXCR/vdjdb/vdjdb-2018-06-04 2/vdjdb.slim.txt", sep = "\t")
vdjdb3 <- read.csv("~/RT2Lab/BC_BILAT_NEO/NGS/RNASeq/MiXCR/vdjdb/vdjdb-2018-06-04 2/vdjdb_full.txt", sep = "\t")
head(vdjdb)

df_all_data$antigen.species <- vdjdb[,c("antigen.species")][match(df_all_data$CDR3.aa,vdjdb[,c("cdr3")])    ]
table(df_all_data$antigen.species)
# CMV              EBV              HCV            HIV-1      HomoSapiens       InfluenzaA              SIV YellowFeverVirus 
# 8               10                2                3                3               10                1                2 

# In unique
seq_CDR3aa             <- df_all_data %>% group_by(CDR3.aa) %>% count() %>% filter(n>1) %>% 
                        select(CDR3.aa) %>% as.matrix() %>% as.character() # 235 sequences shared accross cohort
df_all_data$shared     <- ifelse(df_all_data$CDR3.aa %in% seq_CDR3aa, "Yes","No")
df_all_data$vdjdb_bin  <- ifelse(!is.na(df_all_data$antigen.species) , "Yes","No")
df_all_data %>% select(shared,vdjdb_bin) %>% table() %>% chisq.test() # P = 0.70

# The proportion of samples annotated in VDJdb, a curated database of T-cell receptor sequences of known antigen specificity71, 
# was low (1%) and was not different in sequences unique to a sample and in sequences shared accrss the cohort (8/638 versus 31/3126 respectively, p=0.7). 


# To further investigate the T cell response to NAC and to compare infiltrating TCR repertoires across patients, 
# we extracted TCR beta CDR3 sequences from the RNAseq data using MixCR31, and Immunarch32. 


load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer2/immunarch/immunarchplots/immunarch/div_chao_df.Rdata")
head(div_chao_df)
# chao 1
div_chao_df_merged <- div_chao_df %>% rownames_to_column() %>% 
                                      rename(samplename = rowname) %>% 
                                      left_join(., data_clin_min_new %>% select(samplename,subtype,tumortype,patient_side,tumortype2)) 


p1 <- vis(div_chao) + theme(axis.text.x  = element_text(size = 7))  +  scale_x_discrete(limits=c(new_sample_order))

## By PT / RD 
p_Chao_by_pt_rd <- div_chao_df_merged %>% 
                          ggplot(aes(x = tumortype , y = Estimator, fill = tumortype)) +
                          geom_boxplot(aes(x = tumortype , y = Estimator, fill = tumortype)) +
                          geom_jitter(aes(x = tumortype , y = Estimator, fill = tumortype))+
                          scale_fill_manual(values = c("#3182BD","#DEEBF7")) + theme_bw() +
                          theme(legend.position = "none", axis.ticks.x = element_blank() )+
                          stat_compare_means( label = "p.format", label.x.npc="center", label.y = 1000)+
                          xlab("")+ ylab("Chao1")
p_Chao_by_pt_rd

## By subtype
p_Chao_by_subtype <- div_chao_df_merged %>% 
                          ggplot(aes(x = subtype , y = Estimator, fill = subtype)) +
                          geom_boxplot(aes(x = subtype , y = Estimator, fill = subtype)) +
                          geom_jitter(aes(x = subtype , y = Estimator, fill = subtype)) +
                          scale_fill_manual(values = colsubtype_cvd_pCR) + theme_bw() +
                          theme(legend.position = "none", axis.ticks.x = element_blank() )+
                          stat_compare_means( label = "p.format", label.x.npc="center", label.y = 1000)+
                          xlab("")+ ylab("Chao1")
p_Chao_by_subtype

load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer2/immunarch/immunarchplots/immunarch/div_d50_df.Rdata")

p5     <- vis(div_d50) + theme(axis.text.x  = element_text(size = 7)) +  scale_x_discrete(limits=c(new_sample_order))

div_d50_df_merged <- div_d50_df %>% rownames_to_column() %>% rename(samplename = rowname) %>% left_join(., data_clin_min_new %>% select(samplename,subtype,tumortype,patient_side,tumortype2)) 

## By PT / RD 
p_D50_by_pt_rd <- div_d50_df_merged %>% ggplot(aes(x = tumortype , y = Clones, fill = tumortype)) +
                          geom_boxplot(aes(x = tumortype , y = Clones, fill = tumortype)) +
                          geom_jitter(aes(x = tumortype , y = Clones, fill = tumortype)) +
                          scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                          theme_bw() +
                          theme(legend.position = "none", axis.ticks.x = element_blank() )+
                          stat_compare_means( label = "p.format", label.x.npc="center",label.y = 120)+
                          xlab("")+ ylab("D50")
p_D50_by_pt_rd

## By subtype
p_D50_by_subtype <- div_d50_df_merged %>% 
                          ggplot(aes(x = subtype , y = Clones, fill = subtype)) +
                          geom_boxplot(aes(x = subtype , y = Clones, fill = subtype)) +
                          geom_jitter(aes(x = subtype , y = Clones, fill = subtype)) +
                          scale_fill_manual(values = colsubtype_cvd_pCR) + theme_bw() +
                          theme(legend.position = "none", axis.ticks.x = element_blank() )+
                          stat_compare_means( label = "p.format", label.x.npc="center",label.y = 120)+
                          xlab("")+ ylab("D50")
p_D50_by_subtype

p_compil_tcr_2 <- plot_grid(    sample_clonotypes_p2 + ggtitle("Clonotype (CDR3 AA) Occurrence", subtitle = "Number of samples -- log10 scale"), 
                                p1 + xlab("") + theme (),
                                plot_grid(p_Chao_by_pt_rd    + ggtitle("by tumor type"),
                                          p_Chao_by_subtype  + ggtitle("by BC subtype"), ncol=1),
                                p5+ xlab("") ,
                                plot_grid(p_D50_by_pt_rd     + ggtitle("by tumor type"),
                                          p_D50_by_subtype   + ggtitle("by BC subtype"), ncol=1),
                                labels = "AUTO",
                                rel_widths = c(8,7,3,7,3), nrow = 1)

save_plot(p_compil_tcr_2, file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig10_p_compil_tcr_2.pdf",base_height = 5, base_width = 22)

p_compil_tcr_3 <- plot_grid(    sample_clonotypes_p2 + ggtitle("Clonotype (CDR3 AA) Occurrence", subtitle = "Number of samples -- log10 scale"), 
                                sample_clonotypes_p2 + xlab("") + theme (),
                                plot_grid(p_Chao_by_pt_rd    + ggtitle("by tumor type"),
                                          p_Chao_by_subtype  + ggtitle("by BC subtype"), ncol=1),
                                sample_clonotypes_p2+ xlab("") ,
                                plot_grid(p_D50_by_pt_rd     + ggtitle("by tumor type"),
                                          p_D50_by_subtype   + ggtitle("by BC subtype"), ncol=1),
                                labels = "AUTO",
                                rel_widths = c(8,7,3,7,3), nrow = 1)

save_plot(p_compil_tcr_3, file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig10_p_compil_tcr_3.pdf",base_height = 5, base_width = 22)


# We evaluated the diversity of the TCR using Chao-1 estimator of species richness (Extended Fig10B-D), 
# and the D50 diversity index (Extended Fig10E-G), and they were not different by BC subtype nor PT or RD character of the sample. 

imm_ov1 = repOverlap(immdata$data, .method = "public", .verbose = F)

options(repr.plot.width=12, repr.plot.height=8)
p_rep_overlap_public <- immunarch:::add_class(imm_ov1[new_sample_order,new_sample_order], "immunr_ov_matrix") %>% vis(.text.size=4, .order = new_sample_order) + 
                        ggtitle("Repertoire overlap", subtitle = "Public")+theme(plot.title = element_text(size = 20, face = "bold")) +
                        theme(plot.subtitle = element_text(size = 15, face = "plain")) +theme(axis.text.x = element_text(size=12)) +
                        theme(axis.text.y = element_text(size=12)) +theme(axis.text = element_text(size=12))+
                        theme(axis.title = element_text(size = 20))
p_rep_overlap_public

# save_plot(p_rep_overlap_public, file = "/Users/paulgougis/Desktop/Code/other - general/immunarch ASHP/p_rep_overlap_public.pdf",base_height = 8, base_width = 8)
save_plot(p_rep_overlap_public, file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig10_p_rep_overlap_public.pdf",base_height = 8, base_width = 8)

# To measure repertoire similarity, we calculated the total number of shared clones between samples against “public” clonotypes (Extended Fig10H).
# We found shared TCRs between individuals at a low frequency, while most common sharing relationships were found between PT and RD 
# (yellow bordered squares), and to a smaller extent between left and right tumors (green bordered squares)
# though the median number of shared clonotypes was not statistically significant (20 versus 11, p=0.12). 

p_heatmap_overlap_repertoire <- vis(imm_ov1, "heatmap2") 

# save_plot(p_heatmap_overlap_repertoire, file = "/Users/paulgougis/Desktop/Code/other - general/immunarch ASHP/p_heatmap_overlap_repertoire.pdf",base_height = 8, base_width = 8)
save_plot(p_heatmap_overlap_repertoire, file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig10_p_heatmap_overlap_repertoire.pdf",base_height = 8, base_width = 8)

# Except for 2 samples which showed low sharing with any other sample (PT3_R, PT5_L), clonotypes of the same patients consistently clustered 
# together, either close from the contralateral side or from the corresponding RD/PT, consistent with a systemic effect of TILs (Extended Fig10I).
