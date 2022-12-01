################################################## 
#                       SETUP
################################################## 

# Packages

# Data and masterfiles

load("~/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer2/neoantigens/all_HLA_id_filtered.RData")
all_HLA_id <- all_HLA_id_filtered 

# Functions and colors
source('~/RT2Lab/BC_BILAT_NEO/NGS_new/WES_new/src/patients_tumors_levels_color.R', local = TRUE)

colsubtype_cvd_pCR <- c("luminal" = "#0072B2",
                        "TNBC"    = "#D55E00",
                        "HER2+"   = "#009E73",
                        "pCR"     = "#6ABBEB",
                        "No pCR"  = "#044970",
                        "different pCR status"  = "#F5AE78")

shared_peptides             <- all_HLA_id %>% group_by(HLA_Peptide,patient) %>% count() %>% filter(n>1) %>% ungroup() %>% 
                            select(HLA_Peptide) %>% as.matrix() %>% as.character() 
all_HLA_id$shared_peptides  <- ifelse(all_HLA_id$HLA_Peptide %in% shared_peptides, "yes", "no")

Dat.label <-  all_HLA_id %>% 
  group_by( HLA_type, HLA,patient ) %>% 
  dplyr::summarise(count=n()) %>%  
  group_by(HLA) %>%
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  mutate(n_and_perc   = paste0("n=", count,", (",percent_format,")" ) )   %>%
  ungroup() %>% arrange(HLA_type,count)

Dat.label$HLA <- as.character(Dat.label$HLA)
Dat.label$HLA <- factor(Dat.label$HLA, levels=unique(Dat.label$HLA))

p_HLA_type_pt <- Dat.label %>% ggplot(.) +
                  geom_bar(aes(x=HLA,y=count, fill=patient), stat="identity") +
                  facet_grid(HLA_type ~ ., scales = "free", space = "free") +
                  coord_flip()+
                  scale_fill_manual(name = "", values = collevelpatient)+
                  theme_bw()+
                  ylab("") + xlab("")+     guides(fill = guide_legend(nrow=1,byrow=TRUE))+
                  ggtitle("Number of neoantigen predicted")
p_HLA_type_pt

# Neoantigens 
# We predicted potential neoantigens from somatic mutations using netMHCpan after determining HLA haplotypes with Seq2HLA.


# Most of the antigens were predicted from HLA-C, and the repartition of predicted neoantigens was evenly distributed across patients (FigS5 A).

df_str_tils_neoantigens <- all_HLA_id %>% select(-c(Line:BindLevel)) %>%   unique() %>%
                                          arrange(HLA_Peptide) %>%  group_by(patient, samplename, str_tils) %>% count()

p_str_tils_versus_neoantigens <- ggscatter(df_str_tils_neoantigens, x = "str_tils", y = "n",
                                           add = "reg.line",  shape=21, size = 4,
                                           add.params = list(color = "#81A0D4", fill = "#81A0D4"), conf.int = TRUE) + 
                                          stat_cor(method = "pearson") +   scale_fill_manual(name = "", values = allcollevelpatient)+
                                          theme_bw() + xlab("Str TIL levels") + ylab("Number of neoantigens")
p_str_tils_versus_neoantigens

# The number of neoantigens was positively correlated with the levels of stromal TILs (FigS5B)

df_count_subtype <- all_HLA_id %>% group_by(samplename,subtype) %>% count()

p_neoantigen_by_subtype <- df_count_subtype %>%   ggplot(aes(x = subtype , y = n, fill = subtype)) +
                                                  geom_boxplot(aes(x = subtype , y = n, fill = subtype)) +
                                                  geom_jitter(size=0.5) +
                                                  scale_fill_manual(values = colsubtype_cvd_pCR) +  theme_bw() +
                                                  theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                                  stat_compare_means( label = "p.format", label.x.npc="center")+  xlab("")+ ylab("Number of neoantigens")

# was not associated with the BC subtype (FigS5C)

df_count_pt_rd <- all_HLA_id %>% group_by(samplename,tumortype) %>% count()

p_neoantigen_by_pt_rd <- df_count_pt_rd %>%   ggplot(aes(x = tumortype , y = n, fill = tumortype)) +
                                                  geom_boxplot(aes(x = tumortype , y = n, fill = tumortype)) +
                                                  geom_jitter(size=0.5) +
                                                  scale_fill_manual(values = c("#3182BD","#DEEBF7")) +  theme_bw() +
                                                  theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                                  stat_compare_means( label = "p.format", label.x.npc="center")+  xlab("")+ ylab("Number of neoantigens")

# and was higher in PT samples than in RD samples (FigS5D). No neoantigen was shared across patients, and no neoantigen was shared between 
# the left and the right tumors. Conversely, the RD shared most of the mutations with the corresponding PT (FigS5E, see script elsewhere).

p_boxplot_scatter_plot_neoantigen_1st_row <- plot_grid(p_str_tils_versus_neoantigens,
                                                       p_neoantigen_by_subtype,
                                                       p_neoantigen_by_pt_rd,
                                                       labels = c("","C","D"),
                                                       rel_widths=c(5,3.5,3.5),
                                                       nrow=1,
                                                       align="hv")
# not used in the plot
p_HLA_type_sample <- ggplot(data = all_HLA_id)+ #geom_histogram(stat="count") +
  facet_grid(. ~ patient, scales = "free", space = "free") +
  geom_bar(aes(x=samplename,  fill = shared_peptides)) +
  theme(axis.ticks.x = element_blank())+
  theme_bw()+
  ggtitle(" ") + ylab("")
p_HLA_type_sample

p_right_neoantigen <- plot_grid(p_boxplot_scatter_plot_neoantigen_1st_row,
                                p_HLA_type_sample + theme(legend.position = "bottom"),
                                nrow=2,
                                labels = c("","E"),
                                rel_heights=c(5,6),
                                align="hv")

p_compil_neoantigen <- plot_grid(p_HLA_type_pt + theme(legend.position = "bottom"),
                                 p_right_neoantigen,
                                 ncol=2,
                                 labels = c("A","B"),
                                 rel_widths=c(3,6),
                                 align="hv")

save_plot(p_compil_neoantigen,
          file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/FigS5_p_compil_neoantigen.pdf",
          base_width = 11.25, base_height = 7.7)

