################################################## 
#                       SETUP
################################################## 
rm(list = ls())

# Initial scripts
# Q19_R2_boxplots_pre_post_chimio.r
# Q2_R4_chemotherapy_NET_regimen_and_waterfall_plot_v4.r

# Packages
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggsci)

# Data and masterfiles
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData") # Raw matrix tumor data  (1 line per tumor)
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData") # Raw matrix patients data (1 line per patient) 
load("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/data/raw/dynamics_before_after.RData") 
load("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/data/raw/Dat.label_tils.RData")     
load("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/data/processed/dataf_melt_2.RData") # Data to make waterfall plots 
ls()

# Functions and colors
source('~/RT2Lab/BC_BILAT_NEO/NGS_new/WES_new/src/patients_tumors_levels_color.R')

col_nat_regimen <- c("Taxanes"             = "firebrick3",
                     "Anthracyclines"      = "darkorange",
                     "Anthra-taxanes"      = "deepskyblue3",
                     "Aromatase inhibitor" = "#1A936F",
                     "Tamoxifen"           = "#FFEE93")

col_first_ttt <- c(  "Chemotherapy"      = "deepskyblue3",
                     "Endocrine therapy" = "#1A936F")

colsubtype_cvd_pCR_2 <- c("luminal" = "#0072B2",
                          "TNBC"    = "#D55E00", 
                          "HER2+"   = "#009E73",
                          "pCR"     = "#6ABBEB",
                          "No pCR"  = "#044970",
                          "different pCR status"  = "#F5AE78",
                          "HER2+ Discordant"   = "#009E73",
                          "luminal Concordant" = "#0072B2",
                          "TNBC Concordant"    = "#D55E00",
                          "luminal Discordant" = "#0099ff",
                          "TNBC Discordant"    = "#ff9933") 

col_grade <- c("1" = "#fcbf49",
               "2" = "#f77f00",
               "3" = "#d62828")

# TILs variation before and after neoadjuvant treatment (NAT).

mat_patient_no_is_nac_net <- mat_patient_no_is %>% filter(first_ttt !="Surgery"  )
tumor_bilat_no_is_nac_net <- tumor_bilat_no_is %>% filter(first_ttt !="Surgery"  )

# A total of 70 patients received NAT with chemotherapy 

tumor_bilat_no_is_nac_net %>%   filter(first_ttt == "Endocrine therapy") %>% 
  group_by(NUMDOS,typht) %>% count() %>% as.data.frame() %>% group_by(typht,n) %>% count()

tumor_bilat_no_is_nac_net %>%   filter(first_ttt == "Chemotherapy") %>% 
  group_by(NUMDOS,NAT) %>% count() %>% as.data.frame() %>% group_by(NAT,n) %>% count()

# (n=50, 46 of whom received anthracyclines and taxanes based sequential regimen) 
# or neoadjuvant endocrine therapy (NET, n=20, 18 of whom received aromatase inhibitors).

p_str_tils_changes_before_after_dynamics_and_concordance  <-   ggplot(dynamics_before_after %>% filter(!is.na(concordance_subtype)) ) +
                                        aes	(x=new_NUMDOS,y=tils, color =str_TILs_dynamics_txt, shape = setting_pre_post) +
                                        theme_bw()+    theme(axis.ticks.y = element_blank(),legend.position = "bottom")+
                                        scale_color_manual(name = " ",values = c("#003366", "#ff0000") ) +
                                        scale_linetype_manual(name = " ",values = c("luminal" = "solid","TNBC"    = "dotted","HER2+" = "dashed"))+
                                        facet_grid( concordance_subtype ~ cote, scales = "free", space = "free",#)   +
                                                    labeller = labeller( cote = c("D" = "Right side", "G" = "Left side"),
                                                                         concordance_subtype = c("Concordant" = "Concordant pair",
                                                                                                 "Discordant" = "Discordant pair"))) +
                                        coord_flip()+  geom_line(aes	(x=new_NUMDOS,y=tils,linetype = subtype,group=new_NUMDOS), position = position_dodge(width=0.1),
                                                  arrow = arrow(length=unit(0.2,"cm") ,type = "closed",angle = 20) ,size = 0.4  )	+
                                        ylab("Stromal TIL levels variation") +    xlab(" ") +   guides(color=guide_legend(title = " "))+
                                        ggtitle("Evolution of stromal TIL levels before and after neoadjuvant treatment")

save_plot(p_str_tils_changes_before_after_dynamics_and_concordance,
          file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig3_p_str_tils_changes_before_after_dynamics_and_concordance.pdf",
          base_width = 7, base_height = 9)

# Paired pre and post NAT data on immune infiltration were available for 74 tumors (37 patients) and are displayed on Extended Fig3.

tumor_bilat_no_is %>% group_by(str_TILs_dynamics.f )  %>% filter(!is.na(str_TILs_dynamics.f)) %>% 
  summarise(count=n()) %>% mutate(percent_format = paste0(round(count/sum(count)*100,1), '%')) 
# str_TILs_dynamics.f count percent_format
# <chr>               <int> <chr>         
# 1 No TILS change         18 24.3%         
# 2 TILS decrease          30 40.5%         
# 3 TILS increase          26 35.1%         

# Stromal TIL levels decreased in 30 tumors (40.5%), remained stable in 18 (24.3%) and increased in 26 tumors (35.1%). 

#         Waterfall plots
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Concordance status
p_changes_TILs_abs_val_by_concordance <- ggplot(data = dataf_melt_2 %>%filter(!is.na(concordance_subtype)),
                                                aes(x=numdos7, y=str_TILS_changes_abs_val )    ) +
                                    coord_cartesian(ylim=c(-60, 30)) +
                                    geom_bar(aes(fill=factor(subtype_concordance_subtype)	), stat="identity",  width=0.7,position = position_dodge(width=0.4)) +
                                    theme_bw() +  theme( axis.text.x = element_blank() , axis.ticks.x = element_blank(),
                                           axis.title.y = element_text(angle=90, face="bold"),
                                           panel.spacing = unit(0, "lines"),
                                           plot.background = element_blank(),
                                           legend.title=element_blank(), legend.position = "bottom",
                                           strip.text = element_text(face = "bold"))   +
                                    xlab("")+ ylab("TILs variation (abs. val.)")+
                                    facet_grid(.~concordance_subtype, scales = "free",space="free") +
                                    scale_fill_manual(values =  colsubtype_cvd_pCR_2,
                                      breaks = c('luminal Concordant','TNBC Concordant','luminal Discordant','TNBC Discordant','HER2+ Discordant'), 
                                      labels=  c('luminal Concordant','TNBC Concordant                                                                        ','luminal Discordant','TNBC Discordant','HER2+ Discordant'))  +
                                    guides(fill = guide_legend(nrow=1,byrow=TRUE)) 
print(p_changes_TILs_abs_val_by_concordance)

p_changes_TILs_abs_val_by_concordance_boxplot <- ggplot(data = dataf_melt_2 %>% 
                                                          filter(!is.na(concordance_subtype)), 
                                                        aes(x=concordance_subtype, y=str_TILS_changes_abs_val )    ) +
                                      coord_cartesian(ylim=c(-60, 30)) +
                                      geom_boxplot(aes(fill=factor(concordance_subtype)	))+ 
                                      theme_bw() +theme( 
                                        axis.ticks.x = element_blank(),
                                        axis.title.y = element_text(angle=90, face="bold"),
                                        panel.spacing = unit(0, "lines"),
                                        plot.background = element_blank(),
                                        legend.title=element_blank(), legend.position = "none",
                                        strip.text = element_text(face = "bold"))   + 
                                      xlab("")+ ylab("TILs variation (abs. val.)")+
                                      scale_x_discrete(labels= c("Conc. ","Disc.")) +
                                      stat_compare_means(label = "p.format", label.x.npc="center")
                                      # stat_compare_means(label = "p.signif", label.x.npc="center")
print(p_changes_TILs_abs_val_by_concordance_boxplot)

# By tumor grade 
p_changes_TILs_by_grade <- ggplot(data = dataf_melt_2 %>% filter(!is.na(gradeclasse)), aes(x=numdos7, y=str_TILS_changes_abs_val )    ) +
                              coord_cartesian(ylim=c(-60, 30)) +
                              geom_bar(aes(fill=factor(gradeclasse)	), stat="identity",  width=0.7,
                                       position = position_dodge(width=0.4)) + theme_bw() +
                              theme( axis.text.x = element_blank() , axis.ticks.x = element_blank(),
                                     axis.title.y = element_text(angle=90, face="bold"),
                                     panel.spacing = unit(0, "lines"),
                                     plot.background = element_blank(),
                                     legend.title=element_blank(), legend.position = "bottom",
                                     strip.text = element_text(face = "bold"))   + # , legend.position = "bottom"
                              xlab("")+ ylab("TILs variation (abs. val.)")+ scale_x_discrete(expand = c(0,0))  +
                              facet_grid(.~gradeclasse, scales = "free",space="free",
                                         labeller = labeller(  gradeclasse = c(  "1" = "Grade 1","2" = "Grade 2","3" = "Grade 3"))) +
                              scale_fill_manual(name = "Tumor grade",values = col_grade )
print(p_changes_TILs_by_grade)

# Boxplot
p_changes_TILs_by_grade_boxplot <- ggplot(data = dataf_melt_2 %>% 
                                            filter(!is.na(gradeclasse)), 
                                          aes(x=gradeclasse, y=str_TILS_changes_abs_val )    ) +
                                            coord_cartesian(ylim=c(-60, 30)) +
                                            geom_boxplot(aes(fill=factor(gradeclasse)	))+ 
                                            theme_bw() +theme( 
                                              axis.ticks.x = element_blank(),
                                              axis.title.y = element_text(angle=90, face="bold"),
                                              panel.spacing = unit(0, "lines"),
                                              plot.background = element_blank(),
                                              legend.title=element_blank(), legend.position = "none",
                                              strip.text = element_text(face = "bold"))   + 
                                            xlab("")+ ylab("TILs variation (abs. val.)")+
                                            # stat_compare_means(label = "p.signif", label.x.npc="center")+
                                            stat_compare_means(label = "p.format", label.x.npc="center")+
                                            scale_fill_manual(name = "Tumor grade",values = col_grade )
print(p_changes_TILs_by_grade_boxplot)

# 3.  Waterfall plot by pre-NAC stromal TIL levels
# dataf_sorted <- dataf_melt[order(dataf_melt[,"str_TILS_10_perc"], 
#                                  dataf_melt[,"str_TILS_changes_abs_val"] )  , ]
dataf_sorted <- dataf_melt_2[order(dataf_melt_2[,"str_TILS_10_perc"], 
                                   dataf_melt_2[,"str_TILS_changes_abs_val"] )  , ]
dataf_sorted  <- dataf_sorted[which(!is.na(dataf_sorted$str_TILS_changes_abs_val)),]
dataf_sorted$tmp4  <- c(1:nrow(dataf_sorted) )
dataf_sorted$all <- "all"

p_sBBCs_TILs_sorted <- ggplot(data = dataf_sorted %>% filter(!is.na(str_TILS_10_perc)), aes(x=tmp4,y=str_TILS_changes_abs_val )    ) +
                        geom_bar(aes(fill=factor(str_TILS_10_perc)	), stat="identity",  width=0.7,position = position_dodge(width=0.4)) +
                        coord_cartesian(ylim=c(-60, 30)) +  theme_bw()+  theme( axis.text.x = element_blank() , axis.ticks.x = element_blank(),
                        axis.title.y = element_text(angle=90),legend.position = "bottom")  +
                        xlab("")+ ylab("TILs variation (abs. val.)")+
                        facet_grid(.~all, scales = "free",space="free", 
                                          labeller = labeller(all = c("all"="Pre-NAC TILs levels (by 10% increment)"))) +
                        scale_fill_jco(name = "Pre-NAC str TIL levels")+ guides(fill = guide_legend(nrow=1,byrow=TRUE))

# and boxplot
p_sBBCs_TILs_sorted_boxplot <- ggplot(data = dataf_sorted %>% 
                                        filter(!is.na(str_TILS_10_perc)), 
                                      aes(x=str_TILS_10_perc, y=str_TILS_changes_abs_val )    ) +
                                      coord_cartesian(ylim=c(-60, 30)) +  geom_boxplot(aes(fill=factor(str_TILS_10_perc)	))+ 
                                      # geom_boxplot()+ 
                                      theme_bw() +  theme( axis.ticks.x = element_blank(),axis.title.y = element_text(angle=90, face="bold"),
                                        panel.spacing = unit(0, "lines"),plot.background = element_blank(),
                                        legend.title=element_blank(), legend.position = "none",    strip.text = element_text(face = "bold"))   + 
                                      xlab("")+ ylab("TILs variation (abs. val.)")+
                                      stat_compare_means(label = "p.format", label.x.npc="center")+  scale_fill_jco(name = "Pre-NAC str TIL levels")
                                      # stat_compare_means(label = "p.signif", label.x.npc="center")+  scale_fill_jco(name = "Pre-NAC str TIL levels")
print(p_sBBCs_TILs_sorted_boxplot)

# 4  Waterfall plot by type of neoadjuvant treatment
p_changes_TILs_abs_val_by_NAT_regimen <- ggplot(data = dataf_melt_2, aes(x=numdos7, y=str_TILS_changes_abs_val )    ) +
                                          coord_cartesian(ylim=c(-60, 30)) +
                                          geom_bar(aes(fill=factor(NAT)	), stat="identity",  width=0.7,
                                                   position = position_dodge(width=0.4)) +theme_bw() +
                                          theme( axis.text.x = element_blank() , axis.ticks.x = element_blank(),
                                                 axis.title.y = element_text(angle=90, face="bold"),
                                                 panel.spacing = unit(0, "lines"),
                                                 plot.background = element_blank(),
                                                 legend.title=element_blank(), legend.position = "bottom",
                                                 strip.text = element_text(face = "bold"))   + # , legend.position = "bottom"
                                          xlab("")+ ylab("TILs variation (abs. val.)")+scale_x_discrete(expand = c(0,0))  +
                                          facet_grid(.~first_ttt, scales = "free",space="free") +
                                          scale_fill_manual(name = " ",values = col_nat_regimen ,
                                                            breaks=c('Anthra-taxanes', 'Anthracyclines', 'Taxanes',
                                                                     'Aromatase inhibitor',"Tamoxifen")) 
print(p_changes_TILs_abs_val_by_NAT_regimen)

p_changes_TILs_abs_val_by_NAT_regimen_boxplot <- ggplot(data = dataf_sorted %>% filter(!is.na(first_ttt)), 
                                                        aes(x=first_ttt, y=str_TILS_changes_abs_val )    ) +
                                                coord_cartesian(ylim=c(-60, 30)) +  geom_boxplot(aes(fill=factor(first_ttt)	))+ 
                                                theme_bw() +  theme( axis.ticks.x = element_blank(),
                                                  axis.title.y = element_text(angle=90, face="bold"),
                                                  panel.spacing = unit(0, "lines"),plot.background = element_blank(),
                                                  legend.title=element_blank(), legend.position = "none",
                                                  strip.text = element_text(face = "bold"))   + 
                                                xlab("")+ ylab("TILs variation (abs. val.)")+
                                                scale_x_discrete(labels= c("NAC ","NET")) +
                                                # stat_compare_means( label = "p.signif", label.x.npc="center")+
                                                stat_compare_means( label = "p.format", label.x.npc="center")+
                                                scale_fill_manual(name = " ", values = col_first_ttt )

# 5  Waterfall plot by pCR status
p_changes_TILs_by_pCR <- ggplot(data = dataf_melt_2, aes(x=numdos7, y=str_TILS_changes_abs_val )    ) +
                                        coord_cartesian(ylim=c(-60, 30)) +
                                        geom_bar(aes(fill=factor(pCR)	), stat="identity",  width=0.7,
                                                 position = position_dodge(width=0.4)) +
                                        theme_bw() +
                                        theme( axis.text.x = element_blank() , axis.ticks.x = element_blank(),
                                               axis.title.y = element_text(angle=90, face="bold"),
                                               panel.spacing = unit(0, "lines"),
                                               plot.background = element_blank(),
                                               legend.title=element_blank(), legend.position = "bottom",
                                               strip.text = element_text(face = "bold"))   + # , legend.position = "bottom"
                                        xlab("")+ ylab("TILs variation (abs. val.)")+
                                        scale_x_discrete(expand = c(0,0))  +
                                        facet_grid(.~pCR, scales = "free",space="free") +
                                        scale_fill_manual(name = " ",values = colsubtype_cvd_pCR_2, breaks = c("pCR", "No pCR") ) 
print(p_changes_TILs_by_pCR)

p_changes_TILs_by_pCR_boxplot <- ggplot(data = dataf_sorted %>% filter(!is.na(pCR)), aes(x=pCR, y=str_TILS_changes_abs_val )    ) +
                                  coord_cartesian(ylim=c(-60, 30)) +  geom_boxplot(aes(fill=factor(pCR)	))+ 
                                  theme_bw() +  theme( 
                                    axis.ticks.x = element_blank(),
                                    axis.title.y = element_text(angle=90, face="bold"),
                                    panel.spacing = unit(0, "lines"),
                                    plot.background = element_blank(),
                                    legend.title=element_blank(), legend.position = "none",
                                    strip.text = element_text(face = "bold"))   + 
                                  xlab("")+ ylab("TILs variation (abs. val.)")+
                                  # stat_compare_means(  label = "p.signif", label.x.npc="center")+
                                  stat_compare_means(  label = "p.format", label.x.npc="center")+
                                  scale_fill_manual(name = " ",values = colsubtype_cvd_pCR_2, breaks = c("pCR", "No pCR") ) 

p_all_waterfall_sBBCs <- cowplot::plot_grid(p_changes_TILs_abs_val_by_concordance, 
                                            p_changes_TILs_by_grade,
                                            p_sBBCs_TILs_sorted,
                                            p_changes_TILs_abs_val_by_NAT_regimen,
                                            p_changes_TILs_by_pCR,
                                            nrow = 5, 
                                            # rel_widths = c(3, 2.2), 
                                            labels = "AUTO",
                                            align ="hv")
p_all_waterfall_sBBCs

p_all_boxplot_sBBCs <- cowplot::plot_grid(    p_changes_TILs_abs_val_by_concordance_boxplot , 
                                              p_changes_TILs_by_grade_boxplot , 
                                              p_sBBCs_TILs_sorted_boxplot , 
                                              p_changes_TILs_abs_val_by_NAT_regimen_boxplot , 
                                              p_changes_TILs_by_pCR_boxplot , 
                                              ncol = 5, 
                                              rel_widths = c(2, 3,5,2,2),
                                              labels = c("F","G","H","I","J"),
                                              align ="hv")
p_all_boxplot_sBBCs

p_all_waterfall_sBBCs_with_all_boxplot_sBBCs <- cowplot::plot_grid(p_all_waterfall_sBBCs, p_all_boxplot_sBBCs,
                                                                   nrow = 2, rel_heights = c(15,3), align ="hv")

save_plot(p_all_waterfall_sBBCs_with_all_boxplot_sBBCs, 
          file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig4_p_all_waterfall_sBBCs_with_all_boxplot_sBBCs.pdf",
                    base_width = 11, base_height = 20)

# The decrease of TIL levels was of larger magnitude in tumors belonging to discordant pairs, higher tumor grade, 
# with high pre-NAC stromal TIL levels, and in case of treatment with NAC rather than NET, 
# and the TILs decrease was very strongly associated with the occurrence of a pCR (Extended Fig4). 

head(Dat.label_tils)

# As a whole
p_str_tils_before_after_sBBCs <- Dat.label_tils %>% 
                                    ggplot(aes(x = variable , y = value, fill = variable)) +
                                    geom_boxplot(aes(x = variable , y = value, fill = variable)) +
                                    scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                                    facet_grid(. ~ all) + theme_bw() +
                                    theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                    scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
                                    stat_compare_means( label = "p.format", label.x.npc="center")+
                                    # , method = "wilcox.test") +
                                    xlab("")+ ylab("str TIL levels (%)")

# By NAT
p_str_tils_before_after_sBBCs_by_NAT <- Dat.label_tils %>% 
                                    ggplot(aes(x = variable , y = value, fill = variable)) +
                                    geom_boxplot(aes(x = variable , y = value, fill = variable)) +
                                    scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                                    facet_grid(. ~ first_ttt) +   theme_bw() +
                                    theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                    scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
                                    stat_compare_means( label = "p.format", label.x.npc="center", method = "wilcox.test") +
                                    xlab("")+ ylab("str TIL levels (%)")

# By subtype
p_str_tils_before_after_sBBCs_by_subtype <- Dat.label_tils %>% filter(!is.na(subtype)) %>%
                                    ggplot(aes(x = variable , y = value, fill = variable)) +
                                    geom_boxplot(aes(x = variable , y = value, fill = variable)) +
                                    scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                                    facet_grid(. ~ subtype) + 
                                    theme_bw() +  theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                    scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
                                    stat_compare_means( label = "p.format", label.x.npc="center", method = "wilcox.test") +
                                    xlab("")+ ylab("str TIL levels (%)")

# By concordance
p_str_tils_before_after_sBBCs_by_concordance <- Dat.label_tils %>% filter(!is.na(concordance_subtype)) %>%
                                  ggplot(aes(x = variable , y = value, fill = variable)) +
                                  geom_boxplot(aes(x = variable , y = value, fill = variable)) +
                                  scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                                  facet_grid(. ~ concordance_subtype) + 
                                  theme_bw() +
                                  theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                  scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
                                  stat_compare_means( label = "p.format", label.x.npc="center", method = "wilcox.test") +
                                  xlab("")+ ylab("str TIL levels (%)")

# By grade
p_str_tils_before_after_sBBCs_by_grade <- Dat.label_tils %>% filter(!is.na(gradeclasse)) %>%
                                          ggplot(aes(x = variable , y = value, fill = variable)) +
                                          geom_boxplot(aes(x = variable , y = value, fill = variable)) +
                                          scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                                          facet_grid(. ~ gradeclasse, labeller = labeller(  gradeclasse = c(  "1" = "Grade I",
                                                                                                              "2" = "Grade II",
                                                                                                              "3" = "Grade III"))) + 
                                          theme_bw() +theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                          scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
                                          stat_compare_means( label = "p.format", label.x.npc="center", method = "wilcox.test") +
                                          xlab("")+ ylab("str TIL levels (%)")

# By pCR
p_str_tils_before_after_sBBCs_by_pCR <- Dat.label_tils %>% filter(!is.na(pCR)) %>%
                                        ggplot(aes(x = variable , y = value, fill = variable)) +
                                        geom_boxplot(aes(x = variable , y = value, fill = variable)) +
                                        scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
                                        facet_grid(. ~ pCR) + 
                                        theme_bw() +
                                        theme(legend.position = "none", axis.ticks.x = element_blank() )+
                                        scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
                                        stat_compare_means( label = "p.format", label.x.npc="center", method = "wilcox.test") +
                                        xlab("")+ ylab("str TIL levels (%)")
p_str_tils_before_after_sBBCs_by_pCR

# By str_TILS_10_perc
p_str_tils_before_after_sBBCs_by_pre_NAC_TILs <- Dat.label_tils %>% filter(!is.na(str_TILS_10_perc)) %>%
  ggplot(aes(x = variable , y = value, fill = variable)) +
  geom_boxplot(aes(x = variable , y = value, fill = variable)) +
  scale_fill_manual(values = c("#3182BD","#DEEBF7")) +
  facet_grid(. ~ str_TILS_10_perc) + 
  theme_bw() +
  theme(legend.position = "none", axis.ticks.x = element_blank() )+
  scale_x_discrete(labels= c("Pre-treatment","Post-treatment")) +
  stat_compare_means( label = "p.format", label.x.npc="center", method = "wilcox.test") +
  xlab("")+ ylab("str TIL levels (%)")
p_str_tils_before_after_sBBCs_by_pre_NAC_TILs

p_pre_post_Tils_compil_1st_row <- plot_grid(p_str_tils_before_after_sBBCs + scale_x_discrete(labels= c("Pre \n treatment","Post \n treatment")),
                                            p_str_tils_before_after_sBBCs_by_NAT + scale_x_discrete(labels= c("Pre \n treatment","Post \n treatment")),
                                            p_str_tils_before_after_sBBCs_by_concordance + scale_x_discrete(labels= c("Pre \n treatment","Post \n treatment")),
                                            rel_widths=c(2,3,3),nrow = 1,labels = "AUTO",align="hv")  

p_pre_post_Tils_compil_2nd_row  <- plot_grid(p_str_tils_before_after_sBBCs_by_grade + scale_x_discrete(labels= c("Pre \n treatment","Post \n treatment")),
                                             p_str_tils_before_after_sBBCs_by_pCR+ scale_x_discrete(labels= c("Pre \n treatment","Post \n treatment")),
                                             labels = c("C","D","E"),rel_widths=c(3,2),align="hv")  

p_pre_post_Tils_compil  <-  plot_grid(  p_pre_post_Tils_compil_2nd_row,  p_pre_post_Tils_compil_1st_row,  nrow = 2, align="hv")  

save_plot(p_pre_post_Tils_compil, 
          file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/FigS3_p_pre_post_Tils_compil.pdf", 
          base_width = 10, base_height = 9)

# As a whole, stromal TILs levels were not significantly different before and after NAT, 
# but pre and post-NAT stromal TILs levels were significantly different according to the type of NAT,
# in discordant, grade 3 tumors, and in tumors that reached pCR (FigS3). 
# These findings suggest that neoadjuvant treatment significantly reshapes the immune contexture of sBBCs.
