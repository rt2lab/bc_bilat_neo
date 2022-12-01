################################################## 
#                       SETUP
################################################## 

rm(list = ls())

# Packages

library(cowplot)
library(dplyr)
library(ggplot2)
library(conflicted)
library(ggpubr)
conflict_prefer("filter","dplyr")
conflict_prefer("select","dplyr")

# Data and masterfiles
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData") # 634
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData") # 317
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_NAC_no_NET.RData")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer1/seer_bilat_neo.RData")  
load(file = "/Users/ahamypet/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer1/seer_bilat.RData") 

# Functions and colors
colsubtype_cvd     <- c("luminal" = "#0072B2",
                        "Luminal"="#0072B2",
                        "TNBC"    = "#D55E00", 
                        "HER2+"   = "#009E73")

# Figures concordance regarding BC subtype
## In the Curie cohort (Fig1A)

Dat.label_pat <- mat_patient_no_is %>% 
                  filter(!is.na(concordance_subtype)) %>% arrange(concordance_subtype) %>% 
                  group_by(couple_subtype) %>% 
                  summarise(count=n()) %>%
                  mutate(ypos = cumsum(count) - 0.5*count) %>%
                  mutate(percent_full = count/sum(count)) %>%
                  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>%
                  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full) %>% arrange(count) 

Dat.label_pat$couple_subtype <- factor(Dat.label_pat$couple_subtype, levels=unique(Dat.label_pat$couple_subtype))

Dat.label <- tumor_bilat_no_is %>% 
                filter(!is.na(concordance_subtype),
                       !is.na(couple_subtype)) %>% arrange(concordance_subtype) %>% 
                group_by(concordance_subtype,couple_subtype,subtype) %>% 
                summarise(count=n()) %>%
                mutate(ypos = cumsum(count) - 0.5*count) %>%
                mutate(percent_full = count/sum(count)) %>%
                mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>%
                mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full) %>% 
                arrange(concordance_subtype,count) %>% as.data.frame()

Dat.label$new_count      <- Dat.label_pat[match(Dat.label$couple_subtype,as.character(Dat.label_pat$couple_subtype)),"count"] %>% as.matrix() %>% as.character()
Dat.label$new_pos        <- ifelse(Dat.label$concordance_subtype == "Concordant",Dat.label$count,Dat.label$count*2)
Dat.label$couple_subtype <- factor(Dat.label$couple_subtype, levels=unique(Dat.label$couple_subtype))
values_patients          <- unique(Dat.label$new_count) %>% as.matrix() %>% as.character()

p_bar_curie  <- ggplot(Dat.label) + geom_bar(aes(x=couple_subtype,y=count,fill=subtype), stat="identity") +
                    facet_grid(concordance_subtype~., scales = "free_y", space = "free",
                               labeller = labeller(concordance_subtype = c("Concordant" = "Concordant pairs  (n=256, 85%) \n n=512 tumors", 
                                                                           "Discordant" = "Discordant pairs (n=46, 15%) \n n=92 tumors")) ) +
                    theme(strip.text.y.right = element_text(angle = 180))+
                    geom_text(aes(x=couple_subtype,y=new_pos, label = new_count ), hjust = -0.5, size = 4) +
                    theme_classic()+theme_bw()+
                    coord_flip(ylim = c(0,550)) + 
                    theme(legend.position = "none", axis.ticks = element_blank(),axis.line = element_blank(),axis.text.x = element_blank() ) +
                    xlab("")+ylab("")+ scale_fill_manual(values = c(colsubtype_cvd))+ labs(title ="Concordance between 2 sBBCs")

## In the SEER cohort (Fig1B)

Dat.label_tmp <- seer_bilat %>% 
                  filter(!is.na(concordance_subtype),
                         !is.na(couple_subtype)) %>% arrange(concordance_subtype) %>% #nrow() 
                  group_by(concordance_subtype,couple_subtype,subtype) %>% 
                  summarise(count=n()) %>%
                  arrange(concordance_subtype,count) %>% as.data.frame()

Dat.label      <- Dat.label_tmp  
Dat.label[2,]  <- Dat.label_tmp[1,]
Dat.label[1,]  <- Dat.label_tmp[2,]

Dat.label <- Dat.label %>%   mutate(ypos = cumsum(count) - 0.5*count) %>%
                                mutate(percent_full = count/sum(count)) %>%
                                mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>%
                                mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full) 

Dat.label_pat_tmp <- seer_bilat %>% select(Patient_register, concordance_subtype, couple_subtype) %>% unique() %>% 
  group_by(concordance_subtype,couple_subtype) %>%
  summarise(count=n()) 

Dat.label_pat      <- Dat.label_pat_tmp  
Dat.label_pat[2,]  <- Dat.label_pat_tmp[3,]
Dat.label_pat[3,]  <- Dat.label_pat_tmp[2,]

Dat.label_pat <- Dat.label_pat %>%
                    mutate(ypos = cumsum(count) - 0.5*count) %>%
                    mutate(percent_full = count/sum(count)) %>%
                    mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>%
                    mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  

Dat.label$new_count <- Dat.label_pat[match(Dat.label$couple_subtype,as.character(Dat.label_pat$couple_subtype)),"count"] %>% as.matrix() %>% as.character()
Dat.label$new_pos   <- ifelse(Dat.label$concordance_subtype == "Concordant",Dat.label$count,Dat.label$count*2)
Dat.label[which(Dat.label$new_pos == 25032),"new_pos"] <- 13000 
Dat.label$couple_subtype <- factor(Dat.label$couple_subtype, levels=unique(Dat.label$couple_subtype))
values_patients          <- unique(Dat.label$new_count) %>% as.matrix() %>% as.character()

p_bar_seer  <- ggplot(Dat.label) + 
                    geom_bar(aes(x=couple_subtype,y=count,fill=subtype), stat="identity") +
                    facet_grid(concordance_subtype~.,scales = "free_y", space = "free",
                               labeller = labeller(concordance_subtype = c("concordant" = "Concordant pairs  (n=6833, 82%) \n n=13666 tumors", 
                                                                           "discordant" = "Discordant pairs (n=1534, 18%) \n n=3068 tumors")) ) +
                    theme(strip.text.y.right = element_text(angle = 180))+
                    geom_text(aes(x=couple_subtype,y=new_pos, label = new_count ),hjust = -0.5, size = 4) +
                    theme_classic()+theme_bw()+coord_flip(ylim = c(0,15000)) +
                    theme(legend.position = "none", axis.ticks = element_blank(),axis.line = element_blank(),
                          axis.text.x = element_blank() ) +
                    xlab("")+ylab("")+ scale_fill_manual(values = c(colsubtype_cvd))+ labs(title ="Concordance between 2 sBBCs")
p_bar_seer

p_incidence_concordance_curie_SEER <- plot_grid(   p_bar_curie + labs(title =""),
                                                   p_bar_seer + labs(title =""),ncol = 2,labels = "AUTO",rel_widths = c(5,5))
p_incidence_concordance_curie_SEER

# Figure immune infiltration (TILs) 
## stromal TILs (Fig1C)

with_preNAC_strTILS         <- tumor_bilat_no_is %>% filter(!is.na(str_til_perc)) %>% nrow() # 277
with_preNAC_IT_TILS         <- tumor_bilat_no_is %>% filter(!is.na(it_til_perc)) %>% nrow() # 275
patients_with_str_pairs     <- mat_patient_no_is %>% filter(!is.na(str_til_perc_D),!is.na(str_til_perc_G)) 
patients_with_IT_pairs      <- mat_patient_no_is %>% filter(!is.na(it_til_perc_D),!is.na(it_til_perc_G)) 

dataf           <- tumor_bilat_no_is %>% filter(NUMDOS %in% patients_with_str_pairs$NUMDOS,!is.na(concordance_subtype),!is.na(subtype))

labeller_TILs = c(  "luminal" = "luminal  (n=218)","TNBC"  = "TNBC (n=17)","HER2+" = "HER2+ (n=15)")

p_str_TILs_subtype_concordant_discordant <- ggplot(dataf,aes(x=concordance_subtype,  str_til_perc,fill=concordance_subtype)) +
            geom_boxplot() +theme_bw()+   theme(axis.ticks = element_blank(), legend.position = "none")+   coord_cartesian(ylim= c(0,75))+
            facet_grid( ~ subtype, labeller = labeller(subtype=labeller_TILs), scales = "free", space = "free") +
            stat_compare_means(label = "p.format", label.x.npc="center", method = "wilcox.test") +
            scale_x_discrete(labels= c("Tumor in \n concordant \n pair ","Tumor in \n discordant \n pair"))+
            labs(title="Stromal TIL levels", 
             caption="Pinteraction BC subtype and concordance status on str TILS=0.57",x=" ", y="Str TIL levels (%)")

## IT TILs (Fig1D)
dataf <- tumor_bilat_no_is %>% filter(NUMDOS %in% patients_with_IT_pairs$NUMDOS, !is.na(concordance_subtype),!is.na(subtype))

labeller_TILs = c(  "luminal" = "luminal (n=218)","TNBC"  = "TNBC (n=15)","HER2+" = "HER2+ (n=15)")

p_IT_TILs_subtype_concordant_discordant <- ggplot(dataf,aes(x=concordance_subtype,it_til_perc, fill=concordance_subtype) ) +
                                geom_boxplot() + theme_bw()+ theme(axis.ticks = element_blank(), legend.position = "none")+
                                coord_cartesian(ylim= c(0,65))+
                                scale_x_discrete(labels= c("Tumor in \n concordant \n pair ","Tumor in \n discordant \n pair"))+
                                facet_grid(~subtype, labeller = labeller(subtype=labeller_TILs), scales = "free", space = "free") +
                                stat_compare_means(label = "p.format", label.x.npc="center",method = "wilcox.test") +
                                labs(title="Intratumoral TIL levels", 
                                     caption="Pinteraction BC subtype and concordance status on IT TILS=0.0006",x=" ",
                                     y="Intratumoral TIL levels (%)")

p_str_IT_TILs_subtype_concordant_discordant_no_pcr <- plot_grid( p_str_TILs_subtype_concordant_discordant,
                                                                 p_IT_TILs_subtype_concordant_discordant,
                                                                 ncol = 2,labels = c("C","D"),rel_widths = c(5,5))
# Figures pCR

## plots pCR from Curie  (FigE)

tumor_bilat_no_is_tmp        <- tumor_bilat_no_is %>% filter(pCR %in% c("pCR","No pCR"))
tumor_bilat_no_is_tmp$pCR    <- factor(tumor_bilat_no_is_tmp$pCR, levels=c("pCR","No pCR"))

pCR_same_subtype <- tumor_bilat_no_is_tmp %>%   filter(!is.na(concordance_subtype) ) %>%
                      group_by(subtype,concordance_subtype,pCR) %>%  dplyr::summarise(count=n()) %>%
                      mutate(ypos = cumsum(count) - 0.5*count) %>%
                      mutate(percent_full = count/sum(count)) %>%
                      mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>%
                      mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full) %>% 
                      mutate(n_and_perc   = paste0("n=", count,", (",percent_format,")" ) ) %>% 
                      ungroup() 
pCR_same_subtype$pCR2 <- factor(pCR_same_subtype$pCR, levels= rev(levels(as.factor(pCR_same_subtype$pCR)) ))

tumor_bilat_no_is_tmp %>% filter(subtype == "luminal") %>% select(concordance_subtype, pCR) %>% table() %>% chisq.test(correct = TRUE) # 0.08 # 
tumor_bilat_no_is_tmp %>% filter(subtype == "luminal") %>% select(concordance_subtype, pCR) %>% table() %>% chisq.test(correct = FALSE) # 0.03 # 
tumor_bilat_no_is_tmp %>% filter(subtype == "TNBC")    %>% select(concordance_subtype, pCR) %>% table() %>% chisq.test(correct = TRUE) # 0.43

tumor_bilat_no_is_tmp %>% filter(subtype == "luminal") %>% select(concordance_subtype, pCR) %>% table() %>% fisher.test # 0.05
tumor_bilat_no_is_tmp %>% filter(subtype == "TNBC") %>% select(concordance_subtype, pCR) %>% table() %>% fisher.test # 0.43

Dat.label_pcr_subtype                  <- pCR_same_subtype %>% filter(pCR == "pCR")
Dat.label_pcr_subtype$n_and_perc_total <- c("4/70 \n (6%)","4/18 \n (22%)",
                                            "8/16 \n (50%)","3/11 \n (27%)","3/9 \n (33%)")
Dat.label_pcr_subtype$pvalue           <- c("          p=0.03"," ","          p=0.43"," "," ")

p_pCR_rates_by_concordance_subtype_new2 <- ggplot(Dat.label_pcr_subtype, aes(concordance_subtype, y=percent_full*100  )   ) + 
            geom_bar(stat="identity", position="stack",  width=0.8, aes(fill = concordance_subtype)) +
            geom_text(data= Dat.label_pcr_subtype, aes( y = ypos_percent*100 ,label=n_and_perc_total), size=3.5 ) +
            geom_text(data= Dat.label_pcr_subtype, aes( y = 55 ,label=pvalue), size=3.5 ) +
            scale_x_discrete(labels= c("Tumor in \n concordant \n pair ","Tumor in \n discordant \n pair"))+
            labs(title="Response to neoadjuvant treatment", 
                 caption="Pinteraction between BC subtype and concordance status=0.03",x=" ",y="% pCR")+
            theme_bw() +theme( axis.ticks.x = element_blank(), legend.position = "none") +
            facet_grid(~subtype,labeller = labeller(subtype = c("luminal" = "luminal n=88", "TNBC"    = "TNBC n=27",
                                                       "HER2+"   = "HER2+ n=9") ) )

## plots pCR from the SEER  (FigF)
seer_bilat_neo$pCR         <- factor(seer_bilat_neo$pCR, levels=c("pCR","No pCR"))
seer_bilat_neo$subtype     <- factor(seer_bilat_neo$subtype, levels=c("Luminal","TNBC","HER2+"))
seer_bilat_neo$subtype4    <- factor(seer_bilat_neo$subtype4, levels=c("Luminal","TNBC","HER2+/HR+","HER2+/HR-"))

pCR_same_subtype <- seer_bilat_neo %>% filter(!is.na(concordance_subtype) ) %>%
                    group_by(subtype,concordance_subtype,pCR) %>%dplyr::summarise(count=n()) %>%
                    mutate(ypos = cumsum(count) - 0.5*count) %>%
                    mutate(percent_full = count/sum(count)) %>%
                    mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>%
                    mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full) %>% 
                    mutate(n_and_perc   = paste0("n=", count,", (",percent_format,")" ) ) %>% 
                    ungroup() 
pCR_same_subtype$pCR2 <- factor(pCR_same_subtype$pCR, levels= rev(levels(as.factor(pCR_same_subtype$pCR)) ))
Dat.label_pcr_subtype <- pCR_same_subtype %>% filter(pCR == "pCR")

Dat.label_pcr_subtype$n_and_perc_total <- c("213/452 \n (47%)","107/158 \n (68%)",
                                            "47/76 \n (62%)","41/67 \n (61%)",
                                            "45/68 \n (66%)","62/113 \n (55%)")
Dat.label_pcr_subtype$pvalue           <- c("          p=0.0001"," ","          p=0.99"," ","          p=0.17"," ")

seer_bilat_neo %>% filter(subtype == "Luminal") %>% select(concordance_subtype, pCR) %>% table() %>% chisq.test() # 10-5
seer_bilat_neo %>% filter(subtype == "TNBC")    %>% select(concordance_subtype, pCR) %>% table() %>% chisq.test() # 0.99
seer_bilat_neo %>% filter(subtype == "HER2+")    %>% select(concordance_subtype, pCR) %>% table() %>% chisq.test() # 0.17

p_NAC_SEER <- ggplot(Dat.label_pcr_subtype, aes(concordance_subtype, y=percent_full*100  )   ) + 
                geom_bar(stat="identity", position="stack",  width=0.8, aes(fill = concordance_subtype)) +
                geom_text(data= Dat.label_pcr_subtype, aes( y = ypos_percent*100 ,label=n_and_perc_total), size=3.5 ) +
                geom_text(data= Dat.label_pcr_subtype, aes( y = 75 ,label=pvalue), size=3.5 ) +
                labs(title="Response to NAC (Post-NAC axillar involment)", 
                     caption="Pinteraction between BC subtype and concordance status=0.001",x=" ",y="% axillar pCR")+
                theme_bw() +theme( axis.ticks.x = element_blank(), legend.position = "none") +coord_cartesian(ylim = c(0,80))+
                facet_grid(~subtype,labeller = labeller(subtype = c("Luminal" = "luminal n=610", "TNBC"    = "TNBC n=143","HER2+"   = "HER2+ n=181") ) ) +
                scale_x_discrete(labels= c("Tumor in \n concordant \n pair ","Tumor in \n discordant \n pair"))

p_NAC_curie_SEER <- plot_grid(   p_pCR_rates_by_concordance_subtype_new2, p_NAC_SEER ,
                                 ncol = 2,labels = c("E","F"),rel_widths = c(5,5))

p_concordance_tils_pCR_3_lines <- plot_grid(p_incidence_concordance_curie_SEER,
                                            p_str_IT_TILs_subtype_concordant_discordant_no_pcr,
                                            p_NAC_curie_SEER,
                                            rel_heights = c(5,5,5),
                                            nrow = 3)

save_plot("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/Fig1_p_concordance_tils_pCR_3_lines.pdf", p_concordance_tils_pCR_3_lines, base_height=15, base_width=12)	


