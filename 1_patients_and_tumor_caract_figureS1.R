source('/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/all_Rmd_cliniques/src/setup_BC_BILAT_NEO_clinique.R', local  = TRUE)
install.packages("reshape2")
install.packages("reshape2")
library(reshape)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
conflict_prefer("get_legend", "cowplot")

load(file="~/RT2Lab/BC_BILAT_NEO/clinique/data/processed/unilat_and_bilat.RData")
nrow(unilat_and_bilat) # 18190

# And only for unique patient
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/unilat_and_bilat_unique_pat.Rdata")
nrow(unilat_and_bilat_unique_pat) # 17575

##### Data on patients

# A age

compare_mean <- unilat_and_bilat_unique_pat %>%
  select(age,bilat_or_not) %>%
  filter(complete.cases(.)) %>%
  group_by(bilat_or_not) %>%
  summarise(Mean = mean(age))

p_age_violin   <- ggplot(unilat_and_bilat_unique_pat,aes(x=bilat_or_not , y=age,fill=bilat_or_not)) + 
  geom_violin() +geom_boxplot(width=0.1)+
  stat_summary(fun.y=mean, geom="point",size=1, shape=23, color="red", fill="white")+
  stat_compare_means(  method = "anova",label = "p.format",label.y = 95,label.x = 1.3)+
  scale_fill_manual(values=nice_colors_2)+
  theme_bw() +  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ylab("Age") + xlab("") + ggtitle ("Age")
p_age_violin

# B BMI

p_bmi_violin   <- ggplot(unilat_and_bilat_unique_pat,aes(x=bilat_or_not , y=bmi,fill=bilat_or_not)) + 
  geom_violin() +
  scale_fill_manual(values=nice_colors_2)+
  geom_boxplot(width=0.1)+
  stat_summary(fun.y=mean, geom="point", size=1, shape=23, color="red", fill="white")+
  stat_compare_means(  method = "anova",label = "p.format",
                       label.y = 50,label.x = 1.3)+
  theme_bw() +  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ylab("BMI") + xlab("") + ggtitle ("BMI")
p_bmi_violin

# C. BRCA

unilat_and_bilat_unique_pat$BRCA_mut <- as.character(unilat_and_bilat_unique_pat$BRCA_mut)
unilat_and_bilat_unique_pat$BRCA_mut[which(is.na(unilat_and_bilat_unique_pat$BRCA_mut))] <- "NA"

Dat.label <- unilat_and_bilat_unique_pat %>%
  select(bilat_or_not,BRCA_mut) %>% 
  group_by(bilat_or_not,BRCA_mut) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup()

Dat.label$BRCA_mut2 <- factor(Dat.label$BRCA_mut, levels = rev(levels(as.factor(Dat.label$BRCA_mut)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=BRCA_mut2)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+
  theme_ASHP2+guides(fill=guide_legend(nrow=3,byrow=TRUE,reverse=T))+
  coord_flip()+
  scale_fill_manual(name = "", values=c( "firebrick3","deepskyblue3", "#868686"))+
  xlab("")  + ylab("Percentage of cases")+
  ggtitle("BRCA Mutation")
p_percent_whole

p_BRCA <-   p_percent_whole


# D. First treatment
Dat.label <- unilat_and_bilat_unique_pat %>%
  dplyr :: group_by(bilat_or_not,first_ttt3) %>% 
  dplyr :: summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 

Dat.label$first_ttt3bis <- factor(Dat.label$first_ttt3, levels = rev(levels(as.factor(Dat.label$first_ttt3)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=first_ttt3bis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  coord_flip()+
  theme_bw()+
  theme_ASHP2+guides(fill=guide_legend(nrow=3,byrow=TRUE,reverse=T))+
  scale_fill_manual(name="", values=rev(nice_colors_3))+
  ggtitle("First treatment")
p_percent_whole

p_first_ttt3bis <-   p_percent_whole

# E Breast surgery

Dat.label <- unilat_and_bilat %>%
  select(bilat_or_not,gestchirsein) %>% 
  filter(complete.cases(.))  %>% 
  group_by(bilat_or_not,gestchirsein) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label$gestchirseinbis <- factor(Dat.label$gestchirsein, levels = rev(levels(as.factor(Dat.label$gestchirsein)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=gestchirseinbis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  coord_flip()+
  theme_bw()+theme_ASHP2+
  guides(fill=guide_legend(nrow=3,byrow=TRUE,reverse=T))+
  scale_fill_manual(name="", values=nice_colors_3)+
  ggtitle("Breast surgery")
p_percent_whole

p_gestchirsein <-   p_percent_whole

# F Chemotherapy

Dat.label <- unilat_and_bilat_unique_pat %>%
  select(bilat_or_not,ct) %>% 
  filter(!is.na(ct) )  %>% 
  group_by(bilat_or_not,ct) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label$ctbis <- factor(Dat.label$ct, levels = rev(levels(as.factor(Dat.label$ct)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=ctbis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+theme_ASHP2+guides(fill=guide_legend(nrow=2,byrow=TRUE,reverse=T))+
  coord_flip()+
  scale_fill_manual(name = "", values=c( "firebrick3","deepskyblue3"))+
  # scale_fill_manual(name = "", values=rev(nice_colors_2))+
  ggtitle("Chemotherapy")
p_percent_whole
p_ct <-   p_percent_whole

# G Endocrine therapy
Dat.label <- unilat_and_bilat_unique_pat %>%
  select(bilat_or_not,ht) %>% 
  filter(!is.na(ht) )  %>% 
  group_by(bilat_or_not,ht) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label$htbis <- factor(Dat.label$ht, levels = rev(levels(as.factor(Dat.label$ht)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=htbis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+theme_ASHP2+guides(fill=guide_legend(nrow=2,byrow=TRUE,reverse=T))+
  coord_flip()+
  scale_fill_manual(name = "", values=c( "firebrick3","deepskyblue3"))+
  # scale_fill_manual(name = "", values=rev(nice_colors_2))+
  ggtitle("Endocrine therapy")
p_percent_whole
p_ht <-   p_percent_whole

# arrange plots on patients data

p1_pts    <- plot_grid(p_age_violin,p_bmi_violin,nrow=2,  labels = "AUTO")

pleg_p_BRCA             <- get_legend(p_BRCA)
pleg_p_p_first_ttt3bis  <- get_legend(p_first_ttt3bis)
pleg_p_gestchirsein     <- get_legend(p_gestchirsein)
pleg_p_gestchirgg       <- get_legend(p_gestchirgg)
pleg_p_p_ct             <- get_legend(p_ct)
pleg_p_p_ht             <- get_legend(p_ht)

p_BRCA             <- p_BRCA          +theme(legend.position=("none"))
p_first_ttt3bis    <- p_first_ttt3bis  +theme(legend.position=("none"))
p_gestchirsein     <- p_gestchirsein  +theme(legend.position=("none"))
p_gestchirgg       <- p_gestchirgg    +theme(legend.position=("none"))
p_ht               <- p_ht            +theme(legend.position=("none"))
p_ct               <- p_ct            +theme(legend.position=("none"))

p2_pts                <- plot_grid(p_BRCA,p_first_ttt3bis, p_gestchirsein,p_ct,p_ht, nrow=5, labels = c("C","D","E","F","G")   )                                

p_legends_1_2_pts     <- plot_grid(pleg_p_BRCA,pleg_p_p_first_ttt3bis,pleg_p_gestchirsein,pleg_p_p_ct,pleg_p_p_ht,
                                   align="hv",
                                   # axis = "l",
                                   nrow = 5)
p_pts_charact_by_bilaterality  <- p1_2  <- plot_grid(p1_pts,p2_pts,p_legends_1_2_pts,nrow=1,ncol=3, rel_widths = c(1.6,2.1,1.2)) 


##### Data on tumors

# H ER staining (%)
p_er_boxplot   <- ggplot(unilat_and_bilat,aes(x=bilat_or_not , y=ROPCT,fill=bilat_or_not)) + 
  geom_boxplot() +  scale_fill_manual(values=nice_colors_2)+  theme_bw() +
  stat_compare_means(  method = "anova",label = "p.format",label.y = 115,label.x = 1.3)+
  coord_cartesian(ylim = c(0,130)) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ylab("ER staining (%)") + xlab("") + ggtitle ("ER staining (%)")
p_er_boxplot

# I PR staining (%)
p_pr_boxplot   <- ggplot(unilat_and_bilat,aes(x=bilat_or_not , y=RPPCT,fill=bilat_or_not)) + 
  geom_boxplot() +  scale_fill_manual(values=nice_colors_2)+
  theme_bw() +  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_compare_means(  method = "anova",label = "p.format", label.y = 115,label.x = 1.3)+
  coord_cartesian(ylim = c(0,130)) +
  ylab("PR staining (%)") + xlab("") + ggtitle ("PR staining (%)")
p_pr_boxplot

# J subtype

Dat.label <- unilat_and_bilat %>%
  select(bilat_or_not,subtype) %>% 
  filter(complete.cases(.)) %>%
  group_by(bilat_or_not,subtype) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label$subtype2 <- factor(Dat.label$subtype, levels = rev(levels(as.factor(Dat.label$subtype)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=subtype2)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+coord_flip()+
  theme_ASHP2+guides(fill=guide_legend(nrow=3,byrow=TRUE,reverse=T))+
  scale_fill_manual(name = "", values=c("deepskyblue3","darkorange","firebrick3"))+
  ggtitle("BC subtype")
p_percent_whole
p_subtype <-   p_percent_whole

# K Histology
Dat.label <- unilat_and_bilat %>%
  select(bilat_or_not,typhist) %>% 
  filter(complete.cases(.))  %>% 
  group_by(bilat_or_not,typhist) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label$typhistbis <- factor(Dat.label$typhist, levels = rev(levels(as.factor(Dat.label$typhist)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=typhistbis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+coord_flip()+
  theme_ASHP2+guides(fill=guide_legend(nrow=2,byrow=TRUE,reverse=T))+
  scale_fill_manual(name="", values=rev(nice_colors_4))+
  ggtitle("Histology")
p_percent_whole
p_typhisto_class <-   p_percent_whole

# L Grade
Dat.label <- unilat_and_bilat %>%
  select(bilat_or_not,gradeclasse) %>% 
  filter(complete.cases(.))  %>% 
  group_by(bilat_or_not,gradeclasse) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label$gradeclassebis <- factor(Dat.label$gradeclasse, levels = rev(levels(as.factor(Dat.label$gradeclasse)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=gradeclassebis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+coord_flip()+
  theme_ASHP2+guides(fill=guide_legend(nrow=3,byrow=TRUE,reverse=T))+
  scale_fill_manual(name = "", values=nice_colors_3)+
  ggtitle("Grade") 
p_percent_whole
p_gradeclasse <-   p_percent_whole

# M LVI
unilat_and_bilat$embols <- as.character(unilat_and_bilat$embols)

Dat.label <- unilat_and_bilat %>%
  filter(!is.na(embols))  %>%
  group_by(bilat_or_not,embols) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label

Dat.label$embolsbis <- factor(Dat.label$embols, levels = rev(levels(as.factor(Dat.label$embols)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=embolsbis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+coord_flip()+
  theme_ASHP2+guides(fill=guide_legend(nrow=2,byrow=TRUE,reverse=T))+
  scale_fill_manual(name = "", values=c( "firebrick3","deepskyblue3"))+
  ggtitle("LVI")
p_percent_whole
p_embols <-   p_percent_whole

# N multifocality
Dat.label <- unilat_and_bilat %>%
  group_by(bilat_or_not,multifocality) %>% 
  summarise(count=n()) %>%  
  mutate(ypos = cumsum(count) - 0.5*count) %>% 
  mutate(percent_full = count/sum(count)) %>% 
  mutate(percent_format = paste0(round(count/sum(count)*100), '%')) %>% 
  mutate(ypos_percent = cumsum(percent_full) - 0.5*percent_full)  %>%
  ungroup() 
Dat.label

Dat.label$multifocalitybis <- factor(Dat.label$multifocality, levels = rev(levels(as.factor(Dat.label$multifocality)))         )

p_percent_whole <- ggplot(Dat.label, aes(bilat_or_not, y=percent_full*100, fill=multifocalitybis)) + #
  geom_bar(stat="identity") +
  geom_text(data= Dat.label, aes( y=ypos_percent*100 ,label=percent_format), size=3.5 ) +
  theme_bw()+coord_flip()+
  theme_ASHP2+guides(fill=guide_legend(nrow=3,byrow=TRUE,reverse=T))+
  scale_fill_manual(name = "", values=rev(nice_colors_2))+
  scale_fill_manual(name = "", values=c( "firebrick3", "#868686"))+
  ggtitle("Multifocality")
p_percent_whole
p_multifocality   <-   p_percent_whole

# Assemble tumor factors
p1_tum    <- plot_grid(p_er_boxplot,p_pr_boxplot,nrow=2, labels = c("H","I"))

pleg_p_gradeclasse          <- get_legend(p_gradeclasse)
pleg_p_pT                   <- get_legend(p_pT)
pleg_p_p_subtype            <- get_legend(p_subtype)
pleg_p_p_typhisto_class     <- get_legend(p_typhisto_class)
pleg_p_p_embols             <- get_legend(p_embols)
pleg_p_p_multifocality      <- get_legend(p_multifocality)

p_gradeclasse       <- p_gradeclasse    +theme(legend.position=("none"))
p_subtype           <- p_subtype        +theme(legend.position=("none"))
p_typhisto_class    <- p_typhisto_class +theme(legend.position=("none"))
p_embols            <- p_embols         +theme(legend.position=("none"))
p_multifocality     <- p_multifocality  +theme(legend.position=("none"))

p2_tum                <- plot_grid(p_subtype, p_typhisto_class,p_gradeclasse,p_embols,p_multifocality, nrow=5, labels = c("J","K","L","M","N"))

p_legends_1_2_tum     <- plot_grid(pleg_p_p_subtype       ,
                                   pleg_p_p_typhisto_class,
                                   pleg_p_gradeclasse,
                                   pleg_p_p_embols        ,
                                   pleg_p_p_multifocality ,
                                   align="hv",
                                   # axis = "l",
                                   nrow = 5)                                

p_tumor_charact_by_bilaterality  <- plot_grid(p1_tum,p2_tum,p_legends_1_2_tum,nrow=1,ncol=3,
                                                  rel_widths = c(1.6,2.1,1.2)   )

p_pts_and_tumor_characteristics <- plot_grid(  p_pts_charact_by_bilaterality ,
                                               p_tumor_charact_by_bilaterality,
                                               nrow=2,rel_widths = c(1,1))

# save_plot(file="/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/FigS1_p_pts_and_tumor_characteristics.pdf", 
#           p_pts_and_tumor_characteristics, 
#           base_height=13, base_width=8.5)

save_plot(file="/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/FigS1_p_pts_and_tumor_characteristics.pdf", 
          p_pts_and_tumor_characteristics, 
          base_height=13, base_width=10)

