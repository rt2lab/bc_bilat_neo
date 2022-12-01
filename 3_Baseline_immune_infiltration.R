################################################## 
#                       SETUP
################################################## 

# Packages

# Data and masterfiles
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData") # 626
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData") # 313
# load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_NAC_no_NET.RData")
load("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/data/raw/baseline_L_R.RData")    # Paired data L/R str TILs
load("/Users/ahamypet/RT2Lab/bc_bilat_neo_git/data/raw/baseline_L_R_it.RData") # Paired data L/R IT TILs


# Functions and colors
source('~/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer3/mes_fonctions2.R', local = TRUE)
colsubtype_cvd     <- c("luminal" = "#0072B2","TNBC"    = "#D55E00", "HER2+"   = "#009E73")
colsubtype_cvd_pCR <- c("luminal" = "#0072B2","TNBC"    = "#D55E00", "HER2+"   = "#009E73","pCR"     = "#6ABBEB","No pCR"  = "#044970","different pCR status"  = "#F5AE78") 

# See former script 3_immune_infiltration_clinical

# Immune infiltration levels before treatment were assessed by the presence of a mononuclear cells infiltrate 
# following the recommendations of the international TILs Working Group on hematoxylin and eosin-stained 
# sections in 149 patients (277 tumors). 


with_preNAC_strTILS         <- tumor_bilat_no_is %>% filter(!is.na(str_til_perc)) %>% nrow() # 277
with_preNAC_IT_TILS         <- tumor_bilat_no_is %>% filter(!is.na(it_til_perc)) %>% nrow() # 275

p_str_tils_before_and_concordance_no_facet_L_R_lines  <-  ggplot( data =  baseline_L_R %>% filter(!is.na(concordance_subtype)) )  +
        theme_bw()+    theme(axis.ticks.y = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 90))+
        scale_color_manual(name = " ",values = c("#0072B2", "#D55E00","#009E73") ) +
        scale_shape_manual(name = " ",values = c(17,8), labels= c("Right","Left")) +
        facet_grid( . ~ concordance_subtype, scales = "free", space = "free",#)   +
                    labeller = labeller( concordance_subtype = c("Concordant" = "Concordant pair",
                                                                 "Discordant" = "Discordant pair"))) +
        geom_point( size = 2,   aes	(x=new_NUMDOS,y=str_til_perc ,color = subtype,shape = cote)) +
        geom_line(aes	(x=new_NUMDOS,y=str_til_perc, group = new_NUMDOS),size = 0.4,
                  linetype = "dotted",position = position_dodge(width=0.1)) +#,
        ylab("Stromal TIL levels (%)") +   xlab(" ") +   ggtitle("Comparison of left (L) and right (R) stromal TIL levels ")

p_str_tils_baseline_L_R <- baseline_L_R %>%
        filter(!is.na(concordance_subtype)) %>% #nrow()
        ggplot(aes(x = concordance_subtype , y = delta_str_til_L_R_abs, fill = concordance_subtype)) +
        geom_boxplot(aes(x = concordance_subtype , y = delta_str_til_L_R_abs, fill = concordance_subtype)) +
        theme_bw() +theme(legend.position = "none", axis.ticks.x = element_blank() )+
        stat_compare_means( label = "p.format", label.x.npc="center", label.y = 40 , method = "wilcox.test") + 
        xlab("")+ ylab("Difference str TIL levels (L / R) (%)")+   ggtitle(" ")

p_it_tils_before_and_concordance_no_facet_L_R_lines  <-  ggplot( data =  baseline_L_R_it %>% filter(!is.na(concordance_subtype)) )  +
        theme_bw()+    theme(axis.ticks.y = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 90))+
        scale_color_manual(name = " ",values = c("#0072B2", "#D55E00","#009E73") ) +
        scale_shape_manual(name = " ",values = c(17,8), labels= c("Right","Left")) +
        facet_grid( . ~ concordance_subtype, scales = "free", space = "free",#)   +
                    labeller = labeller( concordance_subtype = c("Concordant" = "Concordant pair",
                                                                 "Discordant" = "Discordant pair"))) +
        geom_point( size = 2,   aes	(x=new_NUMDOS,y=it_til_perc ,color = subtype,shape = cote)) +
        geom_line(aes	(x=new_NUMDOS,y=it_til_perc, group = new_NUMDOS),size = 0.4,linetype = "dotted",position = position_dodge(width=0.1)) +
        ylab("Intratumoral TIL levels (%)") +   xlab(" ") +   ggtitle("Comparison of left (L) and right (R) intratumoral TIL levels ")

p_it_tils_baseline_L_R <- baseline_L_R_it %>% mutate(all = "all") %>% 
        filter(!is.na(concordance_subtype)) %>%
        ggplot(aes(x = concordance_subtype , y = delta_it_til_L_R_abs, fill = concordance_subtype)) +
        geom_boxplot(aes(x = concordance_subtype , y = delta_it_til_L_R_abs, fill = concordance_subtype)) +
        theme_bw() +theme(legend.position = "none", axis.ticks.x = element_blank() )+
        stat_compare_means( label = "p.format", label.x.npc="center", label.y = 40 ) + #, method = "wilcox.test") +
        xlab("")+ ylab("Difference intratumoral TIL levels (L / R) (%)")+   ggtitle(" ")

# Compil
p_indiv_str_it_l_r_compil     <-   plot_grid(p_str_tils_before_and_concordance_no_facet_L_R_lines,
                                             p_str_tils_baseline_L_R+theme(plot.margin = unit(c(0.6,0,2.5,0), "lines")),
                                             p_it_tils_before_and_concordance_no_facet_L_R_lines, 
                                             p_it_tils_baseline_L_R+theme(plot.margin = unit(c(0.6,0,2.5,0), "lines")),
                                             labels = c("AUTO"),  ncol = 2, nrow = 2, rel_widths=c(10,2,10,2), align="hv")

save_plot(p_indiv_str_it_l_r_compil, 
          file = "/Users/ahamypet/RT2Lab/bc_bilat_neo_git/figures/ExtendedFig2_p_indiv_str_it_l_r_compil.pdf",base_width = 17, base_height = 9)

# The difference between the TILs levels from the left and the right tumor was higher in pairs of discordant BC subtypes 
# than in pairs of concordant BC subtypes (ExtendedFig2)  

# At the tumor level, TIL levels were independently associated higher grade tumors (TableS5 and S6), 

# BUILD TABLE S5 and S6 : 
# Lefts columns of the table : counts, median and mean

# Part 1 : patients 
#------------------------------
var_selected            <- c("age_cl_10_2","bmi_4cl" ,"menop","BRCA_mut")
names_var_selected      <- c( "Age class","BMI class","Menopausal status", "BRCA mutation")    
matching_name_file      <- data.frame(variable = var_selected, name_variable = names_var_selected)
variables_to_test       <- var_selected               
names_variables_to_test <- names_var_selected  
variable_to_compare			  <- "str_til_perc"
variable_to_compare_name	<- "Pre-NAC str TILs"
dataf  			              <- tumor_bilat_no_is ; print(dim(dataf))
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/mean_comparison_anova_and_plot_ASHP_v20.R', local = TRUE)
# source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/mean_comparison_anova_and_plot_ASHP_v20_no_round.R', local = TRUE)

mean_pre_NAC_str_TILs_value_patient <- tab_reg_linear %>% as.data.frame()

variable_to_compare			  <- "it_til_perc"
variable_to_compare_name	<- "Pre-NAC IT TILs"
dataf  			              <- tumor_bilat_no_is ; print(dim(dataf))
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/mean_comparison_anova_and_plot_ASHP_v20.R', local = TRUE)
mean_pre_NAC_it_TILs_value_patient <- tab_reg_linear %>% as.data.frame()

# Part 2 : tumors 
#------------------------------
var_selected        <- c("moddiag"          , "T"              ,   "N"               ,  "pT"                ,"subtype"         ,  "er"               ,
"pr"               , "HR"             ,   "HER2"            ,  "ROINT"             ,"RPINT"           ,  "infilt"           ,
 "dcis"            ,  "pnuicc_3cl"    ,    "gradeclasse"    ,   "embols"           , "histo_3cl"      ,   "multifocal_bin"  , "concordance_subtype" )

names_var_selected      <- c("diagnostic modality"           ,     "clinical T stage"                ,   "clinical N stage",
"pathological T stage"          ,     "BC subtype"                      ,   "ER status"                         ,
"PR status"                     ,     "HR status"                       ,   "HER2 status"                       ,
 "intensity of ER positivity"   ,      "intensity of PR positivity"     ,    "Invasive or DCIS"                  ,
 "DCIS component"               ,      "pN status"                      ,    "grade"                             ,
 "lymphovascular invasion"      ,      "histological type"              ,    "multifocality"         ,"Concordance"            )

matching_name_file      <- data.frame(variable = var_selected, name_variable = names_var_selected)
variables_to_test       <- var_selected               
names_variables_to_test <- names_var_selected  

variable_to_compare			  <- "str_til_perc"
variable_to_compare_name	<- "Pre-NAC str TILs"
dataf  			              <- tumor_bilat_no_is ; print(dim(dataf))
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/mean_comparison_anova_and_plot_ASHP_v20_no_round.R', local = TRUE)
mean_pre_NAC_str_TILs_value_tumor <- tab_reg_linear %>% as.data.frame()
write.csv2(mean_pre_NAC_str_TILs_value_tumor, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/mean_pre_NAC_str_TILs_value_tumor.csv")


variable_to_compare			  <- "it_til_perc"
variable_to_compare_name	<- "Pre-NAC IT TILs"
dataf  			              <- tumor_bilat_no_is ; print(dim(dataf))
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/mean_comparison_anova_and_plot_ASHP_v20_no_round.R', local = TRUE)
# source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/mean_comparison_anova_and_plot_ASHP_v20.R', local = TRUE)
mean_pre_NAC_it_TILs_value_tumor <- tab_reg_linear %>% as.data.frame()
write.csv2(mean_pre_NAC_it_TILs_value_tumor, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/mean_pre_NAC_it_TILs_value_tumor.csv")


# Agregate tables
median_pre_NAC_str_TILs_value <- bind_rows(mean_pre_NAC_str_TILs_value_patient,mean_pre_NAC_str_TILs_value_tumor) 
write.csv2(median_pre_NAC_str_TILs_value, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/median_pre_NAC_str_TILs_value.csv")

median_pre_NAC_it_TILs_value <- bind_rows(mean_pre_NAC_it_TILs_value_patient,mean_pre_NAC_it_TILs_value_tumor) 
write.csv2(median_pre_NAC_it_TILs_value, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/median_pre_NAC_it_TILs_value.csv")

# Middle and right panels of the table S5 and S6 :  univariate and multivariate analysis
# from code : Q8_R3_tumor_cellularity_and_TILs

############################################################################################################################
# Uni and multi str TILS
############################################################################################################################

##  Univariate
var_selected      <- tumor_bilat_no_is %>% select(age_cl_10_2, bmi_4cl , 
                                                  menop, BRCA_mut, moddiag,T,N,pT,pnuicc_3cl,dcis,
                                                  tumor_cellularity,
                                                  gradeclasse,embols,histo_3cl,multifocal_bin,subtype,concordance_subtype) %>%
                                          colnames() %>% as.matrix() %>% as.character()
uni               <- univariee(var_selected, FUN = lineaire, formule = str_til_perc ~ 1, data = tumor_bilat_no_is)
uni_lm_str_tils   <- uni
write.csv2(uni_lm_str_tils, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/uni_lm_str_tils.csv")

# Exact p-values
mod0 <-  lm(str_til_perc ~ pnuicc_3cl  , data=tumor_bilat_no_is) ; summary(mod0)
mod0 <-  lm(str_til_perc ~ gradeclasse , data=tumor_bilat_no_is) ; summary(mod0)
mod0 <-  lm(str_til_perc ~ subtype , data=tumor_bilat_no_is) ; summary(mod0)
mod0 <-  lm(str_til_perc ~ concordance_subtype  , data=tumor_bilat_no_is) ; summary(mod0)


# interaction  between stromal TILs / concordance / subtype  NS
mod0 <-  lm(str_til_perc ~ concordance_subtype+subtype  , data=tumor_bilat_no_is) ; summary(mod0)
mod1 <-  lm(str_til_perc ~ concordance_subtype+subtype  + concordance_subtype*subtype, data=tumor_bilat_no_is) ; summary(mod1)
anova(mod0,mod1, test="Chisq")  # P=0.57

## Multivariate

mod0 <-  lm(str_til_perc ~ NULL  , data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype+ tumor_cellularity, test="Chisq")
mod0 <-  lm(str_til_perc ~ gradeclasse  , data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype+ tumor_cellularity, test="Chisq")
mod0 <-  lm(str_til_perc ~ gradeclasse + concordance_subtype , data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype+ tumor_cellularity, test="Chisq")
mod0 <-  lm(str_til_perc ~ gradeclasse + concordance_subtype +moddiag , data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype+ tumor_cellularity, test="Chisq") # Final model

summary(mod0)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      7.734      1.359   5.692 3.34e-08 ***
#   gradeclasse2                     2.079      1.650   1.260  0.20873    
# gradeclasse3                    10.561      2.131   4.955 1.29e-06 ***
#   concordance_subtypeDiscordant    6.065      1.949   3.111  0.00207 ** 
#   moddiagPalpable                  3.042      1.475   2.062  0.04017 *  

# and interestingly, the relationship between TIL levels and breast cancer subtype showed a systemic effect, i.e.,
# it was affected by the subtype of the contralateral tumor. 
# In luminal breast cancers, stromal (Str) TIL levels were lower when the subtype of the contralateral tumor was concordant 
# than when it was discordant, and the same trend was observed for intratumoral TILs (Fig1D-E). 

############################################################################################################################
# Univariate and multivariate IT TILS
############################################################################################################################

## Univariate
uni               <- univariee(var_selected, FUN = lineaire, formule = it_til_perc ~ 1, data = tumor_bilat_no_is)
uni_lm_it_tils    <- uni
write.csv2(uni_lm_it_tils, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/uni_lm_it_tils.csv")

# Interaction very significative
mod0 <-  lm(it_til_perc ~ concordance_subtype+subtype  , data=tumor_bilat_no_is) ; summary(mod0)
mod1 <-  lm(it_til_perc ~ concordance_subtype+subtype  + concordance_subtype*subtype, data=tumor_bilat_no_is) ; summary(mod1)
anova(mod0,mod1, test="Chisq")  # P=0.0006

# Exact p-values
mod0 <-  lm(it_til_perc ~ gradeclasse , data=tumor_bilat_no_is) ; summary(mod0)
mod0 <-  lm(it_til_perc ~ subtype , data=tumor_bilat_no_is) ; summary(mod0)
mod0 <-  lm(it_til_perc ~ concordance_subtype  , data=tumor_bilat_no_is) ; summary(mod0)

## Multivariate with tumor cellularity
mod0 <-  lm(it_til_perc ~ NULL  , data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype + concordance_subtype*subtype + tumor_cellularity, test="Chisq")#+ pnuicc_3cl)#) + embols_post
mod0 <-  lm(it_til_perc ~ gradeclasse  , data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype + concordance_subtype*subtype + tumor_cellularity, test="Chisq")#+ pnuicc_3cl)#) + embols_post
# Interaction comes in 
mod0 <-  lm(it_til_perc ~ gradeclasse  + subtype+ concordance_subtype + concordance_subtype*subtype, data=tumor_bilat_no_is) ; summary(mod0)
add1(mod0,~ age_cl_10_2+ moddiag+ gradeclasse+ histo_3cl+ subtype+ concordance_subtype + concordance_subtype*subtype + tumor_cellularity, test="Chisq")#+ pnuicc_3cl)#) + embols_post

# Final model
mod0 <-  lm(it_til_perc ~ gradeclasse  + subtype+ concordance_subtype + concordance_subtype*subtype, data=tumor_bilat_no_is) ; summary(mod0)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                  3.7953     0.7881   4.816 2.50e-06 ***
#   gradeclasse2                                 1.0241     0.9742   1.051 0.294127    
# gradeclasse3                                 6.8778     1.3122   5.241 3.31e-07 ***
#   subtypeTNBC                                 13.2903     3.4779   3.821 0.000166 ***
#   subtypeHER2+                                 5.4672     3.4194   1.599 0.111072    
# concordance_subtypeDiscordant                3.4139     1.5478   2.206 0.028282 *  
#   subtypeTNBC:concordance_subtypeDiscordant  -18.1955     4.2366  -4.295 2.47e-05 ***
#   subtypeHER2+:concordance_subtypeDiscordant  -4.2209     4.1934  -1.007 0.315096    


# Conversely, in TNBCs, the intratumoral TIL levels were lower when the subtype of the contralateral tumor 
# was concordant than when it was discordant. The interaction test was highly significant (Pinteraction=0.0006), 
# indicating that the impact of tumor subtype on intratumoral immune infiltration was significantly modified by the 
# concordance of the breast cancer subtype of the tumor pair it belonged to. 
# This result was also validated in a third independent cohort from the German breast group (GBG), 
# where the interactions tests were highly significant both for stromal and intratumoral TILs 
# (Pinteraction=0.007 and Pinteraction=0.006 respectively) (FigS2). 
# This suggests that TIL levels are not purely determined by local tumor microenvironment properties. 
