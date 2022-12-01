# aucun package nommé ‘finalfit’ n'est trouvé 
# Remove multivariate analysis from table nadir
# Data seer
# Data GBG

# Response to neoadjuvant treatment 

# load files and setup
df_var_selected_annot <- read_excel("~/RT2Lab/BC_BILAT_NEO/clinique/results/df_var_selected_annot.xls")
source('~/RT2Lab/databases/core/00_common/src/R_functions_Nadir/functions_RT2_Nadir.R', local = TRUE)

load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData")
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData")

tumor_bilat_no_is_nac <- tumor_bilat_no_is %>% filter(first_ttt !="Surgery")
dim(tumor_bilat_no_is_nac) # 140 412
tumor_bilat_no_is_nac$pCR                     <-  as.factor(tumor_bilat_no_is_nac$pCR) 
head(tumor_bilat_no_is_nac)
tumor_bilat_no_is_nac[which(tumor_bilat_no_is_nac$tyct == "No"),"tyct"] <- NA
tumor_bilat_no_is_nac$pCR_bin <- ifelse(tumor_bilat_no_is_nac$pCR == "pCR",1,0)
dim(tumor_bilat_no_is_nac) # 140

tumor_bilat_no_is_tmp <- tumor_bilat_no_is_nac %>% filter(pCR %in% c("pCR","No pCR"))
dim(tumor_bilat_no_is_tmp) # 132

tumor_bilat_no_is_tmp %>% group_by(pCR) %>% count()
# pCR        n
# <fct>  <int>
# 1 No pCR   110
# 2 pCR       22

# Twenty-two tumors out of 140 tumors reached pathological complete response (pCR). 

var_selected <- c("age_cl_10_2"  ,"bmi_4cl"             ,"prev_pregnancy"      ,"menop"             ,  "BRCA_mut"           , "moddiag"            ,
"T"            ,"N"                   ,"subtype"             ,"concordance_subtype" ,"dcis"               , "tumor_cellularity"  ,
 "gradeclasse" ,"it_til_perc"         ,"str_til_perc"        ,"histo_3cl"         ,  "multifocal_bin"     , "tyct"               )

names_var_selected <- c("Age class"            , "BMI class"             ,"Previous pregnancy"   , "Menopausal status"   ,  "BRCA mutation"   ,     
"diagnostic modality"  , "clinical T stage"      ,"clinical N stage"     , "BC subtype"          ,  "Concordance"     ,     
 "DCIS component"      ,  "tumor cellularity (%)", "grade"               ,  "IT TILs (%)"        ,   "str TILs (%)"   ,      
 "histological type"   ,  "multifocality"        , "Chemotherapy regimen" )

table1 <- logisticRegressionTable(     data                           = tumor_bilat_no_is_nac,
                                       data_imputed                   = NA, 
                                       explanatory                    = var_selected , 
                                       nom_explanatory                = names_var_selected, 
                                       var_to_explain                 = "pCR",
                                       level_to_import                = "pCR",
                                       variables_use_multivariable    = NA, 
                                       variables_keep_multivariable   = NA, 
                                       perform_imputation             = F, 
                                       alpha_cut_multivariable        = 0.05, 
                                       alpha_show_multivariable       = 0.05,
                                       all_multivariable_values       = T 
)
table1

# table1_preformat    <- printTables (table1)
uni_multi_pcr_bilat_nac_net <- table1   %>% as.data.frame()
str(uni_multi_pcr_bilat_nac_net)

write.csv2(uni_multi_pcr_bilat_nac_net, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/uni_multi_pcr_bilat_nac_net.csv")

# Pre-NAT stromal TIL levels and breast cancer subtype were independently associated with the occurrence of a pCR (TableS7). 

# Interaction term (p value )
mod1 <-  glm(pCR_bin ~ subtype + concordance_subtype+concordance_subtype*subtype, data=tumor_bilat_no_is_tmp , family="binomial")
summary(mod1)
mod2 <-  glm(pCR_bin ~ subtype + concordance_subtype, data=tumor_bilat_no_is_tmp , family="binomial")
anova(mod1,mod2, test="Chisq") # 0.02487 # interaction  signiticant

# As was seen for TIL levels, the pCR rates showed a systemic effect when the contralateral tumor subtype was discordant. 

# Multivariate
mod1 <-  glm(pCR_bin ~ 1, data=tumor_bilat_no_is_tmp , family="binomial")
add1(mod1,~age_cl_10_2+       subtype+concordance_subtype+str_til_perc+it_til_perc+concordance_subtype*subtype+multifocal_bin, test="Chisq")  
# stromal TILs come int
mod1 <-  glm(pCR_bin ~ str_til_perc, data=tumor_bilat_no_is_tmp , family="binomial")
add1(mod1,~age_cl_10_2+       subtype+concordance_subtype+str_til_perc+it_til_perc+concordance_subtype*subtype+multifocal_bin, test="Chisq")  # Remove BRCA too many missing values
# no other variable in
# 
mod1 <-  glm(pCR_bin ~ str_til_perc + subtype + concordance_subtype+concordance_subtype*subtype, data=tumor_bilat_no_is_tmp , family="binomial")
summary(mod1)
mod2 <-  glm(pCR_bin ~ str_til_perc + subtype + concordance_subtype, data=tumor_bilat_no_is_tmp , family="binomial")
anova(mod1,mod2, test="Chisq") # 0.03245 interaction remains signiticant after multivariate analysis

mod_def <-  mod1 
summary(mod_def)								 
# source("~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/tab_summary_multivariate_regression_ASHP.r")
# 
# # write.xlsx(tab,file="results/prediction/whole population/multivariate_pCR_ASHP_WP.xlsx",row.names=TRUE)
# write.csv2(tab, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/multi_pcr_bilat_NAC_and_NET.csv")

#exact pvalue
mod1 <-  glm(pCR_bin ~ subtype, data=tumor_bilat_no_is_tmp , family="binomial")
summary(mod1)

# In luminal breast cancers, the pCR rate was significantly higher when the contralateral pair was of discordant subtype (22% versus 6%); 
# while no such pattern were found in the other subtypes (Pinteraction =0.03)(Fig1E). 



# Similar results were found in two independent validation cohorts. 
# In the SEER validation cohort, the difference in the rate of axillar pCR in tumors belonging to discordant pairs 
# versus in tumors belonging to concordant pairs was highly significant (68% versus 47%, p=0.00001, respectively) (Fig1F). 


# In the GBG cohort, the pCR rate in luminal breast cancers was significantly higher in tumors belonging to discordant pairs
# versus in tumors belonging to concordant pairs (30% versus 6%, p=0.0002 respectively) (FigS4). 
