#Fix the repreprocessing of time and delays

#
library(survival)

# Survival outcomes

load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData") # 634
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData") # 317

source('~/RT2Lab/BC_BILAT_NEO/article_new/Nature_med/rebuttal_nature/Reviewer4/repreprocess_dates_delays_bilat.R', local = TRUE)

# Part 1 : patients description
#------------------------------

var_selected        <-  c( "age_cl_10_2","bmi_4cl","prev_pregnancy","menop","BRCA_mut")
names_var_selected  <- c("Age class","BMI class","Previous pregnancy","Menopausal status" ,"BRCA mutation")

# Now, see delays.
delrfs_bilat <- tumor_bilat_no_is %>% select (NUMDOS, rfs , 
                                              # rfs2, 
                                              delrfs ) %>% mutate (NUMDOS_rfs = paste0(NUMDOS,rfs)) %>% group_by(NUMDOS,rfs ) %>% dplyr :: summarise(maxdelrfs = max(delrfs))
# tumor_bilat_no_is %>% select (NUMDOS, rfs2 , delrfs2 ) %>% mutate (NUMDOS_rfs = paste0(NUMDOS,rfs2)) %>% group_by(NUMDOS_rfs ) %>% dplyr :: summarise(maxdelrfs = max(delrfs2))

# Add to mat_patient_no_is
mat_patient_no_is <-  mat_patient_no_is  %>% left_join(.,delrfs_bilat) 
head(mat_patient_no_is)
# Many patients with delay <0  
mat_patient_no_is_survival <- mat_patient_no_is %>% filter(maxdelrfs > 0) #%>% nrow() # 291

#### STANDARD COX 
source('~/RT2Lab/databases/core/00_common/src/R_functions_Nadir/functions_RT2_Nadir.R', local = TRUE)

cox_table_lev1 = cox_univariable_multivariable(ev                              = "rfs", 
                                               del                             = "maxdelrfs", 
                                               explanatory                     = var_selected, 
                                               nom_explanatory                 = names_var_selected, 
                                               mydataset                       = mat_patient_no_is, 
                                               variables_use_multivariable     = NA, 
                                               alpha_cut_show_multivariable    = 0.10, 
                                               alpha_cut_multivariable         = .2)

status_RFS_univariate_analysis_bilat_only_patients <- cox_table_lev1[[1]]

write.csv2(status_RFS_univariate_analysis_bilat_only_patients,file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/status_RFS_univariate_analysis_bilat_only_patients.csv")

# => No variable specific to patients is available

# Part 2 : Tumor description
#------------------------------
var_selected <- c("moddiag"            , "T"                ,   "N"                 ,  "subtype",   "concordance_subtype", "dcis"        ,       
"tumor_cellularity"  , "gradeclasse"      ,   "it_til_perc"       ,  "str_til_perc",   "histo_3cl", "multifocal_bin"   ,  
 "tyct")
names_var_selected <- c("diagnostic modality"  , "clinical T stage"      ,"clinical N stage"     , "BC subtype"       ,     "Concordance"    ,      
"DCIS component"       , "tumor cellularity (%)" ,"grade"                , "IT TILs (%)"      ,     "str TILs (%)"   ,      
 "histological type"   ,  "multifocality"        , "Chemotherapy regimen" )

cox_table_lev1 = cox_univariable_multivariable(ev                              = "rfs", 
                                               del                             = "delrfs", 
                                               explanatory                     = var_selected, 
                                               nom_explanatory                 = names_var_selected, 
                                               mydataset                       = tumor_bilat_no_is, 
                                               variables_use_multivariable     = NA, 
                                               alpha_cut_show_multivariable    = 0.10, 
                                               alpha_cut_multivariable         = .2)

status_RFS_univariate_analysis_bilat_only_tumors <- cox_table_lev1[[1]]

write.csv2(status_RFS_univariate_analysis_bilat_only_tumors,file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/status_RFS_univariate_analysis_bilat_only_tumors.csv")

# Check the p-value for interaction subtype - concordance 
c0    <- coxph(Surv(delrfs,rfs)~subtype + concordance_subtype, data=tumor_bilat_no_is_noNA, method="breslow")
c0bis <- coxph(Surv(delrfs,rfs)~concordance_subtype*subtype, data=tumor_bilat_no_is_noNA, method="breslow")
anova(c0,c0bis) # p=0.46

# Exact p-values for tests
c0    <- coxph(Surv(delrfs,rfs)~T, data=tumor_bilat_no_is, method="breslow"); summary(c0)
c0    <- coxph(Surv(delrfs,rfs)~N, data=tumor_bilat_no_is, method="breslow"); summary(c0)
c0    <- coxph(Surv(delrfs,rfs)~gradeclasse, data=tumor_bilat_no_is, method="breslow"); summary(c0)
c0    <- coxph(Surv(delrfs,rfs)~subtype, data=tumor_bilat_no_is, method="breslow"); summary(c0)


# For multivariate, only keep variables significant after univariate and with no missing value.
tumor_bilat_no_is_noNA <-  tumor_bilat_no_is %>% select(delrfs,rfs,subtype,T,N,gradeclasse)   %>% na.omit(.)

c0 <- coxph(Surv(delrfs,rfs)~NULL, data=tumor_bilat_no_is_noNA, method="breslow"); summary(c0)
add1(c0,~ 1 + subtype + N + gradeclasse +T , test="Chisq") 
c0 <- coxph(Surv(delrfs,rfs)~T, data=tumor_bilat_no_is_noNA, method="breslow") ; summary(c0)
add1(c0,~ 1 + subtype + N + gradeclasse +T , test="Chisq") 
c0 <- coxph(Surv(delrfs,rfs)~T + subtype, data=tumor_bilat_no_is_noNA, method="breslow"); summary(c0)
add1(c0,~ 1 + subtype + N + gradeclasse +T , test="Chisq") 
c0 <- coxph(Surv(delrfs,rfs)~T + subtype + gradeclasse, data=tumor_bilat_no_is_noNA, method="breslow") ; summary(c0)
add1(c0,~ 1 + subtype + N + gradeclasse +T , test="Chisq")  # N does not come in


# Write multivariate with coefficients
tmp_var		<- c("subtype","gradeclasse","T")
cox 			<- coxph(Surv(delrfs, rfs)~1+ subtype+gradeclasse+T, data=tumor_bilat_no_is_noNA, method="breslow")
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/Cox_table_loop_multivariate_corrected_ASHP_5.R', local = TRUE)
tmp_tab

# Survival analyses showed that clinical T stage, BC subtype and tumor grade were significantly associated with relapse free survival (TableS8).
 