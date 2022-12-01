################################################## 
#                       SETUP
################################################## 

rm(list = ls())

# Packages
library(corrplot)    
library(matrixcalc)
library(vegan)
library(fmsb)

# Data and masterfiles
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData") # 626
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData") # 313
df_var_selected_annot <- read_excel("~/RT2Lab/BC_BILAT_NEO/clinique/results/df_var_selected_annot.xls")

# Functions and colors
source('~/RT2Lab/BC_BILAT_NEO/clinique/all_Rmd_cliniques/src/setup_BC_BILAT_NEO_clinique.R', local  = TRUE)

# Overall, the 313 paired sBBC tumors shared more common characteristics than expected by chance (TableS4)

var_selected    <-  c("moddiag"                   ,"tclin"                   ,  "T"       ,   "N"   ,
"pT"                        ,"subtype"                 ,  "er"                     ,   "pr"                       ,
"HR"                        ,"HER2"                    ,  "ROINT"                  ,   "ROPCT"                    ,
 "RPINT"                    , "RPPCT"                  ,  
"it_til_perc"             , 
 "str_til_perc"             , "tumor_cellularity"      ,   "mitotic_index"         ,    "perc_stroma"             , 
 "dcis"                     , "tailhist"               ,   "nbggpos"               ,    "pnuicc_3cl"              , 
 "gradeclasse"              , "NBMIT"                  ,   "embols"                ,    "histo_3cl"               , 
######
"gestchirsein"          ,    
#####
 "multifocal_bin"           , "multifocality"          ,   "it_til_perc_postneo"   ,    "str_til_perc_postneo"    , 
 "tumor_cellularity_postneo", "mitotic_index_postneo"  , 
"pCR"                     , 
 "rcb"                      , "rcb_class"              ,   "rcb_class_integer"     ,    "nbggpos_postneo"         , 
 "ypnuicc_3cl"              , "embols_post"            ,   "gestchirsein"          ,    "gestgg" )

names_var_selected  <- c("diagnostic modality"                , "clinical size"                       ,"clinical T stage"                   ,
"clinical N stage"                   , "pathological T stage"                ,"BC subtype"                         ,
"ER status"                          , "PR status"                           ,"HR status"                          ,
 "HER2 status"                       ,  "intensity of ER positivity"         , "percentage of ER positivity"       , 
 "intensity of PR positivity"        ,  "percentage of PR positivity"        , "IT TILs (%)"                       , 
 "str TILs (%)"                      ,  "tumor cellularity (%)"              , "mitotic index"                     , 
 "stroma cellularity (%)"            ,  "DCIS component"                     , "histological size"                 , 
######
"breast surgery"           ,
#####
 "number of positive nodes"          ,  "pN status"                          , "grade"                             , 
 "number mitoses"                    ,  "lymphovascular invasion"            , "histological type"                 , 
 "multifocality"                     ,  "multifocality"                      , "post-NAC IT TILs (%)"              , 
 "post-NAC str TILs (%)"             ,  "post-NAC tumor cellularity (%)"     , "post-NAC mitotic index"            , 
 "pCR status"                        ,  "RCB"                                , "RCB (class)"                       , 
 "RCB (class)"                       ,  "number of positive nodes (post NAC)", "post-NAC node (class)"             , 
 "lymphovascular invasion (post-NAC)",  "breast surgery"                     , "axillar surgery" )

kendall_L_R <- data.frame(  variable     = var_selected,
                            names        = names_var_selected,
                            variable_type = NA,
                            concordance_perc = "",
                            kappa             = "",
                            pval_kappa        = NA,
                            Kendall           = "",
                            pval_Kendall      = NA,
                            pearson           = "",
                            pval_person       = NA,
                            spearman          = "",
                            pval_spearman     = NA)


# pdf("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/all_concordances.pdf", height=10, margin(t=2))
for (i in 1: length(var_selected) ){
  # i=41 
  tmp_var <- var_selected[i]
  print(tmp_var)
  tmp_name_var <- df_var_selected_annot[match(tmp_var,df_var_selected_annot$variable),"name_variable"]
  print(i)
  tmp_var_G <- paste0(tmp_var,"_G")
  tmp_var_D <- paste0(tmp_var,"_D")
  
  tmp_mat <- mat_patient_no_is[,c(tmp_var_G,tmp_var_D)]
  # tmp_mat <- mat_patient[,c(tmp_var_G,tmp_var_D)]
  colnames(tmp_mat)[1:2] <- c("tmp_var_G","tmp_var_D") 
  tmp_mat <- tmp_mat %>% filter(!is.na(tmp_var_G), !is.na(tmp_var_D)) 
  
  # Calculate Kendall                
  concordance_tmp_var      <- kendall.global(tmp_mat, nperm = 999, mult = "holm")
  kendall_L_R[i,"Kendall"] <- round(concordance_tmp_var$Concordance_analysis["W","Group.1"],2) # 0.62
  kendall_L_R[i,"pval_Kendall"]    <- concordance_tmp_var$Concordance_analysis["Prob.F","Group.1"]
  
  # Calculate kappa
  if(class(tmp_mat[,"tmp_var_G"])=="character" ) {
    kendall_L_R[i,"variable_type"] <- "character"
    tmp_levels                  <- c(levels(as.factor(tmp_mat[,1])), levels(as.factor(tmp_mat[,2]))) %>% unique()
    tmp_mat[,1]                 <- factor(tmp_mat[,1], levels = tmp_levels)
    tmp_mat[,2]                 <- factor(tmp_mat[,2], levels = tmp_levels)
    tmp_contingence_table       <- table(tmp_mat[,1],tmp_mat[,2])
    tmp_sum_diagonale           <- sum(diag(tmp_contingence_table))
    tmp_ref_class                <- names(which.max(diag(tmp_contingence_table))) 
    other_class                  <- setdiff(colnames(tmp_contingence_table), tmp_ref_class)
    reordered_vector             <- c(tmp_ref_class,other_class)
    tmp_contingence_table        <- tmp_contingence_table[reordered_vector,reordered_vector]
    tmp_sum_matrix               <- sum(tmp_contingence_table)
    kendall_L_R[i,"concordance_perc"]   <- round((tmp_sum_diagonale*100/tmp_sum_matrix),1) 
    
    # Calculate kappa        
    kappa_tmp_var                 <- Kappa.test(tmp_contingence_table)
    kendall_L_R[i,"kappa"]        <- round(kappa_tmp_var$Result[["estimate"]],2) # 0.62
    # kendall_L_R[i,"pval_kappa"]   <- round(kappa_tmp_var$Result[["p.value"]],2) 
################################################################################################    
    kendall_L_R[i,"pval_kappa"]   <- kappa_tmp_var$Result[["p.value"]] 
################################################################################################    
    
    # Create pseudo correlation matrix 
    tmp_false_diag         <- upper.triangle(tmp_contingence_table) + upper.triangle(t(tmp_contingence_table))  
    diag(tmp_false_diag)   <- diag(tmp_contingence_table)
    tmp_true_diag          <- tmp_false_diag 
    
    # install.packages("ggcorrplot")                    
    #    p.mat <- cor_pmat(tmp_true_diag)
    #    library("ggcorrplot")                    
    #    ggcorrplot(p.mat, hc.order = TRUE, type = "lower",
    # outline.col = "white")
    
    tmp_title <- unname(tmp_name_var)
    col3      <- colorRampPalette(c("white", "darkorange"))
    corrplot(tmp_true_diag, type="upper",   
             col=col3(50),is.corr = FALSE,addCoef.col = TRUE,
             title = tmp_title,
             # tl.srt = 30,
             cl.pos = "n", number.cex=0.75,tl.col="black") # ,) #
  }
  
  if(class(tmp_mat[,"tmp_var_G"]) =="numeric" ) {
    kendall_L_R[i,"variable_type"] <- "numeric"
    corr_test_tmp_var <- cor.test(tmp_mat[,"tmp_var_G"],tmp_mat[,"tmp_var_D"])
    kendall_L_R[i,"pearson"]     <- round(corr_test_tmp_var$estimate,2) # 
    kendall_L_R[i,"pval_person"] <- corr_test_tmp_var$p.value# 
    
    corr_test_spearman_tmp_var       <- cor.test(tmp_mat[,"tmp_var_G"],tmp_mat[,"tmp_var_D"],method="spearman")
    kendall_L_R[i,"spearman"]     <- round(corr_test_spearman_tmp_var$estimate,2) # 
    kendall_L_R[i,"pval_spearman"] <- corr_test_spearman_tmp_var$p.value# 
  }
  
  if(class(tmp_mat[,"tmp_var_G"]) =="integer" ) {
    kendall_L_R[i,"variable_type"] <- "integer"
  }
}
# dev.off()
write.csv2(kendall_L_R, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/kendall_L_R.csv")


kendall_L_R$pval_kappa    <- ifelse(as.numeric(kendall_L_R$pval_kappa   ) <=0.001, "<0.001", round(kendall_L_R$pval_kappa,3)  )
kendall_L_R$pval_Kendall  <- ifelse(as.numeric(kendall_L_R$pval_Kendall ) <=0.001, "<0.001", round(kendall_L_R$pval_Kendall,3)  )
kendall_L_R$pval_person   <- ifelse(as.numeric(kendall_L_R$pval_person  ) <=0.001, "<0.001", round(kendall_L_R$pval_person,3)  )
kendall_L_R$pval_spearman <- ifelse(as.numeric(kendall_L_R$pval_spearman) <=0.001, "<0.001", round(kendall_L_R$pval_spearman,3)  )

kendall_L_R[which(is.na(kendall_L_R$pval_kappa)),"pval_kappa"]          <- ""
kendall_L_R[which(is.na(kendall_L_R$pval_person)),"pval_person"]        <- ""
kendall_L_R[which(is.na(kendall_L_R$pval_spearman)),"pval_spearman"]    <- ""

head(kendall_L_R)
# kendall_L_R$concordant_or_discordant <- NA
kendall_L_R[which(kendall_L_R$pval_kappa),"concordant_or_discordant"]
kendall_L_R <- kendall_L_R %>% 
  mutate(concordant_or_discordant =  case_when(as.numeric(pval_kappa) < 0.051    | pval_kappa == "<0.001" ~ "concordant", 
                                               as.numeric(pval_Kendall) < 0.051  | pval_Kendall == "<0.001" ~ "concordant",
                                               as.numeric(pval_person) < 0.051   | pval_person == "<0.001" ~ "concordant",
                                               as.numeric(pval_spearman) < 0.051 | pval_spearman == "<0.001" ~ "concordant",
                                               TRUE ~ "discordant")  )

# kendall_L_R <-   kendall_L_R %>% select(variable,names,Kendall,pval2,pearson,pval3,spearman,pval4) 
write.csv2(kendall_L_R, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/kendall_L_R.csv")

#the majority (84.7%) of the tumor pairs were concordant regarding clinical and pathological patterns,
#notably regarding breast cancer subtype (Fig1A). #A minority of pairs of tumors belonged to different BC subtypes (discordant pairs: 15.3%) 

#and both the proportion of pairs (18%) and their relative repartition were similar in the validation cohort 
#of 8367 patients with sBBCs from the SEER database (Fig1B). 

# Cf codes figure 1