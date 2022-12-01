source('~/RT2Lab/BC_BILAT_NEO/clinique/all_Rmd_cliniques/src/setup_BC_BILAT_NEO_clinique.R', local  = TRUE)

# preprocessed in preprocess_data_base_curie_bilat_non_bilat.R 
load(file="~/RT2Lab/BC_BILAT_NEO/clinique/data/processed/unilat_and_bilat.RData")
load("/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/unilat_and_bilat_unique_pat.Rdata")
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/tumor_bilat_no_is.RData") # 626
load(file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/clinique/data/processed/mat_patient_no_is.RData") # 313


#  1.1 Unilat	
# ...................
nip_unique <- unilat_and_bilat %>% select(NUMDOS) %>% as.matrix() %>% as.character() %>% unique(); length(nip_unique)  # 17575
nip_unilat <- unilat_and_bilat %>% filter( bilat_or_not=="unilateral") %>% select(NUMDOS) %>% as.matrix() %>% as.character() %>% unique(); length(nip_unilat) 
# 17171

#  1.2 Bilat	
# ...................
nip_bilat <- unilat_and_bilat %>% filter( bilat_or_not=="bilateral") %>% select(NUMDOS) %>% as.matrix() %>% as.character() %>% unique(); length(nip_bilat) # 404
nrow(unilat_and_bilat %>% filter( bilat_or_not=="bilateral")) # 808

unilat_and_bilat %>% summarise(age = median(age, na.rm= TRUE))

ntot <- length(unique(unilat_and_bilat$NUMDOS)) # 17575
perc <- round(length(nip_bilat)/ntot*100,1)

# Figure flow chart

#Slight differences existed in patients and tumor characteristics between patients with unilateral breast cancers and patients with sBBCs 
# (TableS1 and FigS1) 

# Part 1 : Are patients with bilat different from patients with unilat ?
#----------------------------------------------------------------

var_selected <-  c("age"           , "age_cl_10_2"    ,"bmi"            ,"BMI_cl"         ,"bmi_4cl"        ,"menarche"      ,
"prev_pregnancy", "menop"          ,"age_menop"      ,"hrt"            ,"fam_history"    ,"BRCA_screen"   ,
 "BRCA_mut"     ,  "first_ttt3"     ,"ct"             ,"tyct"           ,"ht"             ,"typht"         ,
 "tc"           ,  "typtc"          ,"rt"           )

names_var_selected <- c("Age","Age class","BMI","BMI class","BMI class","Age at menarche",
"Previous pregnancy","Menopausal status","Age at menopause","Hormone replacement therapy",
"Familial history breast / ov. ","Research hereditary predisposition","BRCA mutation","First treatment",
"Chemotherapy","Chemotherapy regimen","Endocrine therapy","Endocrine therapy type","Targeted therapy",
"Targeted therapy type","Radiotherapy")

# Table S1 (patients = part 1)
matching_name_file    <- data.frame(variable = var_selected, name_variable = names_var_selected)
mydataset             <- unilat_and_bilat_unique_pat
Table1                <- CreateTableOne(var_selected ,"bilat_or_not",mydataset) 
table1_preformat      <- print(Table1, quote=TRUE, noSpaces=TRUE,showAllLevels = TRUE, pDigits=3, contDigits=1)
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/format_tableone_v2.R', local  = TRUE)
Table_patient_bilat_no_bilat <- Table1_format
write.csv2(Table1_format, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/TableS1_patient_bilat_no_bilat.csv")

# Legends for exact p-values 
unilat_and_bilat_unique_pat %>% select(bilat_or_not,BRCA_screen) %>% table() %>% chisq.test()
unilat_and_bilat_unique_pat %>% select(bilat_or_not,ht) %>% table() %>% chisq.test()
unilat_and_bilat_unique_pat %>% select(bilat_or_not,typht) %>% table() %>% chisq.test()

#  Part 2. Are the bilateral TUMORS different from the unilat tumors?
#----------------------------------------------------------------
var_selected <- c("moddiag"       ,"tclin"       ,  "T"             ,"N"          ,   "pT"          ,  "subtype"      , "er"            ,"pr"   ,        
"HR"            ,"HER2"        ,  "ROINT"         ,"ROPCT"      ,   "RPINT"       ,  "RPPCT"        , "infilt"        ,"dcis" ,        
 "tailhist"     , "nbggpos"    ,   "gradeclasse"  , "NBMIT"     ,    "embols"     ,   "histo_3cl"   ,  "multifocality", "pCR" ,         
 "gestchirsein" , "gestgg" )
names_var_selected <- c(
  "diagnostic modality"        , "clinical size"               ,"clinical T stage"           , "clinical N stage"           ,
  "pathological T stage"       , "BC subtype"                  ,"ER status"                  , "PR status"                  ,
  "HR status"                  , "HER2 status"                 ,"intensity of ER positivity" , "percentage of ER positivity",
   "intensity of PR positivity",  "percentage of PR positivity", "Invasive or DCIS"          ,  "DCIS component"             ,
   "histological size"         ,  "number of positive nodes"   , "grade"                     ,  "number mitoses"             ,
   "lymphovascular invasion"   ,  "histological type"          , "multifocality"             ,  "pCR status"                 ,
   "breast surgery"            ,  "axillar surgery")            
  
names_var_selected    <- var_tumor_uni_bilat$name_variable
matching_name_file    <- data.frame(variable = var_selected, name_variable = names_var_selected)
mydataset             <- unilat_and_bilat
Table1                <- CreateTableOne(var_selected ,"bilat_or_not",mydataset ) 
table1_preformat      <- print(Table1, quote=TRUE, noSpaces=TRUE,showAllLevels = TRUE, pDigits=3, contDigits=1,nonnormal = c("ROPCT","RPPCT"))
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/format_tableone_v2.R', local  = TRUE)
write.csv2(Table1_format, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/TableS1bis_patient_bilat_no_bilat.csv")

# Legends for exact p-values 
unilat_and_bilat %>% select(bilat_or_not,T) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,subtype) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,er) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,pr) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,HER2) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,typht) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,multifocality) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,gradeclasse) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,embols) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,typhist) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,histo_3cl) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,gestchirsein) %>% table() %>% chisq.test()
unilat_and_bilat %>% select(bilat_or_not,gestgg) %>% table() %>% chisq.test()



# FigS1 => separate code

# Exclude all cases with at least one in situ

nrow(tumor_bilat_no_is)
nrow(mat_patient_no_is)

# Out of 313 patients with invasive sBBCs, most of the tumors were luminal (n=538, 87.6%),
#whereas TNBC (n=44, 7.2%) and HER2-positive breast cancers (n=32, 5.2%) were rare (TableS2). 

# Part 1 : patients description
#------------------------------
var_selected <-  c("age"           , "age_cl_10_2"    ,"bmi",      "bmi_4cl"        ,
                   "prev_pregnancy", "menop"          ,"age_menop"      ,"hrt"            ,"fam_history"    ,"BRCA_screen"   ,
                   "BRCA_mut"     ,  "first_ttt3"     ,"ct"             ,"tyct"           ,"ht"             ,"typht"         ,
                   "tc"           ,"rt"           )

names_var_selected <- c("Age","Age class","BMI","BMI class",
                        "Previous pregnancy","Menopausal status","Age at menopause","Hormone replacement therapy",
                        "Familial history breast / ov. ","Research hereditary predisposition","BRCA mutation","First treatment",
                        "Chemotherapy","Chemotherapy regimen","Endocrine therapy","Endocrine therapy type","Targeted therapy",
                        "Radiotherapy")

names_var_selected <- df_var_selected_annot[match(var_selected,df_var_selected_annot$variable ),"name_variable"]
mydataset          <- mat_patient_no_is
matching_name_file <- data.frame(variable = var_selected, 
                                 name_variable = names_var_selected)
Table1                <- CreateTableOne(var_selected , ,mat_patient_no_is ) 
table1_preformat      <- print(Table1, quote=TRUE, noSpaces=TRUE,showAllLevels = TRUE, pDigits=3, contDigits=1)
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/format_tableone_v2.R', local  = TRUE)
write.csv2(Table1_format, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/TableS2table_pts_caract_no_is.csv")

# Part 2 : Tumor description
#------------------------------
var_selected <- c("moddiag"       ,"tclin"       ,  "T"             ,"N"          ,   "pT"          ,  "subtype"      , "er"            ,"pr"   ,        
                  "HER2"        ,  "ROINT"         ,"ROPCT"      ,   "RPINT"       ,  "RPPCT"        , "infilt"        ,"dcis" ,        
                  "tailhist"     , "nbggpos"    ,   "gradeclasse"  , 
                  "embols"     ,   "histo_3cl"   ,  "multifocality", 
                  "gestchirsein" , "gestgg" )
names_var_selected <- c(
  "diagnostic modality"        , "clinical size"               ,"clinical T stage"           , "clinical N stage"           ,
  "pathological T stage"       , "BC subtype"                  ,"ER status"                  , "PR status"                  ,
  "HER2 status"                 ,"intensity of ER positivity" , "percentage of ER positivity",
  "intensity of PR positivity",  "percentage of PR positivity", "Invasive or DCIS"          ,  "DCIS component"             ,
  "histological size"         ,  "number of positive nodes"   , "grade"                     ,  
  "lymphovascular invasion"   ,  "histological type"          , "multifocality"             ,  
  "breast surgery"            ,  "axillar surgery")            
matching_name_file  <- data.frame(variable = var_selected, 
                                  name_variable = names_var_selected)
mydataset            <- tumor_bilat_no_is;
Table1                <- CreateTableOne(var_selected ,, mydataset) 
table1_preformat      <- print(Table1, quote=TRUE, noSpaces=TRUE,showAllLevels = TRUE, pDigits=3, contDigits=1)
source('~/RT2Lab/databases/core/00_common/src/R_functions_ASHP/format_tableone_v2.R', local  = TRUE)
write.csv2(Table1_format, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/TableS2bis_table_pts_caract_no_is.csv")

#Only 13 patients were carriers of a genetic germline BRCA1 or BRCA2 predisposition. 
# They were significantly younger, and more likely to be diagnosed with large, palpable, high-grade tumors, more frequently of TNBC subtype (Table S3). 

# Part 1 : patients description
#------------------------------
var_selected       <-  c("age" , "age_cl_10_2"    ,"bmi_4cl","prev_pregnancy", "menop")
names_var_selected <- c("Age","Age class","BMI class","Previous pregnancy","Menopausal status")
names_var_selected <- df_var_selected_annot[match(var_selected,df_var_selected_annot$variable ),"name_variable"]
mydataset          <- mat_patient_no_is
mat_patient_no_is$BRCA_mut
matching_name_file <- data.frame(variable = var_selected, 
                                 name_variable = names_var_selected)
data_dictionary <- data.frame(var       = var_selected ,
                              names_var = names_var_selected)
mat_patient_no_is[,var_selected]

source('~/RT2Lab/databases/core/00_common/src/R_functions_Nadir/functions_RT2_Nadir.R', local = TRUE)
table1 <- table1_rt2(mydataset          = mat_patient_no_is,
                     raw_variables      = var_selected,
                     data_dictionary    = data_dictionary,
                     stratif            = "BRCA_mut",
                     missing            = FALSE,
                     perc_by_column     = FALSE,
                     apply_default_test = FALSE)
tableone_brca_patients <- table1[[1]]
write.csv2(tableone_brca_patients, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/TableS3tableone_brca_patients.csv")

# Part 2 : Tumor description
#------------------------------
var_selected <- tumor_bilat_no_is %>% select(  moddiag,  T,  N,  subtype,  concordance_subtype,
  dcis,  tumor_cellularity,  gradeclasse,
  it_til_perc,str_til_perc,  histo_3cl,  multifocal_bin,  tyct) %>% colnames()

names_var_selected <-  c("diagnostic modality","clinical T stage","clinical N stage","BC subtype","Concordance",
                          "DCIS component","tumor cellularity ","grade","IT TILs (%)","str TILs (%)",
                          "histological type","multifocality","Chemotherapy regimen")
data_dictionary    <- data.frame(var       = var_selected ,
                                 names_var = names_var_selected)
source('~/RT2Lab/databases/core/00_common/src/R_functions_Nadir/functions_RT2_Nadir.R', local = TRUE)

table1 <- table1_rt2(mydataset          = tumor_bilat_no_is,
                     raw_variables      = var_selected,
                     data_dictionary    = data_dictionary,
                     stratif            = "BRCA_mut",
                     # stratif_order      = c("Yes","No"),
                     missing            = FALSE,
                     perc_by_column     = FALSE,
                     apply_default_test = FALSE)

tableone_brca_tumors <- table1[[1]]
write.csv2(tableone_brca_tumors, file="/Users/ahamypet/RT2Lab/BC_BILAT_NEO/codes_git/TableS3bistableone_brca_tumors.csv")

# Legends for exact p-values 
tumor_bilat_no_is %>% select(BRCA_mut,subtype) %>% table() %>% chisq.test()
tumor_bilat_no_is %>% select(BRCA_mut,gradeclasse) %>% table() %>% chisq.test()

