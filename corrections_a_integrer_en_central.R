# repreprocess_dates_delays_bilat.R
# => to run

# We consider relapse in all 9; with corresponding delays  
nip_relapse <- tumor_bilat_no_is %>% group_by (NUMDOS, rfs) %>% count() %>% ungroup() %>% filter(n < 2)  %>% select(NUMDOS) %>% as.matrix() %>% as.character()

tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% nip_relapse, c("NUMDOS","rfs","delrfs")]
# NUMDOS rfs   delrfs
# 63  0685406   0 76.78029
# 64  0685406   1 57.82341
# 197 0881266   1 56.70637
# 198 0881266   0 59.23614
# 217 0883466   0 54.96509
# 218 0883466   1 38.89938

tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "0685406" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 57.82341
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "0688544" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 102.73511
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "0881266" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 56.70637
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "0883466" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 38.89938
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "0883712" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 66.26694
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "0980461" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 77.93018
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "1104400" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 49.37988
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "1107440" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 54.73511
tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% "1117756" & tumor_bilat_no_is$rfs == 0, c("delrfs")] <- 41.36345

tumor_bilat_no_is[tumor_bilat_no_is$NUMDOS %in% nip_relapse, c("rfs")] <- 1

tumor_bilat_no_is %>% select (NUMDOS, rfs) %>% unique() %>% nrow() # 313 # All status  are OK


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
# When negative, consider as NA


