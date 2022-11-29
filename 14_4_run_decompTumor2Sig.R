library(musicatk)
variant_table = read.table('bilat_mutations.csv', sep='\t', header=T)
new_variant_table = variant_table[,c('chr', 'start', 'end', 'reference', 'variant', 'sample')]
variants = extract_variants_from_matrix(new_variant_table)
colnames(new_variant_table) = c('Chromosome', 'Start_Position', 'End_Position', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode')
variants = extract_variants_from_matrix(new_variant_table)
head(variants)
g <- select_genome("hg19")
musica <- create_musica(x = variants, genome = g, convert_dbs = FALSE)
musica
build_standard_table(musica, g = g, table_name = "SBS96")
musica
data(cosmic_v3_sbs_sigs_exome)

pred_cosmic <- predict_exposure(musica = musica, table_name = "SBS96",
                               signature_res = cosmic_v3_sbs_sigs_exome,
                               signatures_to_use =  c('SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS9', 'SBS13', 'SBS17a', 'SBS17b', 'SBS18', 'SBS37', 'SBS40', 'SBS41'),
                               algorithm = "decompTumor2Sig")
norm_sig = prop.table(pred_cosmic@exposures, 2)
write.table(norm_sig, 'bilat_decompTumor2Sig.csv', sep='\t', quote=FALSE)