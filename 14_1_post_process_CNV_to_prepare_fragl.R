library(CNTools)

# the objective of this script is to create fragl tables for CNV


all_genes = read.table('mart_export_cytoband_hg19.txt', header=T, sep='\t', stringsAsFactors=F)
ok_chr = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X')
all_genes_restr = all_genes[all_genes$Chromosome.scaffold.name %in% ok_chr,]
colnames(all_genes_restr) = c('Gene.stable.ID', 'Gene.stable.ID.version', 'Transcript.stable.ID', 'Transcript.stable.ID.version', 'chrom', 'start', 'end', 'strand', 'bandid', 'gene_name', 'source_gene_name')
#cytoband$chrom = gsub('chr', '', cytoband$chrom)
sub_all_genes_restr = all_genes_restr[,c('Gene.stable.ID', 'chrom', 'start', 'end', 'strand', 'bandid', 'gene_name', 'source_gene_name')]
all_genes_restr_dedup = sub_all_genes_restr[!duplicated(sub_all_genes_restr), ]



p_folders = list.files('results/superFreq/plots', pattern='P*_[LR]')
for (patient_side in p_folders) {
    seg_list = list.files(paste('results/superFreq/plots', patient_side, 'data', sep='/'),
                          pattern='CNAsegments_*')
    for (filename in seg_list) {
        complete_name = paste('results/superFreq/plots', patient_side, 'data', filename, sep='/')
        cnv_table = read.table(complete_name, sep='\t', header=T, stringsAsFactors=F)
        cnv_table$call <- gsub('?', '', cnv_table$call, fixed=T)
        cnv_table$call <- gsub('CL', '', cnv_table$call)

        cnv_table$seg.mean = nchar(as.character(cnv_table$call))
        #cnv_table$id = paste(cnv_table$chr, cnv_table$start, sep='_')
        cnv_table$id = 'avg_tot_cn'
        sub_cnv_table = cnv_table[,c('chr', 'start', 'end', 'seg.mean', 'id')]
        seq = CNSeg(sub_cnv_table, chromosome="chr", end='end', start='start', segMean='seg.mean', id='id')
        rd_cyto_cn = getRS(seq, by='gene', imput = FALSE, XY = TRUE, geneMap=all_genes_restr_dedup, what='mean')
        head(rs(rd_cyto_cn))
        cnv_table$id = 'clonality'
        sub_cnv_table = cnv_table[,c('chr', 'start', 'end', 'clonality', 'id')]
        seq = CNSeg(sub_cnv_table, chromosome="chr", end='end', start='start', segMean='clonality', id='id')
        rd_cyto_ccf = getRS(seq, by='gene', imput = FALSE, XY = TRUE, geneMap=all_genes_restr_dedup, what='mean')
        head(rs(rd_cyto_ccf))
        df_cn = rs(rd_cyto_cn)
        df_ccf = rs(rd_cyto_ccf)
        df_cn$clonality = df_ccf$clonality
        output_file = paste('results/superFreq/plots', patient_side, 'data', gsub("CNAsegments", "pre_fragl", filename), sep='/')
        write.table(df_cn, output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
}
