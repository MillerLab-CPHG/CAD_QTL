#################### Looking at colocalization using coloc, as well as paintor results
#################### For reference files used in this script, please see Themes_and_baselines.R

setwd("/Users/ch2um/Desktop/QTL")

all_coloc <- fread("QTL_coloc_results.txt")
all_coloc <- all_coloc[gene_name %in% mixq_new[mix_bh<0.05, gene_name]]
coloc_genes <- all_coloc[coloc_evidence=="1", gene_name]
bp <- all_coloc[coloc_evidence=="1" | indy_evidence=="1"]

mixq_coloc <- fread("Annotated_mixQTL_results.txt", select=c("gene_name", "type2", "chr", "start", "meanExp", "mix_bh", "mix_log10p"))
setkey(all_coloc, gene_name)
setkey(mixq_coloc, gene_name)
forplot <- all_coloc[mixq_coloc, nomatch=0]

##### Plotting coloc heatmap in 500 steps #####
forplot.m <- as.data.table(melt(forplot[mix_bh<0.05 & coloc_evidence==1, c("gene_name", "Koyama_CAD_PPH4", "Matsunaga_CAD_PPH4", 
                          "Hartiala_CAD_PPH4", "vdh_CAD_PPH4", "PP_PPH4", "DBP_PPH4", "SBP_PPH4", "LDL_PPH4", "HDL_PPH4", "TC_PPH4",
                          "logTG_PPH4", "CAC_1000G_PPH4", "CHARGE_IMT_PPH4", "CHARGE_Plaque_PPH4", "MVP_TRANS_PPH4","MVP_AFR_PPH4", 
                          "MVP_HIS_PPH4", "MVP_TRANS_PPH4")], 
                          id.vars="gene_name")) #"Topmed_CAC_PPH4", 
forplot.m <- forplot.m[!is.na(gene_name)]

#### complex heatmap plot for clustering algorithm only not for final figure #####
forplot.h <- as.data.frame(forplot[mix_bh<0.05 & gene_type=="protein_coding" & coloc_evidence==1, c("gene_name", "Koyama_CAD_PPH4", "Matsunaga_CAD_PPH4", 
                            "Hartiala_CAD_PPH4", "vdh_CAD_PPH4", "PP_PPH4", "DBP_PPH4", "SBP_PPH4", "LDL_PPH4", "HDL_PPH4", "TC_PPH4",
                            "logTG_PPH4", "CAC_1000G_PPH4", "CHARGE_IMT_PPH4", "CHARGE_Plaque_PPH4", "MVP_TRANS_PPH4", "MVP_AFR_PPH4", 
                            "MVP_HIS_PPH4", "MVP_EUR_meta_PPH4")]) #"Topmed_CAC_PPH4", 
rownames(forplot.h) <- forplot.h$gene_name
forplot.h$gene_name <- NULL
setnames(forplot.h, c("Koyama CAD", "Matsunaga CAD", "Hartiala CAD", "van der Harst CAD", "Evangelou PP", 
                      "Evangelou DBP", "Evangelou SBP", "Graham LDL", "Graham HDL", "Graham TC", "Graham logTG",
                      "1000G CAC", "CHARGE IMT", "CHARGE Plaque", "MVP All", "MVP AFR", "MVP HIS", "MVP EUR")) #"TOPMed CAC", 
mycols <- colorRamp2(breaks = c(0, 0.5, 1), colors = c("#FFFFFF", "#CCE5FF", "#000099"))
ht_opt(heatmap_column_names_gp = gpar(fontsize = 20))
heat <- Heatmap(forplot.h, name = "Colocalization PPH4", col=colorRamp2(breaks = c(0, 0.5, 1), colors = c("#FFFFFF", "#CCFFCC", "#006633")),
        column_title = " ", row_title = " ",
        row_names_gp = gpar(fontsize = 16), cluster_columns=F)

pdf(file="Coloc_heatmap_cluster_Nov2022.pdf", width=12, height=30)
heat
dev.off()

##### List of genes clustered heirarchically using complex heatmap above #####
gene_list <- as.data.table(c("ZNF100","TCF21","ARHGAP42","THSD4","FHL3","ADAMTS3","TBX20","BACH1","SKIV2L", 
                             "SREBF1","PAPPA","CCDC38","SNRPF","TRIM54","NPHP3","SFRP2","SMOC1","SNX31",
                             "SPATA20","HPSE2","POGLUT3","GSTM1","ZFAND2A","CAPG","SCMH1","CENPV",
                             "ZNF641","SLC15A2","STEAP1B","MET","IL5","STEAP2","MANSC4","CMTM3","ZFYVE19",
                             "POLL","QSOX1","ANKDD1B","BTBD16","PGPEP1","LSM4","HLA−C")) #RPL71, XRCC2
setnames(gene_list, "V1", "gene_name"); gene_list[, n:=1:nrow(gene_list)]

##### Single-cell RNAseq cell-type-specific expression #####
sc_names <- c("p_val", "avg_log2FC", "pct1", "pct2", "adj_p", "cell_type", "gene")

##### Alsaigh scRNAseq cell-type-specific expression
al <- fread("/path/Alsaigh_differential_markergenes_by_tabulus_sapien_reference.csv"); setnames(al, sc_names)
al <- al[gene %in% gene_list$gene_name & pct1>0.2 & adj_p<0.05]
setorder(al, adj_p); al <- al[!duplicated(gene)]
gene_list[gene_name %in% al[cell_type=="EC", gene], Alsaigh:="darkviolet"] 
gene_list[gene_name %in% al[cell_type=="FB", gene], Alsaigh:="skyblue"] 
gene_list[gene_name %in% al[cell_type=="Mø" | cell_type=="T Cell" | cell_type=="Mast Cell", gene], Alsaigh:="lightpink"] 
gene_list[gene_name %in% al[cell_type=="SMCs" | cell_type=="Pericyte Cell", gene], Alsaigh:="goldenrod3"] 

#### Hu scRNAseq cell-type-specific expression
hu <- fread("/path/Hu_differential_markergenes_by_tabulus_sapien_reference.csv"); setnames(hu, sc_names)
hu <- hu[gene %in% gene_list$gene_name & pct1>0.2 & adj_p<0.05]
setorder(hu, adj_p); hu <- hu[!duplicated(gene)]
gene_list[gene_name %in% hu[cell_type=="EC", gene], Hu:="darkviolet"] 
gene_list[gene_name %in% hu[cell_type=="FB", gene], Hu:="skyblue"] 
gene_list[gene_name %in% hu[cell_type=="Mø" | cell_type=="T Cell" | cell_type=="Mast Cell", gene], Hu:="lightpink"] 
gene_list[gene_name %in% hu[cell_type=="SMCs" | cell_type=="Pericyte Cell", gene], Hu:="goldenrod3"] 

#### Wirka scRNAseq cell-type-specific expression
wi <- fread("/path/Wirka_differential_markergenes_by_tabulus_sapien_reference.csv"); setnames(wi, sc_names)
wi <- wi[gene %in% gene_list$gene_name & pct1>0.2 & adj_p<0.05]
setorder(wi, adj_p); wi <- wi[!duplicated(gene)]
gene_list[gene_name %in% wi[cell_type=="EC", gene], Wirka:="darkviolet"] 
gene_list[gene_name %in% wi[cell_type=="FB", gene], Wirka:="skyblue"] 
gene_list[gene_name %in% wi[cell_type=="Mø" | cell_type=="T Cell" | cell_type=="Mast Cell", gene], Wirka:="lightpink"] 
gene_list[gene_name %in% wi[cell_type=="SMCs" | cell_type=="Pericyte Cell", gene], Wirka:="goldenrod3"] 

#### Making combined plot dataset #####
setkey(gene_list, gene_name); setkey(forplot.m, gene_name)
pplot <- gene_list[forplot.m]
pplot[is.na(value), value:=0]
pplot <- pplot[!is.na(n)]
studies <- as.data.table(c("Matsunaga_CAD_PPH4", "Koyama_CAD_PPH4", "Hartiala_CAD_PPH4", "vdh_CAD_PPH4", "MVP_TRANS_PPH4",  
                           "DBP_PPH4", "SBP_PPH4", "PP_PPH4", "HDL_PPH4", "LDL_PPH4", "logTG_PPH4", "TC_PPH4", "CAC_1000G_PPH4", 
                           "CHARGE_IMT_PPH4", "CHARGE_plaque_PPH4")) 
setkey(studies, V1); setkey(pplot, variable)
colocplot <- pplot[studies]
colocplot <- colocplot[!is.na(variable) & !is.na(gene_name)]
colocplot[, m:=1:nrow(colocplot)]

colocplot[, gene_name:=reorder(gene_name, -m)]
colocplot[variable=="Koyama_CAD_PPH4", study:="Koyama CAD"]; colocplot[variable=="Matsunaga_CAD_PPH4", study:="Matsunaga CAD"]
colocplot[variable=="Hartiala_CAD_PPH4", study:="Hartiala MI"]; colocplot[variable=="vdh_CAD_PPH4", study:="van der Harst CAD"]
colocplot[variable=="PP_PPH4", study:="Evangelou PP"]; colocplot[variable=="DBP_PPH4", study:="Evangelou DBP"]
colocplot[variable=="SBP_PPH4", study:="Evangelou SBP"]; colocplot[variable=="Topmed_CAC_PPH4", study:="TOPMed CAC"]
colocplot[variable=="CHARGE_IMT_PPH4", study:="Franceschini IMT"]; colocplot[variable=="CHARGE_plaque_PPH4", study:="Franceschini plaque"]
colocplot[variable=="CAC_1000G_PPH4", study:="Kavousi CAC"]; colocplot[variable=="logTG_PPH4", study:="Graham logTG"]
colocplot[variable=="LDL_PPH4", study:="Graham LDL"]; colocplot[variable=="HDL_PPH4", study:="Graham HDL"]
colocplot[variable=="TC_PPH4", study:="Graham TC"]; colocplot[variable=="MVP_TRANS_PPH4", study:="Tcheandjieu CAD"]

##### Finally making plot #####
colocplot[, gene_name:=reorder(gene_name, n)]

heat <-  ggplot(colocplot, aes(x=study, y=gene_name)) + plot_th + 
  geom_tile(fill="black", alpha=colocplot$value) + 
  scale_fill_identity() + scale_alpha_continuous() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(vjust=0.5, hjust=1, size=40, face="plain")
    axis.text.y = element_text(angle = 90, hjust=0.5, vjust=0.5, face="italic", family="sans", size=40),
    plot.margin = margin(t = 25, r = 25, b = 25, l = 25)) 

pdf(file="Coloc_eqtl_heatmap_Nov2022.pdf", width=16, height=26)
heat +  geom_tile(data=colocplot, aes(y=-1, x=gene_name, fill=Alsaigh)) + 
        geom_tile(data=colocplot, aes(y=-2, x=gene_name, fill=Hu)) + 
        geom_tile(data=colocplot, aes(y=-3, x=gene_name, fill=Wirka))
dev.off()


##### Making coloc regional plot for NDUFV3 ##### 
#(Didn't end up including this in paper but format was used for a couple other figures)
ndufv3_qtl <- fread("mixQTL_NDUFV3.txt", select=c(1:3,14,18))
setnames(ndufv3_qtl, c("gene", "gene_name","variant","eQTL_p","method.meta"))
hg38_freq <- fread("hg38_allele_freqs.txt", select=c(1:3))
setnames(hg38_freq, c("chr", "pos", "rsid"))
setkey(hg38_freq, rsid); setkey(ndufv3_qtl, variant)   
ndufv3_eqtl <- ndufv3_qtl[hg38_freq, nomatch=0]
                
ndufv3_dbp_gwas <- fread("DBP_chr21_GWAS.txt", select=c(1,3,4,2))
setnames(ndufv3_dbp_gwas, c("rsid", "chr", "pos", "DBP_p"))
ndufv3_dbp_gwas <- ndufv3_dbp_gwas[rsid %in% ndufv3_qtl$variant]
ndufv3_cad_gwas <- fread("Hartiala2_chr21_GWAS.txt", select=c(2:5))
setnames(ndufv3_cad_gwas, c("rsid", "chr", "pos", "Hart2_p"))
ndufv3_cad_gwas <- ndufv3_cad_gwas[rsid %in% ndufv3_qtl$variant]
ndufv3_cac_gwas <- fread("Topmed_CAC_chr21_GWAS.txt", select=c(2:5))
setnames(ndufv3_cac_gwas, c("rsid", "chr", "pos", "Top_CAC_p"))
ndufv3_cac_gwas <- ndufv3_cac_gwas[rsid %in% ndufv3_qtl$variant]
dfs <- list(ndufv3_dbp_gwas, ndufv3_cad_gwas, ndufv3_cac_gwas)
smr_combo <- dfs %>% reduce(full_join, by=c('rsid', 'chr', 'pos'))

setkey(coloc_combo, rsid, chr, pos)
setkey(ndufv3_eqtl, variant, chr, pos)
coloc_multi <- coloc_combo[ndufv3_eqtl]

coloc_m.m <- as.data.table(melt(coloc_multi[, c(1:6,8,9)], id.vars=c("rsid", "chr", "pos", "gene_name")))
coloc_m.m <- coloc_m.m[!is.na(value) & value<0.05 & pos>42700000]

multi <- ggplot(data=coloc_m.m, aes(x=pos, y=-log10(value), color=variable)) + th + 
                labs(x="Position (Chromosome 21)", y=expression("log10(p-value)")) + scale_color_identity() + 
# Adding PDE9A, WDR4, NDUFV3, PKNOX1, CBS, U2AF1, CRYAA, SIK1
                geom_segment(aes(x=42700000, xend=42775509, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(42734000), y=(-1.5), label="PDE9A"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=42843094, xend=42879531, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(42860000), y=(-1.5), label="WDR4"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=42893268, xend=42913304, y=(-2.5), yend=(-2.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(42903300), y=(-3.5), label="NDUFV3"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=42974562, xend=43033931, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(43003900), y=(-1.5), label="PKNOX1"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=43053880, xend=43072193, y=(-2.5), yend=(-2.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(43063190), y=(-3.5), label="CBS"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=43092956, xend=43107565, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(43100265), y=(-1.5), label="U2AF1"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=43169601, xend=43172294, y=(-2.5), yend=(-2.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(43420500), y=(-3.5), label="CRYAA"), color="#490291", family="sans", size=5.5, fontface=4) + 
                geom_segment(aes(x=43427131, xend=43414483, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5, arrow=arrow(length=unit(0.5, "cm"))) +
                geom_text(aes(x=(43420500), y=(-1.5), label="SIK1"), color="#490291", family="sans", size=5.5, fontface=4) + 
# Plotting GWAS/eQTL results
                geom_point(data=coloc_m.m[variable=="eQTL_p"], cex=6, shape=22, color="gray50", fill="maroon", alpha=0.75, ) +
                geom_point(data=coloc_m.m[variable=="DBP_p"], cex=6, shape=21, color="gray20", fill="gray80") + 
                geom_point(data=coloc_m.m[variable=="Hart2_p"], cex=7, shape=23, color="gray50", fill="green3") + 
                geom_point(data=coloc_m.m[variable=="Top_CAC_p"], cex=6, shape=24, color="gray50", fill="blue3", alpha=0.6) + 
                geom_label_repel(data=coloc_m.m[variable=="DBP_p" & rsid=="rs12627514"], aes(label="rs12627514, SMRp = 4.6e-7"), 
                                fill="white", colour="gray20", size=7, family = "serif", segment.alpha = 0.75, 
                                label.padding=0.5, box.padding=0.3, nudge_y=(-1.5), nudge_x=(-200000)) + 
                geom_label_repel(data=coloc_m.m[variable=="Hart2_p" & rsid=="rs659057"], aes(label="rs659057, SMRp = 0.006"), 
                                fill="white", colour="green4", size=7, family = "serif", segment.alpha = 0.75, 
                                label.padding=0.5, box.padding=0.3, nudge_y=10, nudge_x=(-100000)) + 
                geom_label_repel(data=coloc_m.m[variable=="eQTL_p" & rsid=="rs4148974"], aes(label=rsid), 
                                fill="white", colour="maroon", size=7, family = "serif", segment.alpha = 0.75, 
                                label.padding=0.5, box.padding=0.3, nudge_y=3) + 
                geom_label_repel(data=coloc_m.m[variable=="Top_CAC_p" & rsid=="rs17115411"], aes(label="rs17115411, SMRp = 0.62"), 
                                fill="white", colour="blue3", size=7, family = "serif", segment.alpha = 0.75, 
                                label.padding=0.5, box.padding=0.3, nudge_y=5, nudge_x=(-30000)) 
  
pdf("Coloc_NDUFV3_regional_plot_2022.pdf", h=8.5, w=11)
multi
dev.off()
