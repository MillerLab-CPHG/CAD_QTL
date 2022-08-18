#################### Looking at colocalization using SMR and coloc, as well as paintor results
#################### For reference files used in this script, please see Themes_and_baselines.R

#################### For mixQTL analyses, see mixQTL_figs.R
#################### For local-ancestry-adjusted eQTL analyses, see LA_QTL_figs.R
#################### For differential expression or RNA bulk data plots, see RNAseq_expression_figs.R
#################### For generalization/replication, see "UVA_QTL_vs_GTEx_STARNET_figs.R"

setwd("/Users/ch2um/Desktop/QTL")

all_coloc <- fread("QTL_coloc_results.txt")
all_coloc <- all_coloc[gene_name %in% mixq_new[mix_bh<0.05, gene_name]]
coloc_genes <- all_coloc[coloc_evidence=="1", gene_name]
bp <- all_coloc[coloc_evidence=="1" | indy_evidence=="1"]

mixq_coloc <- fread("Annotated_mixQTL_results.txt", select=c("gene_name", "type2", "chr", "start", "meanExp", "mix_bh", "mix_log10p"))
setkey(all_coloc, gene_name)
setkey(mixq_coloc, gene_name)
forplot <- all_coloc[mixq_coloc, nomatch=0]

##### Heatmap coloc #####
forplot.m <- as.data.table(melt(forplot[mix_bh<0.05 & (BBJ_PPH4>0.8 | Matsunaga_PPH4>0.8 | hart2_PPH4>0.8 | 
                                vdh_PPH4>0.8 | pp_PPH4>0.8 | dbp_PPH4>0.8 | sbp_PPH4>0.8 | top_cac_PPH4>0.8 | imt_PPH4>0.8 | 
                                  CAC_1000G_PPH4>0.8 | CHARGE_plaque_PPH4>0.8), 
                          c("gene_name", "BBJ_PPH4", "Matsunaga_PPH4", "hart2_PPH4", "vdh_PPH4", "pp_PPH4", "dbp_PPH4", 
                          "sbp_PPH4", "top_cac_PPH4", "imt_PPH4", "CAC_1000G_PPH4", "CHARGE_plaque_PPH4")], 
                          id.vars="gene_name")) 
forplot.m <- forplot.m[!is.na(gene_name) & !is.na(value)]

gene_list <- as.data.table(c("gene_1","gene_2","gene_3","gene_4","gene_5","gene_5","gene_7","gene_8"))
setnames(gene_list, "V1", "gene_name")
gene_list[, n:=1:nrow(gene_list)]
setkey(gene_list, gene_name); setkey(forplot.m, gene_name)
qplot <- forplot.m[gene_list]
qplot[is.na(value), value:=0]
qplot <- qplot[!is.na(variable) & !is.na(gene_name)]
setorder(qplot, -n)

qplot[, gene_name:=reorder(gene_name, -n)]
qplot[variable=="BBJ_PPH4", study:="Koyama CAD"]; qplot[variable=="Matsunaga_PPH4", study:="Matsunaga CAD"]
qplot[variable=="hart2_PPH4", study:="Hartiala MI"]; qplot[variable=="vdh_PPH4", study:="van der Harst CAD"]
qplot[variable=="pp_PPH4", study:="Evangelou PP"]; qplot[variable=="dbp_PPH4", study:="Evangelou DBP"]
qplot[variable=="sbp_PPH4", study:="Evangelou SBP"]; qplot[variable=="top_cac_PPH4", study:="TOPMed CAC"]
qplot[variable=="imt_PPH4", study:="Franceschini IMT"]; qplot[variable=="CHARGE_plaque_PPH4", study:="Franceschini plaque"]
qplot[variable=="CAC_1000G_PPH4", study:="1000G CAC"]

heat <-  ggplot(qplot, aes(x=study, y=gene_name)) + 
  geom_tile(fill="forestgreen", alpha=qplot$value) + 
  scale_fill_identity() + scale_alpha_continuous() + th + 
  guides(y = guide_axis(n.dodge = 2)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust = 1, hjust=1, face="italic", family="serif", size=40),
    axis.text.x = element_text(angle = 90, hjust=1, vjust=1, face="bold", family="serif", size=40),
    plot.margin = margin(t = 25, r = 25, b = 25, l = 25))

pdf(file="Coloc_eqtl_heatmap.pdf", width=15, height=30)
heat
dev.off()


##### Making SMR coloc plot for NDUFV3 #####
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
# Adding PDE9A
                geom_segment(aes(x=42700000, xend=42775509, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=42775509, xend=42775509, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_text(aes(x=(42734000), y=(-1.5), label="PDE9A"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding WDR4
                geom_segment(aes(x=42843094, xend=42879531, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=42843094, xend=42843094, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_segment(aes(x=42879531, xend=42879531, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_text(aes(x=(42860000), y=(-1.5), label="WDR4"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding NDUFV3
                geom_segment(aes(x=42893268, xend=42913304, y=(-2.5), yend=(-2.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=42893268, xend=42893268, y=(-2.25), yend=(-2.75), color="#490291"), size=1.5) +
                geom_segment(aes(x=42913304, xend=42913304, y=(-2.25), yend=(-2.75), color="#490291"), size=1.5) +
                geom_text(aes(x=(42903300), y=(-3.5), label="NDUFV3"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding PKNOX1
                geom_segment(aes(x=42974562, xend=43033931, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=42974562, xend=42974562, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_segment(aes(x=43033931, xend=43033931, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_text(aes(x=(43003900), y=(-1.5), label="PKNOX1"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding CBS
                geom_segment(aes(x=43053880, xend=43072193, y=(-2.5), yend=(-2.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=43053880, xend=43053880, y=(-2.25), yend=(-2.75), color="#490291"), size=1.5) +
                geom_segment(aes(x=43072193, xend=43072193, y=(-2.25), yend=(-2.75), color="#490291"), size=1.5) +
                geom_text(aes(x=(43063190), y=(-3.5), label="CBS"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding U2AF1
                geom_segment(aes(x=43092956, xend=43107565, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=43092956, xend=43092956, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_segment(aes(x=43107565, xend=43107565, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_text(aes(x=(43100265), y=(-1.5), label="U2AF1"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding CRYAA
                geom_segment(aes(x=43169601, xend=43172294, y=(-2.5), yend=(-2.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=43169601, xend=43169601, y=(-2.25), yend=(-2.75), color="#490291"), size=1.5) +
                geom_segment(aes(x=43172294, xend=43172294, y=(-2.25), yend=(-2.75), color="#490291"), size=1.5) +
                geom_text(aes(x=(43171500), y=(-3.5), label="CRYAA"), color="#490291", family="serif", size=5.5, fontface=4) + 
# Adding SIK1
                geom_segment(aes(x=43414483, xend=43427131, y=(-0.5), yend=(-0.5), color="#490291"), size=1.5) +
                geom_segment(aes(x=43414483, xend=43414483, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_segment(aes(x=43427131, xend=43427131, y=(-0.75), yend=(-0.25), color="#490291"), size=1.5) +
                geom_text(aes(x=(43420500), y=(-1.5), label="SIK1"), color="#490291", family="serif", size=5.5, fontface=4) + 
# OK finally plotting GWAS/eQTL results
                geom_point(data=coloc_m.m[variable=="eQTL_p"], cex=6, shape=22, color="gray50", fill="maroon", alpha=0.75, ) +
                geom_point(data=coloc_m.m[variable=="DBP_p"], cex=6, shape=21, fill="gray80", color="gray20") + 
                geom_point(data=coloc_m.m[variable=="Hart2_p"], cex=7, shape=23, color="gray50", fill="green3") + 
                geom_point(data=coloc_m.m[variable=="Top_CAC_p"], cex=6, shape=24, fill="blue3", color="gray50", alpha=0.6) + 
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
  
pdf("Coloc_NDUFV3_regional_plot_July2022.pdf", h=8.5, w=11)
multi
dev.off()
