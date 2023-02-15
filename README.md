# CAD_QTL

The set of scripts in this repo are intended to generate the necessary input files and run programs to identify expression QTLs for RNA sequencing data using two programs: mixQTL (https://github.com/liangyy/mixqtl/wiki) and a QTLtools-based (https://qtltools.github.io/qtltools/) pipeline to incorporate local ancestry into the model. If you use these scripts for your own work, please cite our manuscript: https://www.medrxiv.org/content/10.1101/2023.02.09.23285622v1

This is a brief overview of our study design, which involved two methods (inferred local ancestry adjustment with permuted effect estimates generated with QTLtools and a combined haplotype-specific expression and total read count regression estimates generated with mixQTL) to identify expression QTLs in human coronary artery, as well as a separate analysis of splicing QTLs also performed with QTLtools. QTL identification was followed by various fine-mapping and annotation efforts utilizing the scripts in the respectively named folders. Links to tools used should be included in all scripts.

![QTL_cropped_graphical_abstract](https://user-images.githubusercontent.com/48288433/219075077-1a9a251d-335e-4648-b9fa-8b74ad09caac.png)

For questions regarding scripts or implementation, please email Chani at ch2um at virginia dot edu. 
