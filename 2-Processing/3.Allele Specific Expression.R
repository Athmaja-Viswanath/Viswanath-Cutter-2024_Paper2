###Allele Specific Expression Analysis
library(DESeq2)
library(edgeR)
library(ggplot2)
library(cowplot) #add on to ggplot for better themes and figure customization
library(lemon) #to work with legends and axes in ggplot2
library(dplyr)
library(tidyverse)
library(gdata)
library(RColorBrewer)
library(colorBlindness)
library(colorspace)

theme_set(theme_classic())
setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT
orthologs = read.table("new_1to1_orthologgenelist.txt", sep="\t", head=T, comment.char="#")
head(orthologs)
cre_ortho = orthologs[, 1]
clat_ortho = orthologs[ ,2]
cumulative_ortho = orthologs[,3]

nrow(orthologs)

###H1
##DATA PREP
H1_ase= read.table("H1_ase_counts.txt", sep="\t", head=T, row.name=1, comment.char = "#")
H1_ase = H1_ase[cumulative_ortho, ]
names(H1_ase)

H2_ase = read.table("H2_ase_counts.txt", sep="\t", head=T, row.name=1, comment.char = "#")
H2_ase = H2_ase[cumulative_ortho, ]

nrow(H2_ase)
names(H2_ase)

#C.remanei
cre_counts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_counts_ortho = cre_counts[cre_ortho, ]
rownames(cre_counts_ortho) = cumulative_ortho #changing row names
nrow(cre_counts_ortho)

#C.latens
clat_counts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_counts_ortho = clat_counts[clat_ortho, ]
rownames(clat_counts_ortho) = cumulative_ortho
nrow(clat_counts_ortho)


#Combining all data

counts_orthologs = cbind(cre_counts_ortho, clat_counts_ortho, H1_ase[,seq(1,24,2)] + H1_ase[,seq(2,24,2)], H2_ase[,seq(1,30,2)] + H2_ase[,seq(2,30,2)])
head(H1_ase[,seq(1,24,2)])
head(H1_ase[,seq(2,24,2)])
head(H1_ase[,seq(1,24,2)] + H1_ase[,seq(2,24,2)])
head(H2_ase[,seq(1,30,2)])

##coldata Prep
tissue_type.g = c(rep("G", 6)) 
tissue_type.s = c(rep("S", 6))
tissue_type.gs = c(rep("G", 3), rep("S", 3))

Sex.fm = c(rep("F", 3), rep("M", 3))
Sex.m = c(rep("M", 6))
Sex.f = c(rep("F", 6))

batch_wt = c("1", "2", "3", "1", "2", "3")

species_wt = c(rep(c("Cre","Clat"), each=3))

###FG

#############################################################
##Parental expression divergence (copied from other script)##
#############################################################
##This DEG analysis is done in "Differential gene expression analysis" RScript
###C.remanei x C.latens###########

#FG x FG
head(counts_orthologs)
fg.cre.clat = counts_orthologs[ , c(1,2,3,16,17,18)]
colnames(fg.cre.clat)
head(fg.cre.clat)
samples_fg.cre.clat = colnames(fg.cre.clat)
nrow(fg.cre.clat)
 
coldata_fg.cre.clat = data.frame(samples_fg.cre.clat, Sex.f, tissue_type.g, batch_wt, species_wt)
 
dds_fg.cre.clat = DESeqDataSetFromMatrix(countData = fg.cre.clat, colData = coldata_fg.cre.clat, design = ~batch_wt + species_wt)
dds_fg.cre.clat = DESeq(dds_fg.cre.clat)
dds.fg.cre.clat= results(dds_fg.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fg.cre.clat, alpha = 0.05))
dds.fg.cre.clat.res = cbind(as.data.frame(dds.fg.cre.clat),as.data.frame(dds.fg.cre.clat) %>%
                         mutate(fg.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(fg.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


dds.fg.cre.clat.res = dds.fg.cre.clat.res[order(row.names(dds.fg.cre.clat.res)), ]
#View(dds.fg.cre.clat.res)


##FS X FS

colnames(counts_orthologs)
fs.cre.clat = counts_orthologs[ , c(4,5,6,19,20,21)]
colnames(fs.cre.clat)

samples_fs.cre.clat = colnames(fs.cre.clat)
nrow(fs.cre.clat)

coldata_fs.cre.clat = data.frame(samples_fs.cre.clat, Sex.f, tissue_type.s, batch_wt, species_wt)
dds_fs.cre.clat = DESeqDataSetFromMatrix(countData = fs.cre.clat, colData = coldata_fs.cre.clat, design = ~batch_wt + species_wt)
dds_fs.cre.clat = DESeq(dds_fs.cre.clat)
dds.fs.cre.clat = results(dds_fs.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fs.cre.clat, alpha = 0.05))
dds.fs.cre.clat.res = cbind(as.data.frame(dds.fs.cre.clat),as.data.frame(dds.fs.cre.clat) %>%
                              mutate(fs.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                              dplyr::select(fs.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


dds.fs.cre.clat.res = dds.fs.cre.clat.res[order(row.names(dds.fs.cre.clat.res)), ]
#View(dds.fs.cre.clat.res)

##MG X MG
colnames(counts_orthologs)
mg.cre.clat = counts_orthologs[ , c(7,8,9,22,23,24)]
colnames(mg.cre.clat)

samples_mg.cre.clat = colnames(mg.cre.clat)
nrow(mg.cre.clat)

coldata_mg.cre.clat = data.frame(samples_mg.cre.clat, Sex.m, tissue_type.g, batch_wt, species_wt)

# dds_mg.cre.clat = DESeqDataSetFromMatrix(countData = mg.cre.clat, colData = coldata_mg.cre.clat, design = ~ batch + species)
# dds_mg.cre.clat = estimateSizeFactors(dds_mg.cre.clat)
# 
# dds_mg.cre.clat = rowSums(counts(dds_mg.cre.clat, normalized=TRUE) >= 10 ) >= 3
# summary(dds_mg.cre.clat) #11732 genes
# 
# mg.cre.clat_filtered.mg = mg.cre.clat[dds_mg.cre.clat, ]
# nrow(mg.cre.clat_filtered.mg)

dds_mg.cre.clat = DESeqDataSetFromMatrix(countData = mg.cre.clat, colData = coldata_mg.cre.clat, design = ~batch_wt + species_wt)
dds_mg.cre.clat = DESeq(dds_mg.cre.clat)
dds.mg.cre.clat= results(dds_mg.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_mg.cre.clat, alpha = 0.05))
dds.mg.cre.clat.res = cbind(as.data.frame(dds.mg.cre.clat),as.data.frame(dds.mg.cre.clat) %>%
                          mutate(mg.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                          dplyr::select(mg.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.mg.cre.clat.res = dds.mg.cre.clat.res[order(row.names(dds.mg.cre.clat.res)), ]

#View(dds.mg.cre.clat.res)


###MS X MS
colnames(counts_orthologs)
ms.cre.clat = counts_orthologs[ , c(10,11,12,25,26,27)]
colnames(ms.cre.clat)

samples_ms.cre.clat = colnames(ms.cre.clat)
nrow(ms.cre.clat)

coldata_ms.cre.clat = data.frame(samples_ms.cre.clat, Sex.m, tissue_type.s, batch_wt, species_wt)

# dds_ms.cre.clat = DESeqDataSetFromMatrix(countData = ms.cre.clat, colData = coldata_ms.cre.clat, design = ~ batch + species)
# dds_ms.cre.clat = estimateSizeFactors(dds_ms.cre.clat)
# 
# dds_ms.cre.clat = rowSums(counts(dds_ms.cre.clat, normalized=TRUE) >= 10 ) >= 3
# summary(dds_ms.cre.clat) #12720 genes
# 
# ms.cre.clat_filtered.ms = ms.cre.clat[dds_ms.cre.clat, ]
# nrow(ms.cre.clat_filtered.ms)

dds_ms.cre.clat = DESeqDataSetFromMatrix(countData = ms.cre.clat, colData = coldata_ms.cre.clat, design = ~batch_wt + species_wt)
dds_ms.cre.clat = DESeq(dds_ms.cre.clat)
dds.ms.cre.clat= results(dds_ms.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_ms.cre.clat, alpha = 0.05))
dds.ms.cre.clat.res = cbind(as.data.frame(dds.ms.cre.clat),as.data.frame(dds.ms.cre.clat) %>%
                          mutate(ms.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                          dplyr::select(ms.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.ms.cre.clat.res = dds.ms.cre.clat.res[order(row.names(dds.ms.cre.clat.res)), ]
#View(dds.ms.cre.clat.res)

##WM X WM

colnames(counts_orthologs)
wm.cre.clat = counts_orthologs[ , c(13,14,15,28,29,30)]
colnames(wm.cre.clat)

samples_wm.cre.clat = colnames(wm.cre.clat)
nrow(wm.cre.clat)

coldata_wm.cre.clat = data.frame(samples_wm.cre.clat, Sex.m, tissue_type.s, batch_wt, species_wt)

# dds_wm.cre.clat = DESeqDataSetFromMatrix(countData = wm.cre.clat, colData = coldata_wm.cre.clat, design = ~ batch + species)
# dds_wm.cre.clat = estimateSizeFactors(dds_wm.cre.clat)
# 
# dds_wm.cre.clat = rowSums(counts(dds_wm.cre.clat, normalized=TRUE) >= 10 ) >= 3
# summary(dds_wm.cre.clat) #12932 genes
# 
# wm.cre.clat_filtered.wm = wm.cre.clat[dds_wm.cre.clat, ]
# nrow(wm.cre.clat_filtered.wm)

dds_wm.cre.clat = DESeqDataSetFromMatrix(countData = wm.cre.clat, colData = coldata_wm.cre.clat, design = ~batch_wt + species_wt)
dds_wm.cre.clat = DESeq(dds_wm.cre.clat)
dds.wm.cre.clat= results(dds_wm.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_wm.cre.clat, alpha = 0.05))
dds.wm.cre.clat.res = cbind(as.data.frame(dds.wm.cre.clat),as.data.frame(dds.wm.cre.clat) %>%
                          mutate(wm.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                          dplyr::select(wm.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.wm.cre.clat.res = dds.wm.cre.clat.res[order(row.names(dds.wm.cre.clat.res)), ]
#View(dds.wm.cre.clat.res)
table(dds.wm.cre.clat.res$wm.cre.clat)

###########################################
##Allele specific divergence in H1#########
###########################################

Hybrid_allele_counts = cbind(H1_ase[,seq(1,24,2)], H1_ase[,seq(2,24,2)], H2_ase[,seq(1,30,2)], H2_ase[,seq(2,30,2)])

colnames(Hybrid_allele_counts)
colnames(Hybrid_allele_counts) = gsub("X", "", colnames(Hybrid_allele_counts)) ##replacing the last

##Coldata
sample_name = colnames(Hybrid_allele_counts)
tissue = substr(sample_name, 3,3) 
tissue = gsub("M", "W", tissue)
tissue[7] = "W"
tissue[8] = "W"
tissue[9] = "W"
tissue[19] = "W"
tissue[20] = "W"
tissue[21] = "W"
sex = substr(sample_name, 2,2)
sex[7] = "W"
sex[8] = "W"
sex[9] = "W"
sex[19] = "W"
sex[20] = "W"
sex[21] = "W"
species =  c(rep("Cre", 12), rep("Clat", 12),rep("Cre", 15), rep("Clat", 15))
Hybrid = c(rep("H1", 12), rep("H1", 12),rep("H2", 15), rep("H2", 15))
batch = c(rep(c(1, 2 , 3), 18))
coldata_h = data.frame(sample_name, species, Hybrid, sex, tissue, batch)
coldata_h


########Function for getting Trans effects
getTransEffects = function(x){
  x = as.data.frame(x)
  x$T = rep(c("P","H"), each=6)
  x$S = rep(rep(c("Clat","Cre"), each=3), 2)
  fit = lm(x ~ T/S, data=x) #lm(y ~ x/z, data) is just a shortcut for lm(y ~ x + x:z, data)
  fit2 = car::linearHypothesis(fit, c("TH:SCre = TP:SCre"))
  p.value = fit2$`Pr(>F)`[2]
  list(x, fit, fit2, p.value)
}

##Function to classify allele specific expression categories
# classify allele specific expression
# x = "parents hybrids_ASE trans_effects"
get_regulation_type = function(x){
  regulation_type = rep(NA, dim(x)[1])
  regulation_type[ x[,"p.value.sp"] > 0.05 & x[,"p.value.ase"] > 0.05 ] = "conserved"
  regulation_type[ is.na(regulation_type) & x[,"trans_effect"] > 0.05 & x[,"p.value.ase"] < 0.05 & x[,"p.value.sp"] < 0.05 ] = "cis-only"
  regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & x[,"p.value.ase"] > 0.05 & x[,"p.value.sp"] < 0.05 ] = "trans-only"
  regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & x[,"p.value.ase"] < 0.05 & x[,"p.value.sp"] > 0.05 ] = "cis-trans (compensatory)"
  regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & ((x[,"logFC.ase"] > 0 & x[,"logFC.sp"] < 0) | (x[,"logFC.ase"] < 0 & x[,"logFC.sp"] > 0)) ] = "cis x trans (compensatory)" ##check to see if this is compensatory
  regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & ((x[,"logFC.ase"] > 0 & x[,"logFC.sp"] > 0) | (x[,"logFC.ase"] < 0 & x[,"logFC.sp"] < 0)) ] = "cis + trans (enhancing)"
  regulation_type[ is.na(regulation_type) ] = "ambiguous"
  return(regulation_type)
}



#FG
keep = Hybrid == "H1" & tissue == "G"
dds_h1.fg = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h1.fg = DESeq(dds_h1.fg)
dds.h1.fg = results(dds_h1.fg, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h1.fg.res = cbind(as.data.frame(dds.h1.fg),as.data.frame(dds.h1.fg) %>%
                            mutate(h1.fg = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                            dplyr::select(h1.fg)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
#View(dds.h1.fg.res)


##Combining parental and H1 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.fg.h1 = cbindX(cpm(counts(dds_fg.cre.clat), log=T), cpm(counts(dds_h1.fg), log=T)) #these are objects before the result object from Differential gene expression analysis Rscript
df.cpm.fg.h1 = na.omit(df.cpm.fg.h1)

trans_effect_fg.h1 = lapply(split(df.cpm.fg.h1, 1:dim(df.cpm.fg.h1)[1]), getTransEffects) #each gene/row is as an element of the list, lapply says to pply this funciton to each element/gene/row
trans_p.values.fdr.fg.h1 = p.adjust(sapply(trans_effect_fg.h1, function(x) x[[4]] ), method="BH")

###creating a new vector with NAs for all the genes
tmp = rep(NA, dim(df.cpm.fg.h1)[1])
#changing the column names to genes names in the count file
names(tmp) = rownames(df.cpm.fg.h1)
##changing the column names in fdr file to rownames in count file
names(trans_p.values.fdr.fg.h1) = rownames(df.cpm.fg.h1)
#for the genes that have proper fdr p value, replace NAs for those p-values
tmp[ names(trans_p.values.fdr.fg.h1) ] = trans_p.values.fdr.fg.h1
###retain the filtered file as the same fdr variable
trans_p.values.fdr.fg.h1 = tmp



genes=sort(rownames(dds.fg.cre.clat))

dds.fg.cre.clat.res = dds.fg.cre.clat.res[genes,]
dds.h1.fg.res = dds.h1.fg.res[genes, ] ###this makes the order of gene names same in both files
rownames(dds.h1.fg.res) = genes
trans_p.values.fdr.fg.h1 = trans_p.values.fdr.fg.h1[genes] ##this also makes the gene order same
names(trans_p.values.fdr.fg.h1)= genes

fg_h1_ase = data.frame(logFC.sp=dds.fg.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h1.fg.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h1.fg.res$h1.fg) ],
                       species=c("Clat","Cre")[ factor(dds.h1.fg.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.fg.cre.clat.res$padj,
                       p.value.ase=dds.h1.fg.res$padj,
                       genes=rownames(dds.fg.cre.clat.res),
                       trans_effect=trans_p.values.fdr.fg.h1, stringsAsFactors=F)

View(fg_h1_ase)
head(dds.fg.cre.clat)

type.levels = c("conserved","ambiguous","trans-only","cis-only", "cis + trans (enhancing)", "cis x trans (compensatory)", "cis-trans (compensatory)")

fg_h1_ase$type = get_regulation_type(fg_h1_ase)

fg_h1_table = table(fg_h1_ase$type) ##summary table of # of genes in each category

##Filtering genes to be on differnt chromosmes and Labelling autosomes vs X-linked genes

orthologs_chr = read.table(file = "orthologs_chr.txt", sep = "\t", header = TRUE)
orthologs_chr = orthologs_chr[order(orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name), ]

fg_h1_ase_chr = rownames_to_column(fg_h1_ase)
rownames(fg_h1_ase_chr) = fg_h1_ase_chr$rowname

fg_h1_ase_chr = na.omit(cbind(fg_h1_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

fg_h1_ase_chr$A_X = "Autosomes"

fg_h1_ase_chr[(fg_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

##separating information for autosomes and x-chromosome
fg_h1_ase_autosome = fg_h1_ase_chr %>% filter(fg_h1_ase_chr$A_X == "Autosomes")
fg_h1_ase_X = fg_h1_ase_chr %>% filter(fg_h1_ase_chr$A_X == "X")

#write.csv(fg_h1_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fg_h1_ase_chr.csv")


fg_h1_ase_chr_count = as.data.frame(table(fg_h1_ase_chr$type, fg_h1_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fg_h1_ase_chr_count$Sample = "h1"
fg_h1_ase_chr_count$Tissue = "G"
fg_h1_ase_chr_count$Sex = "Female"
fg_h1_ase_chr_count$Name = paste(fg_h1_ase_chr_count$Sample, fg_h1_ase_chr_count$Sex, fg_h1_ase_chr_count$Tissue, fg_h1_ase_chr_count$Chromosome)

#autosomes vs X
fg_h1_a_x = fg_h1_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

##########################################################################################

#H2
keep = Hybrid == "H2" & tissue == "G" & sex == "F"
dds_h2.fg = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h2.fg = DESeq(dds_h2.fg)
dds.h2.fg = results(dds_h2.fg, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h2.fg.res = cbind(as.data.frame(dds.h2.fg),as.data.frame(dds.h2.fg) %>%
                        mutate(h2.fg = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h2.fg)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and H2 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.fg.h2 = cbindX(cpm(counts(dds_fg.cre.clat), log=T), cpm(counts(dds_h2.fg), log=T)) #these are objects before the result object
df.cpm.fg.h2 = na.omit(df.cpm.fg.h2)

trans_effect_fg.h2 = lapply(split(df.cpm.fg.h2, 1:dim(df.cpm.fg.h2)[1]), getTransEffects)
trans_p.values.fdr.fg.h2 = p.adjust(sapply(trans_effect_fg.h2, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.fg.h2)[1])
names(tmp) = rownames(df.cpm.fg.h2)
names(trans_p.values.fdr.fg.h2) = rownames(df.cpm.fg.h2)
tmp[ names(trans_p.values.fdr.fg.h2) ] = trans_p.values.fdr.fg.h2
trans_p.values.fdr.fg.h2 = tmp

genes=sort(rownames(dds.fg.cre.clat))
dds.fg.cre.clat.res = dds.fg.cre.clat.res[genes,]
dds.h2.fg.res = dds.h2.fg.res[genes, ]
rownames(dds.h2.fg.res) = genes
trans_p.values.fdr.fg.h2 = trans_p.values.fdr.fg.h2[genes]
names(trans_p.values.fdr.fg.h2)= genes

fg_h2_ase = data.frame(logFC.sp=dds.fg.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h2.fg.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h2.fg.res$h2.fg) ],
                       species=c("Clat","Cre")[ factor(dds.h2.fg.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.fg.cre.clat.res$padj,
                       p.value.ase=dds.h2.fg.res$padj,
                       genes=rownames(dds.fg.cre.clat.res),
                       trans_effect=trans_p.values.fdr.fg.h2, stringsAsFactors=F)

# classify allele specific expression
type.levels = c("conserved","ambiguous","trans-only","cis-only", "cis + trans (enhancing)", "cis x trans (compensatory)", "cis-trans (compensatory)")

fg_h2_ase$type = get_regulation_type(fg_h2_ase)

fg_h2_table = table(fg_h2_ase$type) ##summary table of # of genes in each category


#####Filtering and labelling chromosomes

fg_h2_ase_chr = rownames_to_column(fg_h2_ase)
rownames(fg_h2_ase_chr) = fg_h2_ase_chr$rowname

fg_h2_ase_chr = na.omit(cbind(fg_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

fg_h2_ase_chr$A_X = "Autosomes"

fg_h2_ase_chr[(fg_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

##separating information for autosomes and x-chromosome
fg_h2_ase_autosome = fg_h2_ase_chr %>% filter(fg_h2_ase_chr$A_X == "Autosomes")
fg_h2_ase_X = fg_h2_ase_chr %>% filter(fg_h2_ase_chr$A_X == "X")

#write.csv(fg_h2_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fg_h2_ase_chr.csv")


fg_h2_ase_chr_count = as.data.frame(table(fg_h2_ase_chr$type, fg_h2_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fg_h2_ase_chr_count$Sample = "h2"
fg_h2_ase_chr_count$Tissue = "G"
fg_h2_ase_chr_count$Sex = "Female"
fg_h2_ase_chr_count$Name = paste(fg_h2_ase_chr_count$Sample, fg_h2_ase_chr_count$Sex, fg_h2_ase_chr_count$Tissue, fg_h2_ase_chr_count$Chromosome)

#autosomes vs X
fg_h2_a_x = fg_h2_ase_chr %>% group_by(A_X, type) %>% dplyr::count()


################################################################

###FS

#H1
keep = Hybrid == "H1" & tissue == "S" & sex == "F"
dds_h1.fs = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h1.fs = DESeq(dds_h1.fs)
dds.h1.fs = results(dds_h1.fs, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h1.fs.res = cbind(as.data.frame(dds.h1.fs),as.data.frame(dds.h1.fs) %>%
                        mutate(h1.fs = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h1.fs)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and H1 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.fs.h1 = cbindX(cpm(counts(dds_fs.cre.clat), log=T), cpm(counts(dds_h1.fs), log=T)) #these are objects before the result object
df.cpm.fs.h1 = na.omit(df.cpm.fs.h1)

trans_effect_fs.h1 = lapply(split(df.cpm.fs.h1, 1:dim(df.cpm.fs.h1)[1]), getTransEffects)
trans_p.values.fdr.fs.h1 = p.adjust(sapply(trans_effect_fs.h1, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.fs.h1)[1])
names(tmp) = rownames(df.cpm.fs.h1)
names(trans_p.values.fdr.fs.h1) = rownames(df.cpm.fs.h1)
tmp[ names(trans_p.values.fdr.fs.h1) ] = trans_p.values.fdr.fs.h1
trans_p.values.fdr.fs.h1 = tmp

genes=sort(rownames(dds.fs.cre.clat))
dds.fs.cre.clat.res = dds.fs.cre.clat.res[genes,]
dds.h1.fs.res = dds.h1.fs.res[genes, ]
rownames(dds.h1.fs.res) = genes
trans_p.values.fdr.fs.h1 = trans_p.values.fdr.fs.h1[genes]
names(trans_p.values.fdr.fs.h1)= genes

fs_h1_ase = data.frame(logFC.sp=dds.fs.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h1.fs.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h1.fs.res$h1.fs) ],
                       species=c("Clat","Cre")[ factor(dds.h1.fs.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.fs.cre.clat.res$padj,
                       p.value.ase=dds.h1.fs.res$padj,
                       genes=rownames(dds.fs.cre.clat.res),
                       trans_effect=trans_p.values.fdr.fs.h1, stringsAsFactors=F)
#Regulatory divergence
type.levels = c("conserved","ambiguous","trans-only","cis-only", "cis + trans (enhancing)", "cis x trans (compensatory)", "cis-trans (compensatory)")

fs_h1_ase$type = get_regulation_type(fs_h1_ase)

fs_h1_table = table(fs_h1_ase$type) ##summary table of # of genes in each category


#####Filtering and labelling chromosomal information

fs_h1_ase_chr = rownames_to_column(fs_h1_ase)
rownames(fs_h1_ase_chr) = fs_h1_ase_chr$rowname

fs_h1_ase_chr = na.omit(cbind(fs_h1_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

fs_h1_ase_chr$A_X = "Autosomes"

fs_h1_ase_chr[(fs_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

##separating information for autosomes and x-chromosome
fs_h1_ase_autosome = fs_h1_ase_chr %>% filter(fs_h1_ase_chr$A_X == "Autosomes")
fs_h1_ase_X = fs_h1_ase_chr %>% filter(fs_h1_ase_chr$A_X == "X")

#write.csv(fs_h1_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fs_h1_ase_chr.csv")


fs_h1_ase_chr_count = as.data.frame(table(fs_h1_ase_chr$type, fs_h1_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fs_h1_ase_chr_count$Sample = "h1"
fs_h1_ase_chr_count$Tissue = "S"
fs_h1_ase_chr_count$Sex = "Female"
fs_h1_ase_chr_count$Name = paste(fs_h1_ase_chr_count$Sample, fs_h1_ase_chr_count$Sex, fs_h1_ase_chr_count$Tissue, fs_h1_ase_chr_count$Chromosome)

#autosomes vs X
fs_h1_a_x = fs_h1_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

###################################################################
##FS
#H2
keep = Hybrid == "H2" & tissue == "S" & sex == "F"
dds_h2.fs = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h2.fs = DESeq(dds_h2.fs)
dds.h2.fs = results(dds_h2.fs, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h2.fs.res = cbind(as.data.frame(dds.h2.fs),as.data.frame(dds.h2.fs) %>%
                        mutate(h2.fs = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h2.fs)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and H2 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.fs.h2 = cbindX(cpm(counts(dds_fs.cre.clat), log=T), cpm(counts(dds_h2.fs), log=T)) #these are objects before the result object
df.cpm.fs.h2 = na.omit(df.cpm.fs.h2)

trans_effect_fs.h2 = lapply(split(df.cpm.fs.h2, 1:dim(df.cpm.fs.h2)[1]), getTransEffects)
trans_p.values.fdr.fs.h2 = p.adjust(sapply(trans_effect_fs.h2, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.fs.h2)[1])
names(tmp) = rownames(df.cpm.fs.h2)
names(trans_p.values.fdr.fs.h2) = rownames(df.cpm.fs.h2)
tmp[ names(trans_p.values.fdr.fs.h2) ] = trans_p.values.fdr.fs.h2
trans_p.values.fdr.fs.h2 = tmp

genes=sort(rownames(dds.fs.cre.clat))
dds.fs.cre.clat.res = dds.fs.cre.clat.res[genes,]
dds.h2.fs.res = dds.h2.fs.res[genes, ]
rownames(dds.h2.fs.res) = genes
trans_p.values.fdr.fs.h2 = trans_p.values.fdr.fs.h2[genes]
names(trans_p.values.fdr.fs.h2)= genes

fs_h2_ase = data.frame(logFC.sp=dds.fs.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h2.fs.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h2.fs.res$h2.fs) ],
                       species=c("Clat","Cre")[ factor(dds.h2.fs.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.fs.cre.clat.res$padj,
                       p.value.ase=dds.h2.fs.res$padj,
                       genes=rownames(dds.fs.cre.clat.res),
                       trans_effect=trans_p.values.fdr.fs.h2, stringsAsFactors=F)

# classify allele specific expression

type.levels = c("conserved","ambiguous","trans-only","cis-only", "cis + trans (enhancing)", "cis x trans (compensatory)", "cis-trans (compensatory)")

fs_h2_ase$type = get_regulation_type(fs_h2_ase)

fs_h2_table = table(fs_h2_ase$type) ##summary table of # of genes in each category


#####Filtering and labelling chromosomal information

fs_h2_ase_chr = rownames_to_column(fs_h2_ase)
rownames(fs_h2_ase_chr) = fs_h2_ase_chr$rowname

fs_h2_ase_chr = na.omit(cbind(fs_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

fs_h2_ase_chr$A_X = "Autosomes"

fs_h2_ase_chr[(fs_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

##separating information for autosomes and x-chromosome
fs_h2_ase_autosome = fs_h2_ase_chr %>% filter(fs_h2_ase_chr$A_X == "Autosomes")
fs_h2_ase_X = fs_h2_ase_chr %>% filter(fs_h2_ase_chr$A_X == "X")

#write.csv(fs_h2_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fs_h2_ase_chr.csv")


fs_h2_ase_chr_count = as.data.frame(table(fs_h2_ase_chr$type, fs_h2_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fs_h2_ase_chr_count$Sample = "h2"
fs_h2_ase_chr_count$Tissue = "S"
fs_h2_ase_chr_count$Sex = "Female"
fs_h2_ase_chr_count$Name = paste(fs_h2_ase_chr_count$Sample, fs_h2_ase_chr_count$Sex, fs_h2_ase_chr_count$Tissue, fs_h2_ase_chr_count$Chromosome)

#autosomes vs X
fs_h2_a_x = fs_h2_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

###################################################3###################################################3
###MS###################################################3###################################################3
###################################################3###################################################3

#H1

keep = Hybrid == "H1" & tissue == "S" & sex == "M"
dds_h1.ms = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h1.ms = DESeq(dds_h1.ms)
dds.h1.ms = results(dds_h1.ms, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h1.ms.res = cbind(as.data.frame(dds.h1.ms),as.data.frame(dds.h1.ms) %>%
                        mutate(h1.ms = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h1.ms)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and H1 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.ms.h1 = cbindX(cpm(counts(dds_ms.cre.clat), log=T), cpm(counts(dds_h1.ms), log=T)) #these are objects before the result object
df.cpm.ms.h1 = na.omit(df.cpm.ms.h1)

trans_effect_ms.h1 = lapply(split(df.cpm.ms.h1, 1:dim(df.cpm.ms.h1)[1]), getTransEffects)
trans_p.values.fdr.ms.h1 = p.adjust(sapply(trans_effect_ms.h1, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.ms.h1)[1])
names(tmp) = rownames(df.cpm.ms.h1)
names(trans_p.values.fdr.ms.h1) = rownames(df.cpm.ms.h1)
tmp[ names(trans_p.values.fdr.ms.h1) ] = trans_p.values.fdr.ms.h1
trans_p.values.fdr.ms.h1 = tmp

genes=sort(rownames(dds.ms.cre.clat))
dds.ms.cre.clat.res = dds.ms.cre.clat.res[genes,]
dds.h1.ms.res = dds.h1.ms.res[genes, ]
rownames(dds.h1.ms.res) = genes
trans_p.values.fdr.ms.h1 = trans_p.values.fdr.ms.h1[genes]
names(trans_p.values.fdr.ms.h1)= genes

ms_h1_ase = data.frame(logFC.sp=dds.ms.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h1.ms.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h1.ms.res$h1.ms) ],
                       species=c("Clat","Cre")[ factor(dds.h1.ms.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.ms.cre.clat.res$padj,
                       p.value.ase=dds.h1.ms.res$padj,
                       genes=rownames(dds.ms.cre.clat.res),
                       trans_effect=trans_p.values.fdr.ms.h1, stringsAsFactors=F)

# classify allele specific expression

ms_h1_ase$type = get_regulation_type(ms_h1_ase) #contains both auto and X-linked (inccrrect X-linked classification)

#####Filtering and labelling chromosomal information

ms_h1_ase_chr = rownames_to_column(ms_h1_ase)
rownames(ms_h1_ase_chr) = ms_h1_ase_chr$rowname

ms_h1_ase_chr = na.omit(cbind(ms_h1_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

ms_h1_ase_chr$A_X = "Autosomes"

ms_h1_ase_chr[(ms_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X" #identifying X-linked genes to remove them 
ms_h1_ase_autosome = ms_h1_ase_chr %>% filter(ms_h1_ase_chr$A_X == "Autosomes")
nrow(ms_h1_ase_autosome) #10814 detectable on autosomes out of 11088 overall on autosomes

table(ms_h1_ase_autosome$type)
#write.csv(ms_h1_ase_autosome, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h1_ase_autosome.csv")

#######Change regulatory divergence for X-linked genes
##H1 have C.remanei X-chromosomes

###Pre classified variable from Classifying Inheritance.R file
###This variable has regulatory divergence for X-linked genes 

dds.ms.h1_X = na.omit(dds.ms.h1[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.ms.h1_X$reg.div = Cre_X_regulatory_divergence(dds.ms.h1_X) #has data only for X-linked genes
dds.ms.h1_X = dds.ms.h1_X %>% select("reg.div") %>% rownames_to_column() 
dds.ms.h1_X$chr_info = "X"
head(dds.ms.h1_X)
table(dds.ms.h1_X$reg.div)
nrow(dds.ms.h1_X) #2624 detectable out of possible 2665 X-linked genes

#write.csv(dds.ms.h1_X, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h1_ase_X_linked.csv")

#################################################################################
#H2
##MS
keep = Hybrid == "H2" & tissue == "S" & sex == "M"
dds_h2.ms = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h2.ms = DESeq(dds_h2.ms)
dds.h2.ms = results(dds_h2.ms, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h2.ms.res = cbind(as.data.frame(dds.h2.ms),as.data.frame(dds.h2.ms) %>%
                        mutate(h2.ms = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h2.ms)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and H2 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.ms.h2 = cbindX(cpm(counts(dds_ms.cre.clat), log=T), cpm(counts(dds_h2.ms), log=T)) #these are objects before the result object
df.cpm.ms.h2 = na.omit(df.cpm.ms.h2)

trans_effect_ms.h2 = lapply(split(df.cpm.ms.h2, 1:dim(df.cpm.ms.h2)[1]), getTransEffects)
trans_p.values.fdr.ms.h2 = p.adjust(sapply(trans_effect_ms.h2, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.ms.h2)[1])
names(tmp) = rownames(df.cpm.ms.h2)
names(trans_p.values.fdr.ms.h2) = rownames(df.cpm.ms.h2)
tmp[ names(trans_p.values.fdr.ms.h2) ] = trans_p.values.fdr.ms.h2
trans_p.values.fdr.ms.h2 = tmp

genes=sort(rownames(dds.ms.cre.clat))
dds.ms.cre.clat.res = dds.ms.cre.clat.res[genes,]
dds.h2.ms.res = dds.h2.ms.res[genes, ]
rownames(dds.h2.ms.res) = genes
trans_p.values.fdr.ms.h2 = trans_p.values.fdr.ms.h2[genes]
names(trans_p.values.fdr.ms.h2)= genes

ms_h2_ase = data.frame(logFC.sp=dds.ms.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h2.ms.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h2.ms.res$h2.ms) ],
                       species=c("Clat","Cre")[ factor(dds.h2.ms.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.ms.cre.clat.res$padj,
                       p.value.ase=dds.h2.ms.res$padj,
                       genes=rownames(dds.ms.cre.clat.res),
                       trans_effect=trans_p.values.fdr.ms.h2, stringsAsFactors=F)


# classify allele specific expression

ms_h2_ase$type = get_regulation_type(ms_h2_ase)

#######Change regulatory divergence for X-linked genes
##H2 have C.latens X-chromosomes
#####Filtering and labelling chromosomal information

ms_h2_ase_chr = rownames_to_column(ms_h2_ase)
rownames(ms_h2_ase_chr) = ms_h2_ase_chr$rowname

ms_h2_ase_chr = na.omit(cbind(ms_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

ms_h2_ase_chr$A_X = "Autosomes"

ms_h2_ase_chr[(ms_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X" #identifying X-linked genes to remove them 
ms_h2_ase_autosome = ms_h2_ase_chr %>% filter(ms_h2_ase_chr$A_X == "Autosomes")
nrow(ms_h2_ase_autosome) #10814 detectable on autosomes out of 11088 overall on autosomes

#write.csv(ms_h2_ase_autosome, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h2_ase_autosome.csv")

table(ms_h2_ase_autosome$type)
###Pre classified variable from Classifying Inheritance.R file
###This variable has regulatory divergence for X-linked genes 

dds.ms.h2_X = na.omit(dds.ms.h2[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.ms.h2_X$reg.div = Clat_X_regulatory_divergence(dds.ms.h2_X)
dds.ms.h2_X = dds.ms.h2_X %>% select("reg.div") %>% rownames_to_column() 
dds.ms.h2_X$chr_info = "X"
table(dds.ms.h2_X$reg.div)
head(dds.ms.h2_X)

nrow(dds.ms.h2_X) #2624 detectable out of possible 2665 X-linked genes

write.csv(dds.ms.h2_X, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h2_ase_X_linked.csv")



###################################################################################################
###WM

#H1
keep = Hybrid == "H1" & tissue == "W" & sex == "W"
dds_h1.wm = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h1.wm = DESeq(dds_h1.wm)
dds.h1.wm = results(dds_h1.wm, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h1.wm.res = cbind(as.data.frame(dds.h1.wm),as.data.frame(dds.h1.wm) %>%
                        mutate(h1.wm = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h1.wm)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and H1 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.wm.h1 = cbindX(cpm(counts(dds_wm.cre.clat), log=T), cpm(counts(dds_h1.wm), log=T)) #these are objects before the result object
df.cpm.wm.h1 = na.omit(df.cpm.wm.h1)

trans_effect_wm.h1 = lapply(split(df.cpm.wm.h1, 1:dim(df.cpm.wm.h1)[1]), getTransEffects)
trans_p.values.fdr.wm.h1 = p.adjust(sapply(trans_effect_wm.h1, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.wm.h1)[1])
names(tmp) = rownames(df.cpm.wm.h1)
names(trans_p.values.fdr.wm.h1) = rownames(df.cpm.wm.h1)
tmp[ names(trans_p.values.fdr.wm.h1) ] = trans_p.values.fdr.wm.h1
trans_p.values.fdr.wm.h1 = tmp

genes=sort(rownames(dds.wm.cre.clat))
dds.wm.cre.clat.res = dds.wm.cre.clat.res[genes,]
dds.h1.wm.res = dds.h1.wm.res[genes, ]
rownames(dds.h1.wm.res) = genes
trans_p.values.fdr.wm.h1 = trans_p.values.fdr.wm.h1[genes]
names(trans_p.values.fdr.wm.h1)= genes

wm_h1_ase = data.frame(logFC.sp=dds.wm.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h1.wm.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h1.wm.res$h1.wm) ],
                       species=c("Clat","Cre")[ factor(dds.h1.wm.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.wm.cre.clat.res$padj,
                       p.value.ase=dds.h1.wm.res$padj,
                       genes=rownames(dds.wm.cre.clat.res),
                       trans_effect=trans_p.values.fdr.wm.h1, stringsAsFactors=F)

# classify allele specific expression

wm_h1_ase$type = get_regulation_type(wm_h1_ase)
table(wm_h1_ase$type)


#####Filtering and labelling chromosomal information

wm_h1_ase_chr = rownames_to_column(wm_h1_ase)
rownames(wm_h1_ase_chr) = wm_h1_ase_chr$rowname

wm_h1_ase_chr = na.omit(cbind(wm_h1_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

wm_h1_ase_chr$A_X = "Autosomes"

wm_h1_ase_chr[(wm_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X" #identifying X-linked genes to remove them 
wm_h1_ase_autosome = wm_h1_ase_chr %>% filter(wm_h1_ase_chr$A_X == "Autosomes")
nrow(wm_h1_ase_autosome) #10824 detectable on autosomes out of 11088 overall on autosomes

table(wm_h1_ase_autosome$type)
write.csv(wm_h1_ase_autosome, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h1_ase_autosome.csv")

#######Change regulatory divergence for X-linked genes
##H1 have C.remanei X-chromosomes

###Pre classified variable from Classifying Inheritance.R file
###This variable has regulatory divergence for X-linked genes 

dds.wm.h1 = rownames_to_column(dds.wm.h1)
rownames(dds.wm.h1) = dds.wm.h1$rowname
dds.wm.h1_X = na.omit(dds.wm.h1[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])[, -1]
dds.wm.h1_X$reg.div = Cre_X_regulatory_divergence(dds.wm.h1_X)
dds.wm.h1_X = dds.wm.h1_X %>% select("reg.div") %>% rownames_to_column() 
dds.wm.h1_X$chr_info = "X"
table(dds.wm.h1_X$reg.div)
head(dds.wm.h1_X)

nrow(dds.wm.h1_X) #2626 detectable out of possible 2665 X-linked genes

write.csv(dds.wm.h1_X, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h1_ase_X_linked.csv")


###############################################################################################
#H2

keep = Hybrid == "H2" & tissue == "W" & sex == "W"
dds_h2.wm = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h2.wm = DESeq(dds_h2.wm)
dds.h2.wm = results(dds_h2.wm, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h2.wm.res = cbind(as.data.frame(dds.h2.wm),as.data.frame(dds.h2.wm) %>%
                        mutate(h2.wm = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h2.wm)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and h2 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.wm.h2 = cbindX(cpm(counts(dds_wm.cre.clat), log=T), cpm(counts(dds_h2.wm), log=T)) #these are objects before the result object
df.cpm.wm.h2 = na.omit(df.cpm.wm.h2)

trans_effect_wm.h2 = lapply(split(df.cpm.wm.h2, 1:dim(df.cpm.wm.h2)[1]), getTransEffects)
trans_p.values.fdr.wm.h2 = p.adjust(sapply(trans_effect_wm.h2, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.wm.h2)[1])
names(tmp) = rownames(df.cpm.wm.h2)
names(trans_p.values.fdr.wm.h2) = rownames(df.cpm.wm.h2)
tmp[ names(trans_p.values.fdr.wm.h2) ] = trans_p.values.fdr.wm.h2
trans_p.values.fdr.wm.h2 = tmp

genes=sort(rownames(dds.wm.cre.clat))
dds.wm.cre.clat.res = dds.wm.cre.clat.res[genes,]
dds.h2.wm.res = dds.h2.wm.res[genes, ]
rownames(dds.h2.wm.res) = genes
trans_p.values.fdr.wm.h2 = trans_p.values.fdr.wm.h2[genes]
names(trans_p.values.fdr.wm.h2)= genes

wm_h2_ase = data.frame(logFC.sp=dds.wm.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h2.wm.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h2.wm.res$h2.wm) ],
                       species=c("Clat","Cre")[ factor(dds.h2.wm.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.wm.cre.clat.res$padj,
                       p.value.ase=dds.h2.wm.res$padj,
                       genes=rownames(dds.wm.cre.clat.res),
                       trans_effect=trans_p.values.fdr.wm.h2, stringsAsFactors=F)

# classify allele specific expression

wm_h2_ase$type = get_regulation_type(wm_h2_ase)

#####Filtering and labelling chromosomal information

wm_h2_ase_chr = rownames_to_column(wm_h2_ase)
rownames(wm_h2_ase_chr) = wm_h2_ase_chr$rowname

wm_h2_ase_chr = na.omit(cbind(wm_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

wm_h2_ase_chr$A_X = "Autosomes"

wm_h2_ase_chr[(wm_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X" #identifying X-linked genes to remove them 
wm_h2_ase_autosome = wm_h2_ase_chr %>% filter(wm_h2_ase_chr$A_X == "Autosomes")
nrow(wm_h2_ase_autosome) #10824 detectable on autosomes out of 11088 overall on autosomes

table(wm_h2_ase_autosome$type)
write.csv(wm_h2_ase_autosome, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h2_ase_autosome.csv")


#######Change regulatory divergence for X-linked genes
##H2 have C.latens X-chromosomes

###Pre classified variable from Classifying Inheritance.R file
###This variable has regulatory divergence for X-linked genes 
dds.wm.h2 = rownames_to_column(dds.wm.h2)
row.names(dds.wm.h2) = dds.wm.h2$rowname
dds.wm.h2_X = na.omit(dds.wm.h2[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])[, -1]
dds.wm.h2_X$reg.div = Clat_X_regulatory_divergence(dds.wm.h2_X)
dds.wm.h2_X = dds.wm.h2_X %>% select("reg.div") %>% rownames_to_column() 
dds.wm.h2_X$chr_info = "X"

table(dds.wm.h2_X$reg.div)
head(dds.wm.h2_X)

nrow(dds.wm.h2_X) #2626 detectable out of possible 2665 X-linked genes

#write.csv(dds.wm.h2_X, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h2_ase_X_linked.csv")


###############################################
####MG
##H2

keep = Hybrid == "H2" & tissue == "G" & sex == "M"
dds_h2.mg = DESeqDataSetFromMatrix(countData = Hybrid_allele_counts[,keep], colData = coldata_h[keep,], design = ~ batch + species)
dds_h2.mg = DESeq(dds_h2.mg)
dds.h2.mg = results(dds_h2.mg, alpha = 0.05) ##setting FDR cutoff to 0.05
dds.h2.mg.res = cbind(as.data.frame(dds.h2.mg),as.data.frame(dds.h2.mg) %>%
                        mutate(h2.mg = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(h2.mg)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1


##Combining parental and h2 cpm data
##cbindX in gdata to combine dataframes with different number of rows
df.cpm.mg.h2 = cbindX(cpm(counts(dds_mg.cre.clat), log=T), cpm(counts(dds_h2.mg), log=T)) #these are objects before the result object
df.cpm.mg.h2 = na.omit(df.cpm.mg.h2)

trans_effect_mg.h2 = lapply(split(df.cpm.mg.h2, 1:dim(df.cpm.mg.h2)[1]), getTransEffects)
trans_p.values.fdr.mg.h2 = p.adjust(sapply(trans_effect_mg.h2, function(x) x[[4]] ), method="BH")

tmp = rep(NA, dim(df.cpm.mg.h2)[1])
names(tmp) = rownames(df.cpm.mg.h2)
names(trans_p.values.fdr.mg.h2) = rownames(df.cpm.mg.h2)
tmp[ names(trans_p.values.fdr.mg.h2) ] = trans_p.values.fdr.mg.h2
trans_p.values.fdr.mg.h2 = tmp

genes=sort(rownames(dds.mg.cre.clat))
dds.mg.cre.clat.res = dds.mg.cre.clat.res[genes,]
dds.h2.mg.res = dds.h2.mg.res[genes, ]
rownames(dds.h2.mg.res) = genes
trans_p.values.fdr.mg.h2 = trans_p.values.fdr.mg.h2[genes]
names(trans_p.values.fdr.mg.h2)= genes

mg_h2_ase = data.frame(logFC.sp=dds.mg.cre.clat.res$log2FoldChange,
                       logFC.ase=dds.h2.mg.res$log2FoldChange,
                       DE=c("sig","no sig","sig")[ factor(dds.h2.mg.res$h2.mg) ],
                       species=c("Clat","Cre")[ factor(dds.h2.mg.res$log2FoldChange > 0) ],
                       sex="female",
                       p.value.sp=dds.mg.cre.clat.res$padj,
                       p.value.ase=dds.h2.mg.res$padj,
                       genes=rownames(dds.mg.cre.clat.res),
                       trans_effect=trans_p.values.fdr.mg.h2, stringsAsFactors=F)


# classify allele specific expression

mg_h2_ase$type = get_regulation_type(mg_h2_ase)

#####Filtering and labelling chromosomal information

mg_h2_ase_chr = rownames_to_column(mg_h2_ase)
rownames(mg_h2_ase_chr) = mg_h2_ase_chr$rowname

mg_h2_ase_chr = na.omit(cbind(mg_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

mg_h2_ase_chr$A_X = "Autosomes"

mg_h2_ase_chr[(mg_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X" #identifying X-linked genes to remove them 
mg_h2_ase_autosome = mg_h2_ase_chr %>% filter(mg_h2_ase_chr$A_X == "Autosomes")
nrow(mg_h2_ase_autosome) #10480 detectable on autosomes out of 11088 overall on autosomes

table(mg_h2_ase_autosome$type)
write.csv(mg_h2_ase_autosome, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/mg_h2_ase_autosome.csv")


#######Change regulatory divergence for X-linked genes
##H2 have C.latens X-chromosomes

###Pre classified variable from Classifying Inheritance.R file
###This variable has regulatory divergence for X-linked genes 

head(dds.mg.h2)
dds.mg.h2_X = na.omit(dds.mg.h2[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.mg.h2_X$reg.div = Clat_X_regulatory_divergence(dds.mg.h2_X)
dds.mg.h2_X = dds.mg.h2_X %>% select("reg.div") %>% rownames_to_column() 
dds.mg.h2_X$chr_info = "X"
table(dds.mg.h2_X$reg.div)

head(dds.mg.h2_X)

nrow(dds.mg.h2_X) #2326 detectable out of possible 2665 X-linked genes

write.csv(dds.mg.h2_X, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/mg_h2_ase_X_linked.csv")


##############################################################################################################
###Making cumulative data 
fg_h1_table.prop = cbind(as.data.frame(fg_h1_table),as.data.frame(fg_h1_table) %>%
                           mutate(Proportion = prop.table(fg_h1_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("G", 7))) %>%
                           mutate(Hybrid = c(rep("H1", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))


fg_h2_table.prop = cbind(as.data.frame(fg_h2_table),as.data.frame(fg_h2_table) %>%
                           mutate(Proportion = prop.table(fg_h2_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("G", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))  

fs_h1_table.prop = cbind(as.data.frame(fs_h1_table),as.data.frame(fs_h1_table) %>%
                           mutate(Proportion = prop.table(fs_h1_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("S", 7))) %>%
                           mutate(Hybrid = c(rep("H1", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))

fs_h2_table.prop = cbind(as.data.frame(fs_h2_table),as.data.frame(fs_h2_table) %>%
                           mutate(Proportion = prop.table(fs_h2_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("S", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))

ms_h1_table.prop = cbind(as.data.frame(ms_h1_table),as.data.frame(ms_h1_table) %>%
                           mutate(Proportion = prop.table(ms_h1_table)) %>%
                           mutate(Sex = c(rep("M", 8))) %>%
                           mutate(Tissue = c(rep("S", 8))) %>%
                           mutate(Hybrid = c(rep("H1", 8))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

ms_h2_table.prop = cbind(as.data.frame(ms_h2_table),as.data.frame(ms_h2_table) %>%
                           mutate(Proportion = prop.table(ms_h2_table)) %>%
                           mutate(Sex = c(rep("M", 8))) %>%
                           mutate(Tissue = c(rep("S", 8))) %>%
                           mutate(Hybrid = c(rep("H2", 8))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 
wm_h1_table.prop = cbind(as.data.frame(wm_h1_table),as.data.frame(wm_h1_table) %>%
                           mutate(Proportion = prop.table(wm_h1_table)) %>%
                           mutate(Sex = c(rep("M", 8))) %>%
                           mutate(Tissue = c(rep("W", 8))) %>%
                           mutate(Hybrid = c(rep("H1", 8))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

wm_h2_table.prop = cbind(as.data.frame(wm_h2_table),as.data.frame(wm_h2_table) %>%
                           mutate(Proportion = prop.table(wm_h2_table)) %>%
                           mutate(Sex = c(rep("M", 8))) %>%
                           mutate(Tissue = c(rep("W", 8))) %>%
                           mutate(Hybrid = c(rep("H2", 8))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))

mg_h2_table.prop = cbind(as.data.frame(mg_h2_table),as.data.frame(mg_h2_table) %>%
                           mutate(Proportion = prop.table(mg_h2_table)) %>%
                           mutate(Sex = c(rep("M", 8))) %>%
                           mutate(Tissue = c(rep("G", 8))) %>%
                           mutate(Hybrid = c(rep("H2", 8))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 
#View(mg_h2_table)
regdiv_all_counts = rbind(fg_h1_table.prop, fg_h2_table.prop, fs_h1_table.prop, fs_h2_table.prop,
                          ms_h1_table.prop, ms_h2_table.prop, mg_h2_table.prop, wm_h1_table.prop, wm_h2_table.prop)

#write.csv(regdiv_all_counts, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/RegDiv_all_counts.csv")



