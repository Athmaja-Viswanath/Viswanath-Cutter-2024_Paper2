library(DESeq2)
library(ggplot2)
library(cowplot) #add on to ggplot for better themes and figure customization
library(lemon) #to work with legends and axes in ggplot2
library(dplyr)
theme_set(theme_classic())
library(gdata)
library(tidyverse)
setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

orthologs = read.table("new_1to1_orthologgenelist.txt", sep="\t", head=T, comment.char="#")
head(orthologs)
cre_ortho = orthologs[, 1]
clat_ortho = orthologs[ ,2]
cumulative_ortho = orthologs[,3]

#C.remanei
cre_counts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_counts_ortho = cre_counts[cre_ortho, ]
rownames(cre_counts_ortho) = cumulative_ortho
nrow(cre_counts_ortho)

# head(cre_ortho)
# head(cre_counts)
# head(cre_counts_ortho)
# cre_counts[c("remanei_genome_00002719"), ]
# cre_counts_ortho[c("remanei_genome_00002719"), ]

#C.latens
clat_counts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_counts_ortho = clat_counts[clat_ortho, ]
rownames(clat_counts_ortho) = cumulative_ortho
nrow(clat_counts_ortho)

##H1 
H1_ase= read.table("H1_ase_counts.txt", sep="\t", head=T, row.name=1, comment.char = "#")
H1_genecounts = cbind(H1_ase[,seq(1,24,2)] + H1_ase[,seq(2,24,2)])#adding allele specific readcounts to get one counts/gen
H1_genecounts = H1_genecounts[cumulative_ortho, ]
head(H1_genecounts)
names(H1_ase)

###H2
H2_ase= read.table("H2_ase_counts.txt", sep="\t", head=T, row.name=1, comment.char = "#")
H2_genecounts = cbind(H2_ase[,seq(1,30,2)] + H2_ase[,seq(2,30,2)])#adding allele specific readcounts to get one counts/gen
H2_genecounts = H2_genecounts[cumulative_ortho, ] #adds togerth but retain first columns name
head(H2_genecounts)


#Combining all data

counts_orthologs = cbind(cre_counts_ortho, clat_counts_ortho, H1_genecounts, H2_genecounts)
colnames(counts_orthologs)

dim(counts_orthologs)


##Prep
tissue_type.g = c(rep("G", 6)) 
tissue_type.s = c(rep("S", 6))
tissue_type.gs = c(rep("G", 3), rep("S", 3))

Sex.fm = c(rep("F", 3), rep("M", 3))
Sex.m = c(rep("M", 6))
Sex.f = c(rep("F", 6))

batch_wt = c("1", "2", "3", "1", "2", "3")

species = c(rep(c("Cre","Clat"), each=3))
species_cre_h1 = c(rep(c("Cre","H1"), each=3))
species_clat_h1 = c(rep(c("Clat","H1"), each=3))
species_cre_h2 = c(rep(c("Cre","H2"), each=3))
species_clat_h2 = c(rep(c("Clat","H2"), each=3))
species_h1_h2 = c(rep(c("H1","H2"), each=3))

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
View(dds.fg.cre.clat.res)
View(dds.fg.cre.clat.res %>% filter(dds.fg.cre.clat.res$fg.cre.clat == 0))
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

table(dds.fs.cre.clat.res$fs.cre.clat)


##MG X MG
colnames(counts_orthologs)
mg.cre.clat = counts_orthologs[ , c(7,8,9,22,23,24)]
colnames(mg.cre.clat)

samples_mg.cre.clat = colnames(mg.cre.clat)
nrow(mg.cre.clat)

coldata_mg.cre.clat = data.frame(samples_mg.cre.clat, Sex.m, tissue_type.g, batch_wt, species_wt)

dds_mg.cre.clat = DESeqDataSetFromMatrix(countData = mg.cre.clat, colData = coldata_mg.cre.clat, design = ~batch_wt + species_wt)
dds_mg.cre.clat = DESeq(dds_mg.cre.clat)
dds.mg.cre.clat= results(dds_mg.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_mg.cre.clat, alpha = 0.05))
dds.mg.cre.clat.res = cbind(as.data.frame(dds.mg.cre.clat),as.data.frame(dds.mg.cre.clat) %>%
                              mutate(mg.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                              dplyr::select(mg.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.mg.cre.clat.res = dds.mg.cre.clat.res[order(row.names(dds.mg.cre.clat.res)), ]

nrow(dds.mg.cre.clat.res)

table(dds.mg.cre.clat.res$mg.cre.clat)

###MS X MS
colnames(counts_orthologs)
ms.cre.clat = counts_orthologs[ , c(10,11,12,25,26,27)]
colnames(ms.cre.clat)

samples_ms.cre.clat = colnames(ms.cre.clat)
nrow(ms.cre.clat)

coldata_ms.cre.clat = data.frame(samples_ms.cre.clat, Sex.m, tissue_type.s, batch_wt, species_wt)

dds_ms.cre.clat = DESeqDataSetFromMatrix(countData = ms.cre.clat, colData = coldata_ms.cre.clat, design = ~batch_wt + species_wt)
dds_ms.cre.clat = DESeq(dds_ms.cre.clat)
dds.ms.cre.clat= results(dds_ms.cre.clat, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_ms.cre.clat, alpha = 0.05))
dds.ms.cre.clat.res = cbind(as.data.frame(dds.ms.cre.clat),as.data.frame(dds.ms.cre.clat) %>%
                              mutate(ms.cre.clat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                              dplyr::select(ms.cre.clat)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.ms.cre.clat.res = dds.ms.cre.clat.res[order(row.names(dds.ms.cre.clat.res)), ]
#View(dds.ms.cre.clat.res)
table(dds.ms.cre.clat.res$ms.cre.clat)

##WM X WM

colnames(counts_orthologs)
wm.cre.clat = counts_orthologs[ , c(13,14,15,28,29,30)]
colnames(wm.cre.clat)

samples_wm.cre.clat = colnames(wm.cre.clat)
nrow(wm.cre.clat)

coldata_wm.cre.clat = data.frame(samples_wm.cre.clat, Sex.m, tissue_type.s, batch_wt, species_wt)

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

##C.remanei x H1

#FG x FG
colnames(counts_orthologs)
fg.cre.H1 = counts_orthologs[ , c(1,2,3,31,32,33)]
colnames(fg.cre.H1)

samples_fg.cre.H1 = colnames(fg.cre.H1)
nrow(fg.cre.H1)

coldata_fg.cre.H1 = data.frame(samples_fg.cre.H1, Sex.f, tissue_type.g, batch_wt, species_cre_h1)

dds_fg.cre.H1 = DESeqDataSetFromMatrix(countData = fg.cre.H1, colData = coldata_fg.cre.H1, design = ~batch_wt + species_cre_h1)
dds_fg.cre.H1 = DESeq(dds_fg.cre.H1)
dds.fg.cre.H1= results(dds_fg.cre.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fg.cre.H1, alpha = 0.05))
dds.fg.cre.H1.res = cbind(as.data.frame(dds.fg.cre.H1),as.data.frame(dds.fg.cre.H1) %>%
                        mutate(fg.cre.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(fg.cre.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.fg.cre.H1.res = dds.fg.cre.H1.res[order(row.names(dds.fg.cre.H1.res)), ]
#View(dds.fg.cre.H1.res)

##FS X FS

colnames(counts_orthologs)
fs.cre.H1 = counts_orthologs[ , c(4,5,6,34,35,36)]
colnames(fs.cre.H1)

samples_fs.cre.H1 = colnames(fs.cre.H1)
nrow(fs.cre.H1)

coldata_fs.cre.H1 = data.frame(samples_fs.cre.H1, Sex.f, tissue_type.s, batch_wt, species_cre_h1)

dds_fs.cre.H1 = DESeqDataSetFromMatrix(countData = fs.cre.H1, colData = coldata_fs.cre.H1, design = ~batch_wt + species_cre_h1)
dds_fs.cre.H1 = DESeq(dds_fs.cre.H1)
dds.fs.cre.H1= results(dds_fs.cre.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fs.cre.H1, alpha = 0.05))
dds.fs.cre.H1.res = cbind(as.data.frame(dds.fs.cre.H1),as.data.frame(dds.fs.cre.H1) %>%
                            mutate(fs.cre.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                            dplyr::select(fs.cre.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.fs.cre.H1.res = dds.fs.cre.H1.res[order(row.names(dds.fs.cre.H1.res)), ]
#View(dds.fs.cre.H1.res)
table(dds.fs.cre.H1.res$fs.cre.H1)

###MS X MS
colnames(counts_orthologs)
ms.cre.H1 = counts_orthologs[ , c(10,11,12,40,41,42)]
colnames(ms.cre.H1)

samples_ms.cre.H1 = colnames(ms.cre.H1)
nrow(ms.cre.H1)

coldata_ms.cre.H1 = data.frame(samples_ms.cre.H1, Sex.m, tissue_type.s, batch_wt, species_cre_h1)

dds_ms.cre.H1 = DESeqDataSetFromMatrix(countData = ms.cre.H1, colData = coldata_ms.cre.H1, design = ~batch_wt + species_cre_h1)
dds_ms.cre.H1 = DESeq(dds_ms.cre.H1)
dds.ms.cre.H1= results(dds_ms.cre.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_ms.cre.H1, alpha = 0.05))
dds.ms.cre.H1.res = cbind(as.data.frame(dds.ms.cre.H1),as.data.frame(dds.ms.cre.H1) %>%
                        mutate(ms.cre.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(ms.cre.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.ms.cre.H1.res = dds.ms.cre.H1.res[order(row.names(dds.ms.cre.H1.res)), ]
#View(dds.ms.cre.H1.res)
table(dds.ms.cre.H1.res$ms.cre.H1)

##WM X WM

colnames(counts_orthologs)
wm.cre.H1 = counts_orthologs[ , c(13,14,15,37,38,39)]
colnames(wm.cre.H1)

samples_wm.cre.H1 = colnames(wm.cre.H1)
nrow(wm.cre.H1)

coldata_wm.cre.H1 = data.frame(samples_wm.cre.H1, Sex.m, tissue_type.s, batch_wt, species_cre_h1)

dds_wm.cre.H1 = DESeqDataSetFromMatrix(countData = wm.cre.H1, colData = coldata_wm.cre.H1, design = ~batch_wt + species_cre_h1)
dds_wm.cre.H1 = DESeq(dds_wm.cre.H1)
dds.wm.cre.H1= results(dds_wm.cre.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_wm.cre.H1, alpha = 0.05))
dds.wm.cre.H1.res = cbind(as.data.frame(dds.wm.cre.H1),as.data.frame(dds.wm.cre.H1) %>%
                        mutate(wm.cre.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(wm.cre.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

table(dds.wm.cre.H1.res$wm.cre.H1)
dds.wm.cre.H1.res = dds.wm.cre.H1.res[order(row.names(dds.wm.cre.H1.res)), ]
#View(dds.wm.cre.H1.res)

##C. latens x H1

#FG x FG
colnames(counts_orthologs)
fg.clat.H1 = counts_orthologs[ , c(16,17,18,31,32,33)]
colnames(fg.clat.H1)

samples_fg.clat.H1 = colnames(fg.clat.H1)
nrow(fg.clat.H1)

coldata_fg.clat.H1 = data.frame(samples_fg.clat.H1, Sex.f, tissue_type.g, batch_wt, species_clat_h1)

dds_fg.clat.H1 = DESeqDataSetFromMatrix(countData = fg.clat.H1, colData = coldata_fg.clat.H1, design = ~batch_wt + species_clat_h1)
dds_fg.clat.H1 = DESeq(dds_fg.clat.H1)
dds.fg.clat.H1= results(dds_fg.clat.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fg.clat.H1, alpha = 0.05))
dds.fg.clat.H1.res = cbind(as.data.frame(dds.fg.clat.H1),as.data.frame(dds.fg.clat.H1) %>%
                         mutate(fg.clat.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(fg.clat.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.fg.clat.H1.res = dds.fg.clat.H1.res[order(row.names(dds.fg.clat.H1.res)), ]

##Clatens x H1
#FS X FS

colnames(counts_orthologs)
fs.clat.H1 = counts_orthologs[ , c(19,20,21,34,35,36)]
colnames(fs.clat.H1)

samples_fs.clat.H1 = colnames(fs.clat.H1)
nrow(fs.clat.H1)
coldata_fs.clat.H1 = data.frame(samples_fs.clat.H1, Sex.f, tissue_type.s, batch_wt, species_clat_h1)

dds_fs.clat.H1 = DESeqDataSetFromMatrix(countData = fs.clat.H1, colData = coldata_fs.clat.H1, design = ~batch_wt + species_clat_h1)
dds_fs.clat.H1 = DESeq(dds_fs.clat.H1)
dds.fs.clat.H1= results(dds_fs.clat.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fs.clat.H1, alpha = 0.05))
dds.fs.clat.H1.res = cbind(as.data.frame(dds.fs.clat.H1),as.data.frame(dds.fs.clat.H1) %>%
                         mutate(fs.clat.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(fs.clat.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.fs.clat.H1.res = dds.fs.clat.H1.res[order(row.names(dds.fs.clat.H1.res)), ]

###MS X MS
colnames(counts_orthologs)
ms.clat.H1 = counts_orthologs[ , c(25,26,27,40,41,42)]
colnames(ms.clat.H1)

samples_ms.clat.H1 = colnames(ms.clat.H1)
nrow(ms.clat.H1)

coldata_ms.clat.H1 = data.frame(samples_ms.clat.H1, Sex.f, tissue_type.g, batch_wt, species_clat_h1)

dds_ms.clat.H1 = DESeqDataSetFromMatrix(countData = ms.clat.H1, colData = coldata_ms.clat.H1, design = ~batch_wt + species_clat_h1)
dds_ms.clat.H1 = DESeq(dds_ms.clat.H1)
dds.ms.clat.H1= results(dds_ms.clat.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_ms.clat.H1, alpha = 0.05))
dds.ms.clat.H1.res = cbind(as.data.frame(dds.ms.clat.H1),as.data.frame(dds.ms.clat.H1) %>%
                         mutate(ms.clat.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(ms.clat.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.ms.clat.H1.res = dds.ms.clat.H1.res[order(row.names(dds.ms.clat.H1.res)), ]

##WM X WM

colnames(counts_orthologs)
wm.clat.H1 = counts_orthologs[ , c(28,29,30,37,38,39)]
colnames(wm.clat.H1)

samples_wm.clat.H1 = colnames(wm.clat.H1)
nrow(wm.clat.H1)

coldata_wm.clat.H1 = data.frame(samples_wm.clat.H1, Sex.m, tissue_type.s, batch_wt, species_clat_h1)

dds_wm.clat.H1 = DESeqDataSetFromMatrix(countData = wm.clat.H1, colData = coldata_wm.clat.H1, design = ~batch_wt + species_clat_h1)
dds_wm.clat.H1 = DESeq(dds_wm.clat.H1)
dds.wm.clat.H1= results(dds_wm.clat.H1, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_wm.clat.H1, alpha = 0.05))
dds.wm.clat.H1.res = cbind(as.data.frame(dds.wm.clat.H1),as.data.frame(dds.wm.clat.H1) %>%
                         mutate(wm.clat.H1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(wm.clat.H1)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.wm.clat.H1.res = dds.wm.clat.H1.res[order(row.names(dds.wm.clat.H1.res)), ]
#View(dds.wm.clat.H1.res)

table(dds.wm.clat.H1.res$wm.clat.H1)
#####Dividing genes into inheritance categories 

classify_inheritance <- function(x){
  y = paste(x[,1], x[,2], x[,3])
  cl = rep(NA, dim(x)[1])
  cl[ y == "-1 0 -1" | y == "1 0 1" ] = "C. remanei dominant" ## 3 number are : Cre x Clat (in terms of Cre, clat baseline), Cre x h (in terms of h), Clat x h (interms of H)
  cl[ y == "-1 1 0" | y == "1 -1 0" ] = "C. latens dominant"
  cl[ y == "0 1 1" | y == "-1 1 1" | y == "1 1 1"] = "Overdominant"
  cl[ y == "0 -1 -1" | y == "-1 -1 -1" | y == "1 -1 -1" ] = "Underdominant"
  cl[ y == "-1 1 -1" | y == "1 -1 1"] = "Additive"
  cl[ y == "0 0 0"] = "no change"
  cl[ is.na(cl) ] = "ambiguous"
  return(cl)
}

###H1
#FG 
##all these columns/DEG analyses are in the RScript "Differential gene expression analysis_new annotations"
dds.fg.h1 = data.frame(dds.fg.cre.clat.res["fg.cre.clat"], dds.fg.cre.H1.res["fg.cre.H1"], dds.fg.clat.H1.res["fg.clat.H1"]) #taking the newly created columns together
dds.fg.h1 = na.omit(dds.fg.h1)
dim(dds.fg.h1)
dds.fg.h1$inheritance.fg = classify_inheritance(dds.fg.h1)
head(dds.fg.h1)
nrow(dds.fg.h1)
fg_h1_inht = table(dds.fg.h1$inheritance.fg)


##FS
dds.fs.h1 = data.frame(dds.fs.cre.clat.res["fs.cre.clat"], dds.fs.cre.H1.res["fs.cre.H1"], dds.fs.clat.H1.res["fs.clat.H1"]) #taking the newly created columns together
dds.fs.h1 = na.omit(dds.fs.h1)
dim(dds.fs.h1)
dds.fs.h1$inheritance.fs = classify_inheritance(dds.fs.h1)
head(dds.fs.h1)
fs_h1_inht = table(dds.fs.h1$inheritance.fs)

##MS

dds.ms.h1 = cbindX(dds.ms.cre.clat.res["ms.cre.clat"], dds.ms.cre.H1.res["ms.cre.H1"], dds.ms.clat.H1.res["ms.clat.H1"]) #taking the newly created columns together
dds.ms.h1 = na.omit(dds.ms.h1)
dim(dds.ms.h1)
dds.ms.h1$inheritance.ms = classify_inheritance(dds.ms.h1)
head(dds.ms.h1)
ms_h1_inht = table(dds.ms.h1$inheritance.ms)

#WM
dds.wm.h1 = cbind(dds.wm.cre.clat.res["wm.cre.clat"], dds.wm.cre.H1.res["wm.cre.H1"], dds.wm.clat.H1.res["wm.clat.H1"]) #taking the newly created columns together
dds.wm.h1 = na.omit(dds.wm.h1)
dim(dds.wm.h1)
dds.wm.h1$inheritance.wm = classify_inheritance(dds.wm.h1)
head(dds.wm.h1)
wm_h1_inht = table(dds.wm.h1$inheritance.wm)



##########################################################################################################################
##H2

##2. Cre x H2
#FG x FG
colnames(counts_orthologs)
fg.cre.H2 = counts_orthologs[ , c(1,2,3,43,44,45)]
colnames(fg.cre.H2)

samples_fg.cre.H2 = colnames(fg.cre.H2)
nrow(fg.cre.H2)

coldata_fg.cre.H2 = data.frame(samples_fg.cre.H2, Sex.f, tissue_type.g, batch_wt, species_cre_h2)

dds_fg.cre.H2 = DESeqDataSetFromMatrix(countData = fg.cre.H2, colData = coldata_fg.cre.H2, design = ~batch_wt + species_cre_h2)
dds_fg.cre.H2 = DESeq(dds_fg.cre.H2)
dds.fg.cre.H2= results(dds_fg.cre.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fg.cre.H2, alpha = 0.05))
dds.fg.cre.H2.res = cbind(as.data.frame(dds.fg.cre.H2),as.data.frame(dds.fg.cre.H2) %>%
                        mutate(fg.cre.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(fg.cre.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.fg.cre.H2.res = dds.fg.cre.H2.res[order(row.names(dds.fg.cre.H2.res)), ]

##FS X FS

colnames(counts_orthologs)
fs.cre.H2 = counts_orthologs[ , c(4,5,6,46,47,48)]
colnames(fs.cre.H2)

samples_fs.cre.H2 = colnames(fs.cre.H2)
nrow(fs.cre.H2)

coldata_fs.cre.H2 = data.frame(samples_fs.cre.H2, Sex.f, tissue_type.s, batch_wt, species_cre_h2)

dds_fs.cre.H2 = DESeqDataSetFromMatrix(countData = fs.cre.H2, colData = coldata_fs.cre.H2, design = ~batch_wt + species_cre_h2)
dds_fs.cre.H2 = DESeq(dds_fs.cre.H2)
dds.fs.cre.H2= results(dds_fs.cre.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fs.cre.H2, alpha = 0.05))
dds.fs.cre.H2.res = cbind(as.data.frame(dds.fs.cre.H2),as.data.frame(dds.fs.cre.H2) %>%
                        mutate(fs.cre.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(fs.cre.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.fs.cre.H2.res = dds.fs.cre.H2.res[order(row.names(dds.fs.cre.H2.res)), ]

###MS X MS
colnames(counts_orthologs)
ms.cre.H2 = counts_orthologs[ , c(10,11,12,52,53,54)]
colnames(ms.cre.H2)

samples_ms.cre.H2 = colnames(ms.cre.H2)
nrow(ms.cre.H2)

coldata_ms.cre.H2 = data.frame(samples_ms.cre.H2, Sex.m, tissue_type.s, batch_wt, species_cre_h2)

dds_ms.cre.H2 = DESeqDataSetFromMatrix(countData = ms.cre.H2, colData = coldata_ms.cre.H2, design = ~batch_wt + species_cre_h2)
dds_ms.cre.H2 = DESeq(dds_ms.cre.H2)
dds.ms.cre.H2= results(dds_ms.cre.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_ms.cre.H2, alpha = 0.05))
dds.ms.cre.H2.res = cbind(as.data.frame(dds.ms.cre.H2),as.data.frame(dds.ms.cre.H2) %>%
                        mutate(ms.cre.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(ms.cre.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.ms.cre.H2.res = dds.ms.cre.H2.res[order(row.names(dds.ms.cre.H2.res)), ]

##WM X WM

colnames(counts_orthologs)
wm.cre.H2 = counts_orthologs[ , c(13,14,15,55,56,57)]
colnames(wm.cre.H2)

samples_wm.cre.H2 = colnames(wm.cre.H2)
nrow(wm.cre.H2)

coldata_wm.cre.H2 = data.frame(samples_wm.cre.H2, Sex.m, tissue_type.s, batch_wt, species_cre_h2)

dds_wm.cre.H2 = DESeqDataSetFromMatrix(countData = wm.cre.H2, colData = coldata_wm.cre.H2, design = ~batch_wt + species_cre_h2)
dds_wm.cre.H2 = DESeq(dds_wm.cre.H2)
dds.wm.cre.H2= results(dds_wm.cre.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_wm.cre.H2, alpha = 0.05))
dds.wm.cre.H2.res = cbind(as.data.frame(dds.wm.cre.H2),as.data.frame(dds.wm.cre.H2) %>%
                        mutate(wm.cre.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(wm.cre.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.wm.cre.H2.res = dds.wm.cre.H2.res[order(row.names(dds.wm.cre.H2.res)), ]
table(dds.wm.cre.H2.res$wm.cre.H2)

##MG X MG
colnames(counts_orthologs)
mg.cre.H2 = counts_orthologs[ , c(7,8,9,49,50,51)]
colnames(mg.cre.H2)

samples_mg.cre.H2 = colnames(mg.cre.H2)
nrow(mg.cre.H2)

coldata_mg.cre.H2 = data.frame(samples_mg.cre.H2, Sex.m, tissue_type.g, batch_wt, species_cre_h2)

dds_mg.cre.H2 = DESeqDataSetFromMatrix(countData = mg.cre.H2, colData = coldata_mg.cre.H2, design = ~batch_wt + species_cre_h2)
dds_mg.cre.H2 = DESeq(dds_mg.cre.H2)
dds.mg.cre.H2= results(dds_mg.cre.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_mg.cre.H2, alpha = 0.05))
dds.mg.cre.H2.res = cbind(as.data.frame(dds.mg.cre.H2),as.data.frame(dds.mg.cre.H2) %>%
                        mutate(mg.cre.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(mg.cre.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.mg.cre.H2.res = dds.mg.cre.H2.res[order(row.names(dds.mg.cre.H2.res)), ]

################################################################################
##C. latens x H2

#FG x FG
colnames(counts_orthologs)
fg.clat.H2 = counts_orthologs[ , c(16,17,18,43,44,45)]
colnames(fg.clat.H2)

samples_fg.clat.H2 = colnames(fg.clat.H2)
nrow(fg.clat.H2)

coldata_fg.clat.H2 = data.frame(samples_fg.clat.H2, Sex.f, tissue_type.g, batch_wt, species_clat_h2)

dds_fg.clat.H2 = DESeqDataSetFromMatrix(countData = fg.clat.H2, colData = coldata_fg.clat.H2, design = ~batch_wt + species_clat_h2)
dds_fg.clat.H2 = DESeq(dds_fg.clat.H2)
dds.fg.clat.H2= results(dds_fg.clat.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fg.clat.H2, alpha = 0.05))
dds.fg.clat.H2.res = cbind(as.data.frame(dds.fg.clat.H2),as.data.frame(dds.fg.clat.H2) %>%
                         mutate(fg.clat.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(fg.clat.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.fg.clat.H2.res = dds.fg.clat.H2.res[order(row.names(dds.fg.clat.H2.res)), ]

#FS X FS
colnames(counts_orthologs)
fs.clat.H2 = counts_orthologs[ , c(19,20,21,46,47,48)]
colnames(fs.clat.H2)

samples_fs.clat.H2 = colnames(fs.clat.H2)
nrow(fs.clat.H2)

coldata_fs.clat.H2 = data.frame(samples_fs.clat.H2, Sex.f, tissue_type.s, batch_wt, species_clat_h2)

dds_fs.clat.H2 = DESeqDataSetFromMatrix(countData = fs.clat.H2, colData = coldata_fs.clat.H2, design = ~batch_wt + species_clat_h2)
dds_fs.clat.H2 = DESeq(dds_fs.clat.H2)
dds.fs.clat.H2= results(dds_fs.clat.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_fs.clat.H2, alpha = 0.05))
dds.fs.clat.H2.res = cbind(as.data.frame(dds.fs.clat.H2),as.data.frame(dds.fs.clat.H2) %>%
                         mutate(fs.clat.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(fs.clat.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.fs.clat.H2.res = dds.fs.clat.H2.res[order(row.names(dds.fs.clat.H2.res)), ]

###MS X MS
colnames(counts_orthologs)
ms.clat.H2 = counts_orthologs[ , c(25,26,27,52,53,54)]
colnames(ms.clat.H2)

samples_ms.clat.H2 = colnames(ms.clat.H2)
nrow(ms.clat.H2)

coldata_ms.clat.H2 = data.frame(samples_ms.clat.H2, Sex.m, tissue_type.s, batch_wt, species_clat_h2)

dds_ms.clat.H2 = DESeqDataSetFromMatrix(countData = ms.clat.H2, colData = coldata_ms.clat.H2, design = ~batch_wt + species_clat_h2)
dds_ms.clat.H2 = DESeq(dds_ms.clat.H2)
dds.ms.clat.H2= results(dds_ms.clat.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_ms.clat.H2, alpha = 0.05))
dds.ms.clat.H2.res = cbind(as.data.frame(dds.ms.clat.H2),as.data.frame(dds.ms.clat.H2) %>%
                         mutate(ms.clat.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(ms.clat.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.ms.clat.H2.res = dds.ms.clat.H2.res[order(row.names(dds.ms.clat.H2.res)), ]

##WM X WM

colnames(counts_orthologs)
wm.clat.H2 = counts_orthologs[ , c(28,29,30,55,56,57)]
colnames(wm.clat.H2)

samples_wm.clat.H2 = colnames(wm.clat.H2)
nrow(wm.clat.H2)

coldata_wm.clat.H2 = data.frame(samples_wm.clat.H2, Sex.m, tissue_type.s, batch_wt, species_clat_h2)

dds_wm.clat.H2 = DESeqDataSetFromMatrix(countData = wm.clat.H2, colData = coldata_wm.clat.H2, design = ~batch_wt + species_clat_h2)
dds_wm.clat.H2 = DESeq(dds_wm.clat.H2)
dds.wm.clat.H2= results(dds_wm.clat.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_wm.clat.H2, alpha = 0.05))
dds.wm.clat.H2.res = cbind(as.data.frame(dds.wm.clat.H2),as.data.frame(dds.wm.clat.H2) %>%
                         mutate(wm.clat.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(wm.clat.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1

dds.wm.clat.H2.res = dds.wm.clat.H2.res[order(row.names(dds.wm.clat.H2.res)), ]
table(dds.wm.clat.H2.res$wm.clat.H2)

##MG X MG
colnames(counts_orthologs)
mg.clat.H2 = counts_orthologs[ , c(22,23,24,49,50,51)]
colnames(mg.clat.H2)

samples_mg.clat.H2 = colnames(mg.clat.H2)
nrow(mg.clat.H2)

coldata_mg.clat.H2 = data.frame(samples_mg.clat.H2, Sex.m, tissue_type.g, batch_wt, species_clat_h2)

dds_mg.clat.H2 = DESeqDataSetFromMatrix(countData = mg.clat.H2, colData = coldata_mg.clat.H2, design = ~batch_wt + species_clat_h2)
dds_mg.clat.H2 = DESeq(dds_mg.clat.H2)
dds.mg.clat.H2= results(dds_mg.clat.H2, alpha = 0.05) ##setting FDR cutoff to 0.05
summary(results(dds_mg.clat.H2, alpha = 0.05))
dds.mg.clat.H2.res = cbind(as.data.frame(dds.mg.clat.H2),as.data.frame(dds.mg.clat.H2) %>%
                         mutate(mg.clat.H2 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                         dplyr::select(mg.clat.H2)) #added another column where if padj > 0.05, value is 0, if not look at log2foldchange, if log2fold > 0 then 1 else -1
dds.mg.clat.H2.res = dds.mg.clat.H2.res[order(row.names(dds.mg.clat.H2.res)), ]

#####Categorization of inheritance for H2############################################################
#FG 

dds.fg.h2 = data.frame(dds.fg.cre.clat.res["fg.cre.clat"], dds.fg.cre.H2.res["fg.cre.H2"], dds.fg.clat.H2.res["fg.clat.H2"]) #taking the newly created columns together
dds.fg.h2 = na.omit(dds.fg.h2)
dim(dds.fg.h2)
dds.fg.h2$inheritance.fg = classify_inheritance(dds.fg.h2)
head(dds.fg.h2)
fg_h2_inht = table(dds.fg.h2$inheritance.fg)

##FS
dds.fs.h2 = data.frame(dds.fs.cre.clat.res["fs.cre.clat"], dds.fs.cre.H2.res["fs.cre.H2"], dds.fs.clat.H2.res["fs.clat.H2"]) #taking the newly created columns together
dds.fs.h2 = na.omit(dds.fs.h2)
dim(dds.fs.h2)
dds.fs.h2$inheritance.fs = classify_inheritance(dds.fs.h2)
head(dds.fs.h2)
fs_h2_inht = table(dds.fs.h2$inheritance.fs)

#MG 
dds.mg.h2 = cbindX(dds.mg.cre.clat.res["mg.cre.clat"], dds.mg.cre.H2.res["mg.cre.H2"], dds.mg.clat.H2.res["mg.clat.H2"]) #taking the newly created columns together
dds.mg.h2 = na.omit(dds.mg.h2)
dim(dds.mg.h2)
dds.mg.h2$inheritance.mg = classify_inheritance(dds.mg.h2)
head(dds.mg.h2)
mg_h2_inht = table(dds.mg.h2$inheritance.mg)

##MS
dds.ms.h2 = cbindX(dds.ms.cre.clat.res["ms.cre.clat"], dds.ms.cre.H2.res["ms.cre.H2"], dds.ms.clat.H2.res["ms.clat.H2"]) #taking the newly created columns together
dds.ms.h2 = na.omit(dds.ms.h2)
dim(dds.ms.h2)
dds.ms.h2$inheritance.ms = classify_inheritance(dds.ms.h2)
head(dds.ms.h2)
ms_h2_inht = table(dds.ms.h2$inheritance.ms)

#WM
dds.wm.h2 = cbind(dds.wm.cre.clat.res["wm.cre.clat"], dds.wm.cre.H2.res["wm.cre.H2"], dds.wm.clat.H2.res["wm.clat.H2"]) #taking the newly created columns together
dds.wm.h2 = na.omit(dds.wm.h2)
dim(dds.wm.h2)
dds.wm.h2$inheritance.wm = classify_inheritance(dds.wm.h2)
head(dds.wm.h2)
wm_h2_inht = table(dds.wm.h2$inheritance.wm)

###Making cumulative data 
fg_h1_inht.prop = cbind(as.data.frame(fg_h1_inht),as.data.frame(fg_h1_inht) %>%
                          mutate(Proportion = prop.table(fg_h1_inht)) %>%
                          mutate(Sex = c(rep("F", 7))) %>%
                          mutate(Tissue = c(rep("G", 7))) %>%
                          mutate(Hybrid = c(rep("H1", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid))


fg_h2_inht.prop = cbind(as.data.frame(fg_h2_inht),as.data.frame(fg_h2_inht) %>%
                          mutate(Proportion = prop.table(fg_h2_inht)) %>%
                          mutate(Sex = c(rep("F", 7))) %>%
                          mutate(Tissue = c(rep("G", 7))) %>%
                          mutate(Hybrid = c(rep("H2", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid))  

fs_h1_inht.prop = cbind(as.data.frame(fs_h1_inht),as.data.frame(fs_h1_inht) %>%
                          mutate(Proportion = prop.table(fs_h1_inht)) %>%
                          mutate(Sex = c(rep("F", 7))) %>%
                          mutate(Tissue = c(rep("S", 7))) %>%
                          mutate(Hybrid = c(rep("H1", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid))

fs_h2_inht.prop = cbind(as.data.frame(fs_h2_inht),as.data.frame(fs_h2_inht) %>%
                          mutate(Proportion = prop.table(fs_h2_inht)) %>%
                          mutate(Sex = c(rep("F", 7))) %>%
                          mutate(Tissue = c(rep("S", 7))) %>%
                          mutate(Hybrid = c(rep("H2", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid))

ms_h1_inht.prop = cbind(as.data.frame(ms_h1_inht),as.data.frame(ms_h1_inht) %>%
                          mutate(Proportion = prop.table(ms_h1_inht)) %>%
                          mutate(Sex = c(rep("M", 7))) %>%
                          mutate(Tissue = c(rep("S", 7))) %>%
                          mutate(Hybrid = c(rep("H1", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

ms_h2_inht.prop = cbind(as.data.frame(ms_h2_inht),as.data.frame(ms_h2_inht) %>%
                          mutate(Proportion = prop.table(ms_h2_inht)) %>%
                          mutate(Sex = c(rep("M", 7))) %>%
                          mutate(Tissue = c(rep("S", 7))) %>%
                          mutate(Hybrid = c(rep("H2", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid)) 
wm_h1_inht.prop = cbind(as.data.frame(wm_h1_inht),as.data.frame(wm_h1_inht) %>%
                          mutate(Proportion = prop.table(wm_h1_inht)) %>%
                          mutate(Sex = c(rep("M", 7))) %>%
                          mutate(Tissue = c(rep("W", 7))) %>%
                          mutate(Hybrid = c(rep("H1", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

wm_h2_inht.prop = cbind(as.data.frame(wm_h2_inht),as.data.frame(wm_h2_inht) %>%
                          mutate(Proportion = prop.table(wm_h2_inht)) %>%
                          mutate(Sex = c(rep("M", 7))) %>%
                          mutate(Tissue = c(rep("W", 7))) %>%
                          mutate(Hybrid = c(rep("H2", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid))

mg_h2_inht.prop = cbind(as.data.frame(mg_h2_inht),as.data.frame(mg_h2_inht) %>%
                          mutate(Proportion = prop.table(mg_h2_inht)) %>%
                          mutate(Sex = c(rep("M", 7))) %>%
                          mutate(Tissue = c(rep("G", 7))) %>%
                          mutate(Hybrid = c(rep("H2", 7))) %>%
                          dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

inheritance_all_counts = rbind(fg_h1_inht.prop, fg_h2_inht.prop, fs_h1_inht.prop, fs_h2_inht.prop,
                               ms_h1_inht.prop, ms_h2_inht.prop, wm_h1_inht.prop, wm_h2_inht.prop, mg_h2_inht.prop)

# write.csv(inheritance_all_counts, file = "Allele Specific Expression/Hybrid Analysis/Basic Results/inheritance_all_counts.csv")
# 
# write.csv(inheritance_all_counts, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/inheritance_all_counts.csv")


##################################################################################################################
#####CLASSIFYING REGULATORY DIVERGENCE FOR X-LINKED GENES OF F1 HYBRID MALES############################
#########################################################################################################

#####Classifying Regulatory divergecne for X-linked genes in Males

##For males from H1 where X-chromosomes if from C. remanei mother
Cre_X_regulatory_divergence <- function(x){
  y = paste(x[,1], x[,2], x[,3])
  cl = rep(NA, dim(x)[1])
  cl[ y == "-1 0 -1" | y == "1 0 1" ] = "cis-only" ## 3 number are : Cre x Clat (in terms of Cre, clat baseline), Cre x h (in terms of h), Clat x h (interms of H)
  cl[ y == "-1 1 0" | y == "1 -1 0" ] = "trans-only"
  cl[ y == "0 1 1" | y == "0 -1 -1" | y == "0 -1 1" | y == "0 1 -1"] = "cis-trans (compensatory)"
  cl[ y == "0 0 0"] = "conserved"
  cl[ is.na(cl) ] = "Other"
  return(cl)
}

##For males from H2 where X-chromosomes if from C. latens mother
Clat_X_regulatory_divergence <- function(x){
  y = paste(x[,1], x[,2], x[,3])
  cl = rep(NA, dim(x)[1])
  cl[ y == "-1 0 -1" | y == "1 0 1" ] = "trans-only" ## 3 number are : Cre x Clat (in terms of Cre, clat baseline), Cre x h (in terms of h), Clat x h (interms of H)
  cl[ y == "-1 1 0" | y == "1 -1 0" ] = "cis-only"
  cl[ y == "0 1 1" | y == "0 -1 -1" | y == "0 -1 1" | y == "0 1 -1"] = "cis-trans (compensatory)"
  cl[ y == "0 0 0"] = "conserved"
  cl[ is.na(cl) ] = "Other"
  return(cl)
}

##########################################################################################
##adding chromosomal information to the orthologous gene information
chrom_info = read.table("Cremanie_genenames_and_chromosome.txt", fill = T) ##no X-chromosome in C.latens so using C.r emanei X-shormosome for both
chrom_info_X = subset(chrom_info, V1 =="X")
chrom_info_I = subset(chrom_info, V1 =="I")
chrom_info_II = subset(chrom_info, V1 =="II")
chrom_info_III = subset(chrom_info, V1 =="III")
chrom_info_IV = subset(chrom_info, V1 =="IV")
chrom_info_V = subset(chrom_info, V1 =="V")

###get X-linked orthologous genes
ortho_X = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_X$V2)
ortho_X$chromosome = "X"
ortho_I = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_I$V2)
ortho_I$chromosome = "I"
ortho_II = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_II$V2)
ortho_II$chromosome = "II"
ortho_III = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_III$V2)
ortho_III$chromosome = "III"
ortho_IV = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_IV$V2)
ortho_IV$chromosome = "IV"
ortho_V = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_V$V2)
ortho_V$chromosome = "V"
nrow(ortho_X)

orthologs_chr = rbind(ortho_I, ortho_II, ortho_III, ortho_IV, ortho_V, ortho_X)
nrow(orthologs_chr) ##3loss of genes on scaffolds
nrow(orthologs)
nrow(ortho_X)

# write.table(orthologs_chr, file = "orthologs_chr.txt", sep = "\t",quote = FALSE,
#            row.names = FALSE)


##########################################################################################################
###########H2 MALE SAMPLES WHICH HAVE C.LATENS X-CHROMOSOME###############################################
##########################################################################################################

##########################################################################################################
#####H2 MS
##########################################################################################################

dim(dds.ms.h2)

dds.ms.h2_X = na.omit(dds.ms.h2[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.ms.h2_X$reg.div = Clat_X_regulatory_divergence(dds.ms.h2_X)
table(dds.ms.h2_X$reg.div)
names(dds.ms.h2_X)

##########################################################################################################
#####H2 MG
##########################################################################################################
##changing rownames to new column
dim(dds.mg.h2)

dds.mg.h2_X = na.omit(dds.mg.h2[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.mg.h2_X$reg.div = Clat_X_regulatory_divergence(dds.mg.h2_X)
table(dds.mg.h2_X$reg.div)


##########################################################################################################
#####H2 WM
##########################################################################################################
##changing rownames to new column
head(dds.wm.h2)
dds.wm.h2_X = na.omit(dds.wm.h2[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.wm.h2_X$reg.div = Clat_X_regulatory_divergence(dds.wm.h2_X)
table(dds.wm.h2_X$reg.div)


##########################################################################################################
###########H1 MALE SAMPLES WHICH HAVE C.REMANEI X-CHROMOSOME###############################################
##########################################################################################################

##########################################################################################################
#####H1 MS
##########################################################################################################

##changing rownames to new column
head(dds.ms.h1)

dds.ms.h1_X = na.omit(dds.ms.h1[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.ms.h1_X$reg.div = Cre_X_regulatory_divergence (dds.ms.h1_X)
table(dds.ms.h1_X$reg.div)

##########################################################################################################
#####H1 WM
##########################################################################################################
##changing rownames to new column
head(dds.wm.h1)
dds.wm.h1_X = na.omit(dds.wm.h1[ortho_X$C..remanei.Gene.name_C.latens.Gene.name, ])
dds.wm.h1_X$reg.div = Cre_X_regulatory_divergence(dds.wm.h1_X)
table(dds.wm.h1_X$reg.div)

####Move to "Allele Specific expression" RScript

