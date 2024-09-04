###Differential gene expression analysis of ortholgous genes
library(DESeq2)
library(ggplot2)
library(cowplot) #add on to ggplot for better themes and figure customization
library(lemon) #to work with legends and axes in ggplot2
library(dplyr)
library(gdata)
library(RColorBrewer)
library(colorBlindness)
library(colorspace)
library(tidyverse)

theme_set(theme_classic())
setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT
orthologs = read.table("new_1to1_orthologgenelist.txt", sep="\t", head=T, comment.char="#")
head(orthologs)
cre_ortho = orthologs[, 1]
clat_ortho = orthologs[ ,2]
cumulative_ortho = orthologs[,3]

head(cre_ortho)
nrow(orthologs)
head(clat_ortho)

############
##DATA PREP#
############

#C.remanei
cre_counts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_counts_ortho = cre_counts[cre_ortho, -c(13:15)] ##not including Whole animal samples
rownames(cre_counts_ortho) = cumulative_ortho #changing row names
nrow(cre_counts_ortho)
colnames(cre_counts_ortho)

#C.latens
clat_counts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_counts_ortho = clat_counts[clat_ortho, -c(13:15)]#removing WM samples
rownames(clat_counts_ortho) = cumulative_ortho
nrow(clat_counts_ortho)
names(clat_counts_ortho)
nrow(clat_counts)
#Combining all data
wt_counts_orthologs = cbind(cre_counts_ortho, clat_counts_ortho)
colnames(wt_counts_orthologs)
colnames(wt_counts_orthologs) = gsub("X", "", colnames(wt_counts_orthologs)) ##searching and replacing the column names 

#Coldata prep
sample_name = colnames(wt_counts_orthologs)
tissue = substr(sample_name, 3,3) 
sex = substr(sample_name, 2,2)
species =  c(rep("Cre", 12), rep("Clat", 12))
batch = c(rep(c(1, 2 , 3), 8))
coldata = data.frame(sample_name, species, sex, tissue, batch)
coldata$species = factor(coldata$species)
coldata$tissue = factor(coldata$tissue)
coldata$sex = factor(coldata$sex)
coldata$batch = factor(coldata$batch)
coldata ##WM not included

#Run deseq2
dds_wt = DESeqDataSetFromMatrix(countData = wt_counts_orthologs, colData = coldata, design = ~ species +
                                  tissue + sex + species:tissue + sex:tissue + species:sex )
dds_wt = DESeq(dds_wt)
resultsNames(dds_wt)

# get the model matrix
mod_mat_wt <- model.matrix(design(dds_wt), colData(dds_wt))

############################################################################################
#1.Effect of Species on gene expression 
############################################################################################
#Define coefficient vectors for each condition
cre = colMeans(mod_mat_wt[dds_wt$species == "Cre", ])
clat = colMeans(mod_mat_wt[dds_wt$species == "Clat", ])


crexclat = results(dds_wt, contrast = cre - clat, alpha = 0.05) #baseline clat
crexclat_res <- cbind(as.data.frame(crexclat), as.data.frame(crexclat) %>%
                        mutate(crexclat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(crexclat))
crexclat_res = na.omit(crexclat_res)
crexclat_res1 = crexclat_res %>% filter(crexclat_res$crexclat != 0)
crexclat_res2 = crexclat_res %>% filter(crexclat_res$crexclat == 0)

summary(crexclat)##upregulated in Cre, reference level is Clatens
nrow(crexclat_res)
table(crexclat_res$crexclat) ###Number of species-biased and conserved genes


#############################################################################################################################
#########CONVERTING GENE NAMES FOR GPROFILER#################################################################################
#############################################################################################################################
####Saving species-biased genes names to view in GPROFILER 

###Loading cremanie and clatens gene names for older genomes (from Daniel)
# remanei_gene_names = read.table("DESeq2/remanei_gene_names.txt", fill = T)
# latens_gene_names = read.table("DESeq2/latens_gene_names.txt", fill = T)
# head(latens_gene_names)
# 
# ###getting species-biased and unbiased genes from DESeq2 (above) analysis
# crexclat_up = crexclat_res %>% filter(crexclat_res$crexclat == 1) #Higher in C. remanei
# crexclat_down = crexclat_res %>% filter(crexclat_res$crexclat == -1) #Higher in C. latens
# crexclat_nobias = crexclat_res %>% filter(crexclat_res$crexclat == 0) #Unbiased
# crexclat_speciesbiased = crexclat_res %>% filter(crexclat_res$crexclat != 0) #Unbiased
# crexclat_speciesbiased =rownames_to_column(crexclat_speciesbiased)
# crexclat_speciesbiased = separate(crexclat_speciesbiased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
# crexclat_speciesbiased = separate(crexclat_speciesbiased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
# head(crexclat_speciesbiased)
# 
# ##Saving species-baised and unbiased genes
# # write.table(crexclat_up, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/crexclat_up.txt",
# #             row.names = TRUE, col.names = TRUE, quote = FALSE)
# # 
# # write.table(crexclat_down, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/crexclat_down.txt",
# #             row.names = TRUE, col.names = TRUE, quote = FALSE)
# # 
# # write.table(crexclat_speciesbiased, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/species_biased.txt",
# #             row.names = TRUE, col.names = TRUE, quote = FALSE)
# 
# ###Getting gene names for genes with higher expressison in C. remanei
# cre_biased = rownames_to_column(crexclat_up)
# cre_biased = separate(cre_biased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
# cre_biased = separate(cre_biased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2] #run it right after rprevious step
# cre_biased_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% cre_biased$Crenames)[,2] #Cre old gene names
# cre_biased_clat = subset(latens_gene_names, latens_gene_names$V4 %in% cre_biased$Clatnames)[,2] #corresponding gene name sin C. latens
# 
# 
# # write.table(cre_biased_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/cre_biased_cregenes.txt", 
# #             row.names = FALSE, col.names = FALSE, quote = FALSE)
# # 
# # write.table(cre_biased_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/cre_biased_clatgenes.txt", 
# #             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# ###Getting gene names for genes with higher expresison in C. latens
# clat_biased = rownames_to_column(crexclat_down)
# clat_biased = separate(clat_biased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
# clat_biased = separate(clat_biased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2] #run it right after rprevious step
# clat_biased_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% clat_biased$Crenames)[,2] #Cre old gene names
# clat_biased_clat = subset(latens_gene_names, latens_gene_names$V4 %in% clat_biased$Clatnames)[,2] #corresponding gene name sin C. latens
# 
# nrow(clat_biased)
# length(clat_biased_clat)
# # write.table(clat_biased_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/clat_biased_cregenes.txt",
# #             row.names = FALSE, col.names = FALSE, quote = FALSE)
# # 
# # write.table(clat_biased_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/clat_biased_clatgenes.txt",
# #             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# 
# ###Getting gene names for genes with no diffrence in expression across species
# sp_unbiased = rownames_to_column(crexclat_nobias)
# sp_unbiased = separate(sp_unbiased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
# sp_unbiased = separate(sp_unbiased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2] #run it right after rprevious step
# sp_unbiased_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% sp_unbiased$Crenames)[,2] #Cre old gene names
# sp_unbiased_clat = subset(latens_gene_names, latens_gene_names$V4 %in% sp_unbiased$Clatnames)[,2] #corresponding gene name sin C. latens
# length(sp_unbiased_cre)
# nrow(sp_unbiased)
# length(sp_unbiased_clat)
# head(sp_unbiased)
# write.table(sp_unbiased_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/sp_unbiased_cregenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(sp_unbiased_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/sp_unbiased_clatgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(sp_unbiased, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/sp_unbiased.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

####################################################################################################################################
####PART I##########################################################################################################################
####################################################################################################################################


####FIGURE 1A-VISUALIZING SPECIES-BIASED GENES #####################################################################################
####################################################################################################################################

table(crexclat_res$crexclat)

sp_genes_count = data.frame(gene_category = c("Conserved", "Conserved", "DEG","DEG"),
                            genes_count = c(6865, 0, 3339, 3224), #creating a proxy for conserved deg cetagory
                            DEG_category = c("Clat-b", "Cre-b","Clat-b","Cre-b"))

ggplot(sp_genes_count, aes(x = sp_genes_count$gene_category, y = sp_genes_count$genes_count, fill = sp_genes_count$DEG_category)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,7500,500))+
  ggtitle("Species-biased genes")


####################################################################################################################################
####FIGURE 1B-VISUALIZING SPECIES-BIASED GENES ON CHROMOSOMES#######################################################################
####################################################################################################################################

##adding chromosomal information to the dataframe
chrom_info = read.table("Cremanie_genenames_and_chromosome.txt", fill = T)
head(chrom_info)

t##adding chromosomal information to the dataframe
chrom_info = read.table("Cremanie_genenames_and_chromosome.txt", fill = T)
head(chrom_info)
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

##changing rownames to new column
head(crexclat_res)
crexclat_res_a = na.omit(crexclat_res[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ])
crexclat_res_a = rownames_to_column(crexclat_res_a)
nrow(crexclat_res_a)
head(orthologs_chr)
rownames(orthologs_chr) = orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name  
crexclat_res_chr = (na.omit(orthologs_chr[crexclat_res_a$rowname, ]) %>% cbind(crexclat_res_a))
#View(crexclat_res_chr)

spgenes_freq = as.data.frame(table(crexclat_res_chr$chromosome, crexclat_res_chr$crexclat))
#View(spgenes_freq)


####################################################################################################
#####Plotting species-biased genes across chromosomes###############################################
####################################################################################################
ggplot(spgenes_freq, aes(x = spgenes_freq$Var1, y = spgenes_freq$Freq, fill = spgenes_freq$Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,6500,500))+
  ggtitle("Species-biased genes")

#ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 1/species_biased_chr.pdf", dpi = 300, width = 18, height = 8 )

####################################################################################################################################
#####Proportion of species-biased genes on chromosomes###################################################3

##plotting proportions across chromosomes
library(ggstats)

ggplot(spgenes_freq) +
  aes(x = spgenes_freq$Var1, fill = spgenes_freq$Var2, weight = spgenes_freq$Freq, by = spgenes_freq$Var1) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))

#ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 1/fig1b_chr_degs.pdf", dpi = 300, width = 18, height = 8 )

################################################################################
#FINDING GENES THAT ARE TISSUE-BIASED (INCLUDING THOSE THAT ARE DUE TO INTERATION TERM) 
################################################################################

#Define coefficient vectors for each condition
gonad = colMeans(mod_mat_wt[dds_wt$tissue == "G", ])
soma = colMeans(mod_mat_wt[dds_wt$tissue == "S", ])

gxs = results(dds_wt, contrast = soma - gonad, alpha = 0.05) #baseline clat
gxs_res <- cbind(as.data.frame(gxs), as.data.frame(gxs) %>%
                   mutate(gxs = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                   dplyr::select(gxs))
nrow(gxs_res)
summary(gxs)##upregulated in Soma, reference level is gonad
gxs_res = na.omit(gxs_res)
gxs_res1 = gxs_res %>% filter(gxs_res$gxs != 0)
gxs_res2 = gxs_res %>% filter(gxs_res$gxs == 0)

table(gxs_res$gxs)
design(dds_wt)
#write.csv(gxs_res, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 1/gxs_res_genes.csv")

#####################################################################################################################
##VISUALIZING TISSUE-BIASED GENES################################################################################
ts_genes_count = data.frame(gene_category = c("Conserved", "Conserved", "DEG","DEG"),
           genes_count = c(3469, 0, 4494, 5465), #creating a proxy for conserved deg cetagory
           DEG_category = c("Soma-biased", "Gonad-biased","Soma-biased","Gonad-biased"))

ggplot(ts_genes_count, aes(x = ts_genes_count$gene_category, y = ts_genes_count$genes_count, fill = ts_genes_count$DEG_category)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,10500,500))+
  ggtitle("Tissue-biased genes")

##############################################################################################
#####PLOTTING TISSUE-BIASED GENES ON CHROMOSOMES############################################
###############################################################################################

##adding chromosomal information to the dataframe
chrom_info = read.table("Cremanie_genenames_and_chromosome.txt", fill = T)
head(chrom_info)

table(chrom_info$V1)

##changing rownames to new column

#gxs_res3 = rownames_to_column(gxs_res)

gxs_res_a = na.omit(gxs_res[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ])
gxs_res_a = rownames_to_column(gxs_res_a)
nrow(gxs_res_a)
head(orthologs_chr)
rownames(orthologs_chr) = orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name  
gxs_res_chr = (na.omit(orthologs_chr[gxs_res_a$rowname, ]) %>% cbind(gxs_res_a))
View(gxs_res_chr)

tsgenes_freq = as.data.frame(table(gxs_res_chr$chromosome, gxs_res_chr$gxs))
View(tsgenes_freq)


####################################################################################################
#####Plotting tissue-biased genes across chromosomes###############################################
####################################################################################################
ggplot(tsgenes_freq, aes(x = tsgenes_freq$Var1, y = tsgenes_freq$Freq, fill = tsgenes_freq$Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,6500,500))+
  ggtitle("Tissue-biased genes")

####################################################################################################################################
#####Proportion of tissue-biased genes on chromosomes###################################################3

##plotting proportions across chromosomes
library(ggstats)

ggplot(tsgenes_freq) +
  aes(x = tsgenes_freq$Var1, fill = tsgenes_freq$Var2, weight = tsgenes_freq$Freq, by = tsgenes_freq$Var1) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Tissue-biased genes across chromosomes")

###################################################################################################################



################################################################################
#FINDING GENES THAT ARE SEX-BIASED (INCLUDING THOSE THAT ARE DUE TO INTERATION TERM) 
################################################################################

#Define coefficient vectors for each condition
male = colMeans(mod_mat_wt[dds_wt$sex == "M", ])
female = colMeans(mod_mat_wt[dds_wt$sex == "F", ])

mxf = results(dds_wt, contrast = male - female, alpha = 0.05) #baseline clat
mxf_res <- cbind(as.data.frame(mxf), as.data.frame(mxf) %>%
                   mutate(mxf = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                   dplyr::select(mxf))
nrow(mxf_res)##upregulated in Cre, reference level is Clatens
summary(mxf)
mxf_res = na.omit(mxf_res)
mxf_res1 = mxf_res %>% filter(mxf_res$mxf != 0)
mxf_res2 = mxf_res %>% filter(mxf_res$mxf == 0)

table(mxf_res$mxf)
#write.csv(mxf_res, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 1/mxf_res_genes.csv")

#####################################################################################################################
##VISUALIZING SEX-BIASED GENES################################################################################
mf_genes_count = data.frame(gene_category = c("Conserved", "Conserved", "DEG","DEG"),
                            genes_count = c(4111, 0, 5184, 4133), #creating a proxy for conserved deg cetagory
                            DEG_category = c("Male-biased", "Female-biased","Male-biased","Female-biased"))

ggplot(mf_genes_count, aes(x = mf_genes_count$gene_category, y = mf_genes_count$genes_count, fill = mf_genes_count$DEG_category)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,10500,500))+
  ggtitle("Sex-biased genes")


##############################################################################################
#####PLOTTING SEX-BIASED GENES ON CHROMOSOMES############################################
###############################################################################################

##adding chromosomal information to the dataframe
chrom_info = read.table("Cremanie_genenames_and_chromosome.txt", fill = T)
head(chrom_info)

table(chrom_info$V1)

##changing rownames to new column

mxf_res3 = rownames_to_column(mxf_res)

mxf_res_a = na.omit(mxf_res[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ])
mxf_res_a = rownames_to_column(mxf_res_a)
nrow(mxf_res_a)
head(orthologs_chr)
rownames(orthologs_chr) = orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name  
mxf_res_chr = (na.omit(orthologs_chr[mxf_res_a$rowname, ]) %>% cbind(mxf_res_a))
View(mxf_res_chr)

sexgenes_freq = as.data.frame(table(mxf_res_chr$chromosome, mxf_res_chr$mxf))
View(sexgenes_freq)

####################################################################################################
#####Plotting sex-biased genes across chromosomes###############################################
####################################################################################################
ggplot(sexgenes_freq, aes(x = sexgenes_freq$Var1, y = sexgenes_freq$Freq, fill = sexgenes_freq$Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,6500,500))+
  ggtitle("Sex-biased genes")

####################################################################################################################################
#####Proportion of tissue-biased genes on chromosomes###################################################3

##plotting proportions across chromosomes
library(ggstats)

ggplot(sexgenes_freq) +
  aes(x = sexgenes_freq$Var1, fill = sexgenes_freq$Var2, weight = sexgenes_freq$Freq, by = sexgenes_freq$Var1) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Sex-biased genes across chromosomes")

###################################################################################################################



