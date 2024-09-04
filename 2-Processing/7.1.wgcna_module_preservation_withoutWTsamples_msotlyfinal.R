library(WGCNA)
library(DESeq2)
#library(GEOquery) for data from NCBI
library(tidyverse)
library(ggplot2)
library(CorLevelPlot)
library(gridExtra)
options(stringsAsFactors = FALSE)

setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT
orthologs = read.table("new_1to1_orthologgenelist.txt", sep="\t", head=T, comment.char="#")
head(orthologs)
cre_ortho = orthologs[, 1]
clat_ortho = orthologs[ ,2]
cumulative_ortho = orthologs[,3]

head(cre_ortho)
nrow(orthologs)

###H1
##DATA PREP
H1_ase= read.table("H1_ase_counts.txt", sep="\t", head=T, row.name=1, comment.char = "#")

#H1_ase = H1_ase[,seq(1,24,2)] + H1_ase[,seq(2,24,2)]
H1_ase = H1_ase[cumulative_ortho, ]

colnames(H1_ase)
nrow(H1_ase)

H2_ase = read.table("H2_ase_counts.txt", sep="\t", head=T, row.name=1, comment.char = "#")
#H2_ase = H2_ase[,seq(1,30,2)] + H2_ase[,seq(2,30,2)]
H2_ase = H2_ase[cumulative_ortho, ]

ncol(H2_ase)
colnames(H2_ase)
nrow(H2_ase)

#C.remanei
cre_counts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_counts_ortho = cre_counts[cre_ortho, ]
rownames(cre_counts_ortho) = cumulative_ortho #changing row names
nrow(cre_counts_ortho)

#C.latens
clat_counts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_counts_ortho = clat_counts[clat_ortho, ]
rownames(clat_counts_ortho) = cumulative_ortho
nrow(clat_counts)

###Combining data
counts_orthologs = cbind(cre_counts_ortho, clat_counts_ortho, H1_ase[,seq(1,24,2)] + H1_ase[,seq(2,24,2)], H2_ase[,seq(1,30,2)] + H2_ase[,seq(2,30,2)])

colnames(counts_orthologs)

colnames(counts_orthologs) = gsub("_Cre", "", colnames(counts_orthologs)) ##searching and replacing the column names 
colnames(counts_orthologs) = gsub("X", "", colnames(counts_orthologs)) ##searching and replacing the column names 
colnames(counts_orthologs)

#Coldata prep

sample_name = colnames(counts_orthologs)
tissue = substr(sample_name, 3,3) 
tissue = gsub("M", "W", tissue)
tissue[37] = "W"
tissue[38] = "W"
tissue[39] = "W"
sex = substr(sample_name, 2,2)
sex[37] = "W"
sex[38] = "W"
sex[39] = "W"
species =  c(rep("Cre", 15), rep("Clat", 15),rep("H1", 12), rep("H2", 15))
batch = c(rep(c(1, 2 , 3), 19))
coldata = data.frame(sample_name, species, sex, tissue, batch)
#View(coldata)


head(counts_orthologs)


###QUALITY CONTROL TO DETECT OUTLIERS ACROSS ALL SAMPLES

gsg = goodSamplesGenes(t(counts_orthologs))
summary(gsg)

gsg$allOK ##gives outliers, if TRUE = all genes and samples are good and have no outliers

table(gsg$goodGenes) #number of genes that are outliers = FALSE
table(gsg$goodSamples)


counts_orthologs = counts_orthologs[gsg$goodGenes == TRUE, ] ##filtering good genes
nrow(counts_orthologs)


##Normalization across all the samples
###Create a deseq2 dataset

#create coldata if needed

coldata2 = coldata[, -1]     ##changing column to rownames
rownames(coldata2) = coldata[ , 1]

all(rownames(coldata2) %in% colnames(counts_orthologs))
all(rownames(coldata2) == colnames(counts_orthologs))

dds = DESeqDataSetFromMatrix(countData = counts_orthologs, 
                             colData = coldata2, 
                             design = ~1) ##not specifiying a model
#View(coldata2)
####removing low count genes

dds75 = dds[rowSums(counts(dds) >= 15) >= 42, ]

nrow(dds75)#11935 genes 


##variance stabilize transformaiton 

dds_norm = vst(dds75)


#get normalized counts
#transfromed data for downstream analysis

norma.counts = assay(dds_norm) %>% #transposed and normalized counts for all samples
  t()

#View(norma.counts)

##########################################################################################
####MODULE PRESERVATION ANALYSIS########################################################
################################################################################################

##H2
###Getting normalized count for each data separately
norma.counts_cre=norma.counts[1:15, ]
norma.counts_cre=norma.counts_cre[-(13:15), ]
norma.counts_clat=norma.counts[16:30, ]
norma.counts_clat=norma.counts_clat[-(13:15), ]
norma.counts_h2=norma.counts[43:57, ]
norma.counts_h2=norma.counts_h2[-(13:15), ]


row.names(norma.counts_h2)
# We now set up the multi-set expression data

setLabels = c("Cre", "H2", "Clat")
multiExpr=list(Cre=list(data=norma.counts_cre),
               H2=list(data=norma.counts_h2),
               Clat = list(data=norma.counts_clat))


###QUALITY CONTROL TO DETECT OUTLIERS FROM THE MULITSET DATA

###checking samples for mutliExpr data

gsg_t = goodSamplesGenesMS(multiExpr, verbose = 3)

gsg_t$allOK

table(gsg_t$goodGenes) #number of genes that are outliers = FALSE
table(gsg_t$goodSamples)

##Filtering base genes from both datasets
norma.counts_h2 = norma.counts_h2[,gsg_t$goodGenes == TRUE] ##filtering good genes
norma.counts_cre = norma.counts_cre[,gsg_t$goodGenes == TRUE]
norma.counts_clat = norma.counts_clat[,gsg_t$goodGenes == TRUE]


##remaking the data with good genes
multiExpr=list(Cre=list(data=norma.counts_cre),
               H2=list(data=norma.counts_h2),
               Clat = list(data=norma.counts_clat))


#####Network construction

##Choose a set of soft threshold powers

power = c(c(1:10), seq(from = 12, to = 50, by = 2))

sft_cre_h2 = pickSoftThreshold(norma.counts_cre,
                               powerVector = power,
                               networkType = "signed",
                               verbose = 5)

sft_clat_h2 = pickSoftThreshold(norma.counts_clat,
                                powerVector = power,
                                networkType = "signed",
                                verbose = 5)

sft.data_cre_h2 = sft_cre_h2$fitIndices ##we will use the R suqared value and mean connectivity , max R^2 and min mean connectivity
sft.data_clat_h2 = sft_clat_h2$fitIndices


###Visualization to pick the right power
a1_cre_h2 = ggplot(sft.data_cre_h2, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2_cre_h2 = ggplot(sft.data_cre_h2, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1_cre_h2, a2_cre_h2, nrow = 2)  ###24 with WT and 22 without WT

#picking soft threshold for clat network

a1_clat_h2 = ggplot(sft.data_clat_h2, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2_clat_h2 = ggplot(sft.data_clat_h2, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1_clat_h2, a2_clat_h2, nrow = 2) #22 with WT, 20 without WT

###converting matrix to numeric

norma.counts_cre[] = sapply(norma.counts_cre, as.numeric)  ##different step than official tutorial
norma.counts_clat[] = sapply(norma.counts_clat, as.numeric)

soft_power_cre_h2 = 22
soft_power_clat_h2 = 20

temp_cor = cor #to prevent WGCNa from using other cor function

cor = WGCNA::cor


bwnet_cre_h2 = blockwiseModules(norma.counts_cre, 
                                maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                                TOMType = "signed",
                                networkType = "signed",
                                power = soft_power_cre_h2,
                                mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                                numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                                randomSeed = 1234,
                                verbose = 3) 

bwnet_clat_h2 = blockwiseModules(norma.counts_clat, 
                                 maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                                 TOMType = "signed",
                                 networkType = "signed",
                                 power = soft_power_clat_h2,
                                 mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                                 numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                                 randomSeed = 1234,
                                 verbose = 3)

cor = temp_cor

module_colour_cre_h2 = bwnet_cre_h2$colors
module_colour_clat_h2 = bwnet_clat_h2$colors

module_eigengenes_cre = bwnet_cre_h2$MEs
module_eigengenes_clat = bwnet_clat_h2$MEs
#module_colour_cre = labels2colors(module_colour_cre)
sex.cre = binarizeCategoricalColumns(coldata2$sex[1:12],
                                     includePairwise = TRUE,
                                     includeLevelVsAll = FALSE)

row.names(sex.cre) = row.names(coldata2[1:12, ])

sex.clat = binarizeCategoricalColumns(coldata2$sex[16:27],
                                      includePairwise = TRUE,
                                      includeLevelVsAll = FALSE)

row.names(sex.clat) = row.names(coldata2[16:27, ])


##2. Tissues
tissues_cre = binarizeCategoricalVariable(coldata2$tissue[1:12],
                                          includePairwise = TRUE,
                                          includeLevelVsAll = FALSE) ##o is the reference category
row.names(tissues_cre) = row.names(coldata2[1:12, ]) #need to change rownames when binarizing categorical variables


tissue.clat = binarizeCategoricalColumns(coldata2$tissue[16:27],
                                         includePairwise = TRUE,
                                         includeLevelVsAll = FALSE)

row.names(tissue.clat) = row.names(coldata2[16:27, ])
traits_cre = cbind(sex.cre, tissues_cre)

traits_clat = cbind(sex.clat, tissue.clat)
##define number of genes and samples

nSample_cre = nrow(norma.counts_cre)
nGenes_cre = ncol(norma.counts_cre)

nSample_clat = nrow(norma.counts_clat)
nGenes_clat = ncol(norma.counts_clat)


##Extra steps if module number is used instead of module colour
# MEs0 = moduleEigengenes(norma.counts, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)

##Cre
modules.trait.correlation.sex.cre = cor(module_eigengenes_cre, sex.cre, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue.cre = cor(module_eigengenes_cre, tissues_cre, use = "p") #correlating eingengenes adn traits

modules.trait.correlation.cre = cor(module_eigengenes_cre, traits_cre, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sex.cre = corPvalueStudent(modules.trait.correlation.sex.cre, nSample_cre)
modules.trait.corr.pvals.tissue.cre = corPvalueStudent(modules.trait.correlation.tissue.cre, nSample_cre)

modules.trait.corr.pvals.cre = corPvalueStudent(modules.trait.correlation.cre, nSample_cre)

#clat

modules.trait.correlation.sex.clat = cor(module_eigengenes_clat, sex.clat, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue.clat = cor(module_eigengenes_clat, tissue.clat, use = "p") #correlating eingengenes adn traits

modules.trait.correlation.clat = cor(module_eigengenes_clat, traits_clat, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sex.clat = corPvalueStudent(modules.trait.correlation.sex.clat, nSample_clat)
modules.trait.corr.pvals.tissue.clat = corPvalueStudent(modules.trait.correlation.tissue.clat, nSample_clat)

modules.trait.corr.pvals.clat = corPvalueStudent(modules.trait.correlation.clat, nSample_clat)

##############################################################################################################
#visualize module trait association as a heatmap
##############################################################################################################

heatmap_data_cre = merge(module_eigengenes_cre, traits_cre, by = "row.names")
heatmap_data_cre = heatmap_data_cre %>% 
  column_to_rownames(var = "Row.names")

colnames(heatmap_data_cre)
CorLevelPlot(heatmap_data_cre,
             x = names(heatmap_data_cre)[15:17], 
             y = names(heatmap_data_cre)[1:14],
             col = c("red", "pink","white", "pink", "red"),
             main = "SEX")



CorLevelPlot(heatmap_data_cre,
             x = names(heatmap_data_cre)[15], 
             y = names(heatmap_data_cre)[1:12],
             col = c("#004b8d","skyblue", "white", "skyblue", "#004b8d" ),
             main = "TISSUES")

CorLevelPlot(heatmap_data_cre,
             x = names(heatmap_data_cre)[13:14], 
             y = names(heatmap_data_cre)[1:12],
             col = c("blue1", "skyblue", "white", "pink", "red" ))



##CLAT

heatmap_data_clat = merge(module_eigengenes_clat, traits_clat, by = "row.names")
heatmap_data_clat = heatmap_data_clat %>% 
  column_to_rownames(var = "Row.names")
colnames(heatmap_data_clat)

CorLevelPlot(heatmap_data_clat,
             x = names(heatmap_data_clat)[19:20], 
             y = names(heatmap_data_clat)[1:18],
             col = c("blue1", "skyblue", "white", "pink", "red" ))

###Modules in red/blue are significantly associated with one of the sex/species over the other
colnames(heatmap_data)

options(scipen=999)
?CorLevelPlot

cor.test(heatmap_data_cre$data.M.vs.F, heatmap_data_cre$MEbrown, method = "pearson", alternative = "two.sided")

#gives me pariwise correalation 
View(cor(heatmap_data_cre, method = "pearson"))

#correlation values that are 
correaltion_cre = as.data.frame(round(cor(heatmap_data_cre, method = "pearson")[(1:12), (13:14)], digits = 2))
correaltion_cre = rownames_to_column(correaltion_cre)
correaltion_cre = correaltion_cre[order(correaltion_cre$rowname), ]

correaltion_clat = as.data.frame(round(cor(heatmap_data_clat, method = "pearson")[(1:18), (19:20)], digits = 2))
correaltion_clat = rownames_to_column(correaltion_clat)
correaltion_clat = correaltion_clat[order(correaltion_clat$rowname), ]

##Combining correlation adn module size

#after running module preservation
Z.PreservationStats_a #cre x h2
Z.PreservationStats_b #crexclat

Z.PreservationStats_a1 = rownames_to_column(Z.PreservationStats_a)
Z.PreservationStats_a1 = Z.PreservationStats_a1[order(Z.PreservationStats_a1$rowname), ][-4, (1:3)]
Z.PreservationStats_a1 = Z.PreservationStats_a1 %>% rename(modulesize_crexh2 = moduleSize, Zsummary_crexh2 = Zsummary.pres)

Z.PreservationStats_b1 = rownames_to_column(Z.PreservationStats_b)
Z.PreservationStats_b1 = Z.PreservationStats_b1[order(Z.PreservationStats_b1$rowname), ][-4, (1:3)]
Z.PreservationStats_b1 = Z.PreservationStats_b1 %>% rename(modulesize_crexclat = moduleSize, Zsummary_crexclat = Zsummary.pres)


Z.PreservationStats_c1 = rownames_to_column(Z.PreservationStats_c) 
Z.PreservationStats_c1 = Z.PreservationStats_c1[order(Z.PreservationStats_c1$rowname), ][-5, (1:3)]
Z.PreservationStats_c1 = Z.PreservationStats_c1 %>% rename(modulesize_clatxh2 = moduleSize, Zsummary_clatxh2 = Zsummary.pres)

Z.PreservationStats_d1 = rownames_to_column(Z.PreservationStats_d)
Z.PreservationStats_d1 = Z.PreservationStats_d1[order(Z.PreservationStats_d1$rowname), ][-5, (1:3)]
Z.PreservationStats_d1 = Z.PreservationStats_d1 %>% rename(modulesize_clatxcre = moduleSize, Zsummary_clatxcre = Zsummary.pres)

cumulative_cre_h2 = cbind(correaltion_cre, Z.PreservationStats_b1, Z.PreservationStats_a1) 

cumulative_clat_h2 = cbind(correaltion_clat, Z.PreservationStats_c1, Z.PreservationStats_d1) 


#cumulative_cre_h2 = (cumulative_cre_h2 %>% gather(avsb, values, data.M.vs.F, S.vs.G, Zsummary_crexclat, Zsummary_crexh2, modulesize_crexclat))
cumulative_cre_h2 = cumulative_cre_h2[, -(1:2)]

# write.csv(cumulative_cre_h2, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/module_preservation_h2_cre.csv",
#            row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.csv(cumulative_clat_h2, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/module_preservation_h2_clat.csv",
#           row.names = FALSE, col.names = FALSE, quote = FALSE)

# ggplot(cumulative_cre_h2, aes(cumulative_cre_h2$avsb, cumulative_cre_h2$rowname)) + 
#   geom_tile()+
#   geom_text(aes(label = values), color = "white", size = 4)

# We now set up the multi-set expression data
# and corresponding module colors:
setLabels = c("Cre", "H2", "Clat")
multiExpr_h2=list(Cre=list(data=norma.counts_cre),
                  H2=list(data=norma.counts_h2),
                  Clat = list(data=norma.counts_clat))

multiColor_h2=list(Cre=module_colour_cre_h2, Clat=module_colour_clat_h2)


# set the random seed of the permutation test analysis
set.seed(1)
system.time({
  mp_h2 = modulePreservation(multiExpr_h2, multiColor_h2,
                             referenceNetworks = c(1, 3), nPermutations = 200, networkType = "signed",
                             randomSeed = 1, quickCor = 0, verbose = 3)
})


# # Save the results of the module preservation analysis
#save(mp, file = "Allele Specific Expression/Hybrid Analysis/modulePreservation_cre_h2.RData")
# # If needed, reload the data:
# load(file = "modulePreservation_cre_h2.RData")

# specify the reference and the test networks
#crexh2
ref=1; test = 2

str(mp)
Obs.PreservationStats_a= mp_h2$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_a=mp_h2$preservation$Z[[ref]][[test]]

#crexclat
ref=1; test = 3

str(mp)
Obs.PreservationStats_b= mp_h2$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_b=mp_h2$preservation$Z[[ref]][[test]]

#Clat as reference set
#clat x h2
ref=2; test = 2

Obs.PreservationStats_c= mp_h2$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_c=mp_h2$preservation$Z[[ref]][[test]]

#clat x cre
ref=2; test = 1

Obs.PreservationStats_d= mp_h2$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_d=mp_h2$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats_a
Obs.PreservationStats_b
Obs.PreservationStats_c
Obs.PreservationStats_d

# Z statistics from the permutation test analysis
# View(Z.PreservationStats)
# View(mp$preservation)
# View(mp)
# str(Obs.PreservationStats)

##DATA VISUALIZATION
#cre as ref
modColors_a = rownames(Obs.PreservationStats_a) #cre x h2
moduleSize_a = Obs.PreservationStats_a$moduleSize

modColors_b = rownames(Obs.PreservationStats_b) #cre x clat
moduleSize_b = Obs.PreservationStats_b$moduleSize

#clat as ref

modColors_c = rownames(Obs.PreservationStats_c) #clat x h2
moduleSize_c = Obs.PreservationStats_c$moduleSize

modColors_d = rownames(Obs.PreservationStats_d) #clat x cre
moduleSize_d = Obs.PreservationStats_d$moduleSize

# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label_a = modColors_a[selectModules]
point.label_a0= modColors_a

point.label_b = modColors_b[selectModules]
point.label_b0= modColors_b

point.label_c = modColors_c[selectModules]
point.label_c0= modColors_c

point.label_d = modColors_d[selectModules]
point.label_d0= modColors_d

#Composite preservation statistics
medianRank_a=Obs.PreservationStats_a$medianRank.pres
Zsummary_a=Z.PreservationStats_a$Zsummary.pres

medianRank_b=Obs.PreservationStats_b$medianRank.pres
Zsummary_b=Z.PreservationStats_b$Zsummary.pres

medianRank_c=Obs.PreservationStats_c$medianRank.pres
Zsummary_c=Z.PreservationStats_c$Zsummary.pres

medianRank_d=Obs.PreservationStats_d$medianRank.pres
Zsummary_d=Z.PreservationStats_d$Zsummary.pres

# par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# # plot medianRank versus module size
# plot(moduleSize[selectModules],medianRank[selectModules],col=1,
#      bg=modColors[selectModules],pch = 21,main="medianRank Preservation",
#      cex = 2, ylab ="medianRank",xlab="Module size", log="x");
# labelPoints(moduleSize[selectModules],medianRank[selectModules], point.label, cex=1, offs=0.03)

##cre x clat

# plot(moduleSize_a[selectModules],medianRank_a[selectModules],col=1,
#      bg=modColors_a[selectModules],pch = 21,main="medianRank Preservation",
#      cex = 2, ylab ="medianRank",xlab="Module size", log="x");
# labelPoints(moduleSize_a[selectModules],medianRank_a[selectModules], point.label_a, cex=1, offs=0.03)

combined_cre_h2 = cbind(Obs.PreservationStats_a, Z.PreservationStats_a, point.label_a0)
combined_cre_h2 = combined_cre_h2[ ,-1]

combined_cre_clat_2 = cbind(Obs.PreservationStats_b, Z.PreservationStats_b, point.label_b0)
combined_cre_clat_2 = combined_cre_clat_2[ ,-1]

combined_clat_h2 = cbind(Obs.PreservationStats_c, Z.PreservationStats_c, point.label_c0)
combined_clat_h2 = combined_clat_h2[ ,-1]

combined_clat_cre_2 = cbind(Obs.PreservationStats_d, Z.PreservationStats_d, point.label_d0)
combined_clat_cre_2 = combined_clat_cre_2[ ,-1]


##median rank vs mdoule size using ggplot
mrp_cre_h2 = ggplot(combined_cre_h2, aes(x=combined_cre_h2$moduleSize, y=combined_cre_h2$medianRank.pres, colour = point.label_a0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_a0)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Cre x H2")+
  theme_bw()
# geom_hline(yintercept=0, linetype=2) +
# geom_vline(xintercept=0, linetype=2) +

mrp_cre_clat_2 = ggplot(combined_cre_clat_2, aes(x=combined_cre_clat_2$moduleSize, y=combined_cre_clat_2$medianRank.pres, colour = point.label_b0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_b0)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Cre (ref) x Clat (H2)")+
  theme_bw()
# geom_hline(yintercept=0, linetype=2) +
# geom_vline(xintercept=0, linetype=2) +

mrp_clat_h2 = ggplot(combined_clat_h2, aes(x=combined_clat_h2$moduleSize, y=combined_clat_h2$medianRank.pres, colour = point.label_c0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_c0)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Clat x H2")+
  theme_bw()

mrp_clat_cre_2 = ggplot(combined_clat_cre_2, aes(x=combined_clat_cre_2$moduleSize, y=combined_clat_cre_2$medianRank.pres, colour = point.label_d0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_d0)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Clat (ref) x Cre (H2)")+
  theme_bw()


medianrank_vs_modulesize_h2 = grid.arrange(mrp_cre_clat_2, mrp_clat_cre_2, mrp_cre_h2, mrp_clat_h2, nrow = 2)

ggsave(medianrank_vs_modulesize_h2, filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/medianrank_vs_modulesize_h2.pdf", dpi = 300, width = 18, height = 12 )



# plot Zsummary versus module size
# plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
#      bg=modColors[selectModules],pch = 21,main="Zsummary Preservation",
#      cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
# labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
# # Add threshold lines for Zsummary
# abline(h=0);abline(h=2, col = "blue", lty = 2)
# abline(h=10, col = "red", lty = 2)


#Cre x clat
# plot(moduleSize_a[selectModules],Zsummary_a[selectModules], col = 1,
#      bg=modColors_a[selectModules],pch = 21,main="Zsummary Preservation",
#      cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
# labelPoints(moduleSize_a[selectModules],Zsummary_a[selectModules],point.label_a,cex=1,offs=0.03)
# # Add threshold lines for Zsummary
# abline(h=0);abline(h=2, col = "blue", lty = 2)
# abline(h=10, col = "red", lty = 2)

# Caption: Preservation of female mouse liver modules in male livers. The right panel shows the composite statistic medianRank (Eq.9.20) versus module size. 
# The higher the medianRank the less preserved is the module relative to other modules. Since medianRank is based on the observed preservation statistics 
# (as opposed to Z statistics or p-values) we find that it is much less dependent on module size. The upper right panel shows the composite statistic Zsummary (Eq.9.1). 
# If Zsummary> 10 there is strong evidence that the module is preserved (Langfelder et al 2011). If Zsummary<2, there is no evidence that the module preserved. 
# Note that Zsummary shows a strong dependence on module size. The lightyellow female module shows no evidence of preservation in the male liver network.
zsum_cre_h2 = ggplot(combined_cre_h2, aes(x=combined_cre_h2$moduleSize, y=combined_cre_h2$Zsummary.pres, colour = point.label_a0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_a0)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Cre x H2")+
  geom_hline(yintercept=2, linetype=2, colour = "red")

zsum_cre_clat_2 = ggplot(combined_cre_clat_2, aes(x=combined_cre_clat_2$moduleSize, y=combined_cre_clat_2$Zsummary.pres, colour = point.label_b0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_b0)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Cre (ref) x Clat (H2)")+
  geom_hline(yintercept=2, linetype=2, colour = "red")



zsum_clat_h2 = ggplot(combined_clat_h2, aes(x=combined_clat_h2$moduleSize, y=combined_clat_h2$Zsummary.pres, colour = point.label_c0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_c0)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Clat x H2")+
  geom_hline(yintercept=2, linetype=2, colour = "red")

zsum_clat_cre_2 = ggplot(combined_clat_cre_2, aes(x=combined_clat_cre_2$moduleSize, y=combined_clat_cre_2$Zsummary.pres, colour = point.label_d0)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_d0)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Clat (ref) x Cre (H2)")+
  geom_hline(yintercept=2, linetype=2, colour = "red")

zsummary_vs_modulesize_h2 = grid.arrange(zsum_cre_clat_2, zsum_clat_cre_2, zsum_cre_h2, zsum_clat_h2, nrow = 2)

ggsave(zsummary_vs_modulesize_h2, filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/zsummary_vs_modulesize_h2.pdf", dpi = 300, width = 18, height = 12 )



############################################################################################################################################################
####H1 MODULE PRESERVATION ANALYSIS##########################################################################################################################
#############################################################################################################################################################

###Getting normalized count for each data separately
norma.counts_cre=norma.counts[1:15, ]
norma.counts_cre_1 = norma.counts_cre[-(7:9), ] #removing MG data
norma.counts_cre_1 = norma.counts_cre_1[-(10:12), ]
norma.counts_clat=norma.counts[16:30, ]
norma.counts_clat_1 = norma.counts_clat[-(7:9), ] #removing MG data
norma.counts_clat_1 = norma.counts_clat_1[-(10:12), ]
norma.counts_h1=norma.counts[31:42, ]
norma.counts_h1=norma.counts_h1[-(7:9), ]

rownames(norma.counts_h1)
rownames(norma.counts_clat_1)
rownames(norma.counts_cre_1)


# We now set up the multi-set expression data

setLabels = c("Cre", "H1", "Clat")
multiExpr=list(Cre=list(data=norma.counts_cre_1),
               H1=list(data=norma.counts_h1),
               Clat = list(data=norma.counts_clat_1))

View(multiExpr)
###QUALITY CONTROL TO DETECT OUTLIERS FROM THE MULITSET DATA

###checking samples for mutliExpr data

gsg_t = goodSamplesGenesMS(multiExpr, verbose = 3)

gsg_t$allOK

table(gsg_t$goodGenes) #number of genes that are outliers = FALSE
table(gsg_t$goodSamples)

##Filtering base genes from both datasets
norma.counts_h1 = norma.counts_h1[,gsg_t$goodGenes == TRUE] ##filtering good genes
norma.counts_cre_1 = norma.counts_cre_1[,gsg_t$goodGenes == TRUE]
norma.counts_clat_1 = norma.counts_clat_1[,gsg_t$goodGenes == TRUE]


##remaking the data with good genes
multiExpr=list(Cre=list(data=norma.counts_cre_1),
               H1=list(data=norma.counts_h1),
               Clat = list(data=norma.counts_clat_1))



#####Network construction

##Choose a set of soft threshold powers

power = c(c(1:10), seq(from = 12, to = 50, by = 2))

sft_cre_h1 = pickSoftThreshold(norma.counts_cre_1,
                               powerVector = power,
                               networkType = "signed",
                               verbose = 5)

sft_clat_h1 = pickSoftThreshold(norma.counts_clat_1,
                                powerVector = power,
                                networkType = "signed",
                                verbose = 5)

sft.data_cre_h1 = sft_cre_h1$fitIndices ##we will use the R suqared value and mean connectivity , max R^2 and min mean connectivity
sft.data_clat_h1 = sft_clat_h1$fitIndices


###Visualization to pick the right power
a1_cre_h1 = ggplot(sft.data_cre_h1, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2_cre_h1 = ggplot(sft.data_clat_h1, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()


grid.arrange(a1_cre_h1, a2_cre_h1, nrow = 2)  ##30 with WT and without wt

a1_clat_h1 = ggplot(sft.data_clat_h1, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2_clat_h1 = ggplot(sft.data_clat_h1, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1_clat_h1, a2_clat_h1, nrow = 2) #28 with WT and 30


###converting matric to numeric

norma.counts_cre_1[] = sapply(norma.counts_cre_1, as.numeric)  ##different step than official tutorial
norma.counts_clat_1[] = sapply(norma.counts_clat_1, as.numeric)  ##different step than official tutorial

soft_power_cre_h1 = 30
soft_power_clat_h1 = 30

temp_cor = cor #to prevent WGCNa from using other cor function

cor = WGCNA::cor


bwnet_cre_h1 = blockwiseModules(norma.counts_cre_1, 
                                maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                                TOMType = "signed",
                                networkType = "signed",
                                power = soft_power_cre_h1,
                                mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                                numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                                randomSeed = 1234,
                                verbose = 3) 

bwnet_clat_h1 = blockwiseModules(norma.counts_clat_1, 
                                 maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                                 TOMType = "signed",
                                 networkType = "signed",
                                 power = soft_power_clat_h1,
                                 mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                                 numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                                 randomSeed = 1234,
                                 verbose = 3)

cor = temp_cor

module_colour_cre_h1 = bwnet_cre_h1$colors
module_colour_clat_h1 = bwnet_clat_h1$colors
#module_colour_cre = labels2colors(module_colour_cre)

module_eigengenes_cre_h1 = bwnet_cre_h1$MEs
module_eigengenes_clat_h1 = bwnet_clat_h1$MEs
#module_colour_cre = labels2colors(module_colour_cre)

coldata3 = rbind(coldata2[1:12, ], coldata2[16:27, ])
coldata3 = coldata3[ -c(7,8,9,19,20,21), ]
sex.cre_h1 = binarizeCategoricalColumns(coldata3$sex[1:9],
                                     includePairwise = TRUE,
                                     includeLevelVsAll = FALSE)

row.names(sex.cre_h1) = row.names(coldata3[1:9, ])

sex.clat_h1 = binarizeCategoricalColumns(coldata3$sex[10:18],
                                      includePairwise = TRUE,
                                      includeLevelVsAll = FALSE)

row.names(sex.clat_h1) = row.names(coldata3[10:18, ])


##2. Tissues
tissues_cre_h1 = binarizeCategoricalVariable(coldata3$tissue[1:9],
                                          includePairwise = TRUE,
                                          includeLevelVsAll = FALSE) ##o is the reference category
row.names(tissues_cre_h1) = row.names(coldata3[1:9, ]) #need to change rownames when binarizing categorical variables


tissue.clat_h1 = binarizeCategoricalColumns(coldata3$tissue[10:18],
                                         includePairwise = TRUE,
                                         includeLevelVsAll = FALSE)

row.names(tissue.clat_h1) = row.names(coldata3[10:18, ])
traits_cre_h1 = cbind(sex.cre_h1, tissues_cre_h1)

traits_clat_h1 = cbind(sex.clat_h1, tissue.clat_h1)
##define number of genes and samples

nSample_cre_h1 = nrow(norma.counts_cre_1)
nGenes_cre_h1 = ncol(norma.counts_cre_1)

nSample_clat_h1 = nrow(norma.counts_clat_1)
nGenes_clat_h1 = ncol(norma.counts_clat_1)


##Extra steps if module number is used instead of module colour
# MEs0 = moduleEigengenes(norma.counts, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)

##Cre
modules.trait.correlation.sex.cre_h1 = cor(module_eigengenes_cre_h1, sex.cre_h1, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue.cre_h1 = cor(module_eigengenes_cre_h1, tissues_cre_h1, use = "p") #correlating eingengenes adn traits

modules.trait.correlation.cre_h1 = cor(module_eigengenes_cre_h1, traits_cre_h1, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sex.cre_h1 = corPvalueStudent(modules.trait.correlation.sex.cre_h1, nSample_cre_h1)
modules.trait.corr.pvals.tissue.cre_h1 = corPvalueStudent(modules.trait.correlation.tissue.cre_h1, nSample_cre_h1)

modules.trait.corr.pvals.cre_h1 = corPvalueStudent(modules.trait.correlation.cre_h1, nSample_cre_h1)

#clat

modules.trait.correlation.sex.clat_h1 = cor(module_eigengenes_clat_h1, sex.clat_h1, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue.clat_h1 = cor(module_eigengenes_clat_h1, tissue.clat_h1, use = "p") #correlating eingengenes adn traits

modules.trait.correlation.clat_h1 = cor(module_eigengenes_clat_h1, traits_clat_h1, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sex.clat_h1 = corPvalueStudent(modules.trait.correlation.sex.clat_h1, nSample_clat_h1)
modules.trait.corr.pvals.tissue.clat_h1 = corPvalueStudent(modules.trait.correlation.tissue.clat_h1, nSample_clat_h1)

modules.trait.corr.pvals.clat_h1 = corPvalueStudent(modules.trait.correlation.clat_h1, nSample_clat_h1)

##############################################################################################################
#visualize module trait association as a heatmap
##############################################################################################################

heatmap_data_cre_h1 = merge(module_eigengenes_cre_h1, traits_cre_h1, by = "row.names")
heatmap_data_cre_h1 = heatmap_data_cre_h1 %>% 
  column_to_rownames(var = "Row.names")

colnames(heatmap_data_cre_h1)

CorLevelPlot(heatmap_data_cre_h1,
             x = names(heatmap_data_cre_h1)[12:13], 
             y = names(heatmap_data_cre_h1)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red" ))



##CLAT

heatmap_data_clat_h1 = merge(module_eigengenes_clat_h1, traits_clat_h1, by = "row.names")
heatmap_data_clat_h1 = heatmap_data_clat_h1 %>% 
  column_to_rownames(var = "Row.names")
colnames(heatmap_data_clat_h1)

CorLevelPlot(heatmap_data_clat_h1,
             x = names(heatmap_data_clat_h1)[15:16], 
             y = names(heatmap_data_clat_h1)[1:14],
             col = c("blue1", "skyblue", "white", "pink", "red" ))

###Modules in red/blue are significantly associated with one of the sex/species over the other


options(scipen=999)


#correlation values that are 
correaltion_cre_h1 = as.data.frame(round(cor(heatmap_data_cre_h1, method = "pearson")[(1:11), (12:13)], digits = 2))
correaltion_cre_h1 = rownames_to_column(correaltion_cre_h1)
correaltion_cre_h1 = correaltion_cre_h1[order(correaltion_cre_h1$rowname), ]

correaltion_clat_h1 = as.data.frame(round(cor(heatmap_data_clat_h1, method = "pearson")[(1:14), (15:16)], digits = 2))
correaltion_clat_h1 = rownames_to_column(correaltion_clat_h1)
correaltion_clat_h1 = correaltion_clat_h1[order(correaltion_clat_h1$rowname), ]

##Combining correlation adn module size

#after running module preservation
Z.PreservationStats_1 #cre x h1
Z.PreservationStats_2 #crexclat

Z.PreservationStats_11 = rownames_to_column(Z.PreservationStats_1)
Z.PreservationStats_11 = Z.PreservationStats_11[order(Z.PreservationStats_11$rowname), ][-4, (1:3)]
Z.PreservationStats_11 = Z.PreservationStats_11 %>% rename(modulesize_crexh1 = moduleSize, Zsummary_crexh1 = Zsummary.pres)

Z.PreservationStats_21 = rownames_to_column(Z.PreservationStats_2)
Z.PreservationStats_21 = Z.PreservationStats_21[order(Z.PreservationStats_21$rowname), ][-4, (1:3)]
Z.PreservationStats_21 = Z.PreservationStats_21 %>% rename(modulesize_crexclat = moduleSize, Zsummary_crexclat = Zsummary.pres)


Z.PreservationStats_31 = rownames_to_column(Z.PreservationStats_3) 
Z.PreservationStats_31 = Z.PreservationStats_31[order(Z.PreservationStats_31$rowname), ][-4, (1:3)]
Z.PreservationStats_31 = Z.PreservationStats_31 %>% rename(modulesize_clatxh1 = moduleSize, Zsummary_clatxh1 = Zsummary.pres)

Z.PreservationStats_41 = rownames_to_column(Z.PreservationStats_4)
Z.PreservationStats_41 = Z.PreservationStats_41[order(Z.PreservationStats_41$rowname), ][-4, (1:3)]
Z.PreservationStats_41 = Z.PreservationStats_41 %>% rename(modulesize_clatxcre = moduleSize, Zsummary_clatxcre = Zsummary.pres)

cumulative_cre_h1 = cbind(correaltion_cre_h1, Z.PreservationStats_21, Z.PreservationStats_11) 

cumulative_clat_h1 = cbind(correaltion_clat_h1, Z.PreservationStats_31, Z.PreservationStats_41) 


# write.csv(cumulative_cre_h1, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/module_preservation_h1_cre.csv",
#            row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.csv(cumulative_clat_h1, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/module_preservation_h1_clat.csv",
#           row.names = FALSE, col.names = FALSE, quote = FALSE)









###############################################
# We now set up the multi-set expression data
# and corresponding module colors:
setLabels = c("Cre", "H1", "Clat")
multiExpr_h1=list(Cre=list(data=norma.counts_cre_1),
                  H1=list(data=norma.counts_h1),
                  Clat = list(data=norma.counts_clat_1))

multiColor_h1=list(Cre=module_colour_cre_h1, Clat=module_colour_clat_h1)

# set the random seed of the permutation test analysis
set.seed(1)
system.time({
  mp_h1 = modulePreservation(multiExpr_h1, multiColor_h1,
                             referenceNetworks = c(1, 3), nPermutations = 200, networkType = "signed",
                             randomSeed = 1, quickCor = 0, verbose = 3)
})

# specify the reference and the test networks
#crexh1
ref=1; test = 2

Obs.PreservationStats_1= mp_h1$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_1=mp_h1$preservation$Z[[ref]][[test]]

#crexclat
ref=1; test = 3

Obs.PreservationStats_2= mp_h1$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_2=mp_h1$preservation$Z[[ref]][[test]]

#Clat as reference set
#clat x h1
ref=2; test = 2

Obs.PreservationStats_3= mp_h1$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_3=mp_h1$preservation$Z[[ref]][[test]]

#clat x cre
ref=2; test = 1

Obs.PreservationStats_4= mp_h1$preservation$observed[[ref]][[test]] ##ref is first element of the list, test is element of the chosen list
Z.PreservationStats_4=mp_h1$preservation$Z[[ref]][[test]]

# Look at the observed preservation statistics
Obs.PreservationStats_1
Obs.PreservationStats_2
Obs.PreservationStats_3
Obs.PreservationStats_4

##DATA VISUALIZATION
#cre as ref
modColors_1 = rownames(Obs.PreservationStats_1) #cre x h1
moduleSize_1 = Obs.PreservationStats_1$moduleSize

modColors_2 = rownames(Obs.PreservationStats_2) #cre x clat
moduleSize_2 = Obs.PreservationStats_2$moduleSize

#clat as ref
modColors_3 = rownames(Obs.PreservationStats_3) #clat x h1
moduleSize_3 = Obs.PreservationStats_3$moduleSize

modColors_4 = rownames(Obs.PreservationStats_4) #clat x cre
moduleSize_4 = Obs.PreservationStats_4$moduleSize

# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label_1 = modColors_1[selectModules]
point.label_10= modColors_1

point.label_2 = modColors_2[selectModules]
point.label_20= modColors_2

point.label_3 = modColors_3[selectModules]
point.label_30= modColors_3

point.label_4 = modColors_4[selectModules]
point.label_40= modColors_4

#Composite preservation statistics
medianRank_1=Obs.PreservationStats_1$medianRank.pres
Zsummary_1=Z.PreservationStats_1$Zsummary.pres

medianRank_2=Obs.PreservationStats_2$medianRank.pres
Zsummary_2=Z.PreservationStats_2$Zsummary.pres

medianRank_3=Obs.PreservationStats_3$medianRank.pres
Zsummary_3=Z.PreservationStats_3$Zsummary.pres

medianRank_4=Obs.PreservationStats_4$medianRank.pres
Zsummary_4=Z.PreservationStats_4$Zsummary.pres


###Combining data for ggplot

combined_cre_h1 = cbind(Obs.PreservationStats_1, Z.PreservationStats_1, point.label_10)
combined_cre_h1 = combined_cre_h1[ ,-1]

combined_cre_clat_1 = cbind(Obs.PreservationStats_2, Z.PreservationStats_2, point.label_20)
combined_cre_clat_1 = combined_cre_clat_1[ ,-1]

combined_clat_h1 = cbind(Obs.PreservationStats_3, Z.PreservationStats_3, point.label_30)
combined_clat_h1 = combined_clat_h1[ ,-1]

combined_clat_cre_1 = cbind(Obs.PreservationStats_4, Z.PreservationStats_4, point.label_40)
combined_clat_cre_1 = combined_clat_cre_1[ ,-1]

##median rank vs mdoule size using ggplot
mrp_cre_h1 = ggplot(combined_cre_h1, aes(x=combined_cre_h1$moduleSize, y=combined_cre_h1$medianRank.pres, colour = point.label_10)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_10)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Cre x H1")+
  theme_bw()
# geom_hline(yintercept=0, linetype=2) +
# geom_vline(xintercept=0, linetype=2) +

mrp_cre_clat_1 = ggplot(combined_cre_clat_1, aes(x=combined_cre_clat_1$moduleSize, y=combined_cre_clat_1$medianRank.pres, colour = point.label_20)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_20)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Cre (ref) x Clat (H1)")+
  theme_bw()
# geom_hline(yintercept=0, linetype=2) +
# geom_vline(xintercept=0, linetype=2) +

mrp_clat_h1 = ggplot(combined_clat_h1, aes(x=combined_clat_h1$moduleSize, y=combined_clat_h1$medianRank.pres, colour = point.label_30)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_30)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Clat x H1")+
  theme_bw()

mrp_clat_cre_1 = ggplot(combined_clat_cre_1, aes(x=combined_clat_cre_1$moduleSize, y=combined_clat_cre_1$medianRank.pres, colour = point.label_40)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 20)+
  geom_text(label = point.label_40)+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  geom_hline(yintercept=20, linetype=2, color = "green") +
  ggtitle("Clat (ref) x Cre (H1)")+
  theme_bw()


medianrank_vs_modulesize_h1 = grid.arrange(mrp_cre_clat_1, mrp_clat_cre_1, mrp_cre_h1, mrp_clat_h1, nrow = 2)

ggsave(medianrank_vs_modulesize_h1, filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/medianrank_vs_modulesize_h1.pdf", dpi = 300, width = 18, height = 12 )


####Z-SUMMARY VS MODULE SIZE
zsum_cre_h1 = ggplot(combined_cre_h1, aes(x=combined_cre_h1$moduleSize, y=combined_cre_h1$Zsummary.pres, colour = point.label_10)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_10)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Cre x H1")+
  geom_hline(yintercept=2, linetype=2, colour = "red")

zsum_cre_clat_1 = ggplot(combined_cre_clat_1, aes(x=combined_cre_clat_1$moduleSize, y=combined_cre_clat_1$Zsummary.pres, colour = point.label_20)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_20)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Cre (ref) x Clat (H1)")+
  geom_hline(yintercept=2, linetype=2, colour = "red")


zsum_clat_h1 = ggplot(combined_clat_h1, aes(x=combined_clat_h1$moduleSize, y=combined_clat_h1$Zsummary.pres, colour = point.label_30)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_30)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Clat x H1")+
  geom_hline(yintercept=2, linetype=2, colour = "red")

zsum_clat_cre_1 = ggplot(combined_clat_cre_1, aes(x=combined_clat_cre_1$moduleSize, y=combined_clat_cre_1$Zsummary.pres, colour = point.label_40)) +
  geom_point()+
  xlim(0, 1500) +
  ylim(0, 60)+
  geom_text(label = point.label_40)+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_continuous(breaks=seq(0,1500,250))+
  theme_bw()+
  geom_hline(yintercept=10, linetype=2, color = "blue") +
  geom_hline(yintercept=60, linetype=2, color = "green") +
  geom_vline(xintercept=0, linetype=2, color = "green") +
  ggtitle("Clat (ref) x Cre (H1)")+
  geom_hline(yintercept=2, linetype=2, colour = "red")

zsummary_vs_modulesize_h1 = grid.arrange(zsum_cre_clat_1, zsum_clat_cre_1, zsum_cre_h1, zsum_clat_h1, nrow = 2)

ggsave(zsummary_vs_modulesize_h1, filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 6/zsummary_vs_modulesize_h1.pdf", dpi = 300, width = 18, height = 12 )





