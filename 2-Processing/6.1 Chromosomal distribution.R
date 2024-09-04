#first run 2.classifying inheritance.R and 3. Allele specific expression.R to get the categories

##Chromosomal distribution of gene expression inheritance pattern

orthologs_chr = read.table(file = "orthologs_chr.txt", sep = "\t", header = TRUE)
orthologs_chr = orthologs_chr[order(orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name), ]


###FG
head(dds.fg.h1)
dds.fg.h1_chr = na.omit(cbind(dds.fg.h1[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
           orthologs_chr$chromosome))

#write.csv(dds.fg.h1_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fg_h1_inht_chr.csv")


dds.fg.h1_chr_count = as.data.frame(table(dds.fg.h1_chr$inheritance.fg, dds.fg.h1_chr$`orthologs_chr$chromosome`)) 
dds.fg.h1_chr_count = dds.fg.h1_chr_count %>% rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.fg.h1_chr_count$Sample = "H1"
dds.fg.h1_chr_count$Tissue = "G"
dds.fg.h1_chr_count$Sex = "Female"
dds.fg.h1_chr_count$Name = paste(dds.fg.h1_chr_count$Sample, dds.fg.h1_chr_count$Sex, dds.fg.h1_chr_count$Tissue, dds.fg.h1_chr_count$Chromosome)


dds.fg.h2_chr = na.omit(cbind(dds.fg.h2[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.fg.h2_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fg_h2_inht_chr.csv")

dds.fg.h2_chr_count = as.data.frame(table(dds.fg.h2_chr$inheritance.fg, dds.fg.h2_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.fg.h2_chr_count$Sample = "H2"
dds.fg.h2_chr_count$Tissue = "G"
dds.fg.h2_chr_count$Sex = "Female"
dds.fg.h2_chr_count$Name = paste(dds.fg.h2_chr_count$Sample, dds.fg.h2_chr_count$Sex, dds.fg.h2_chr_count$Tissue, dds.fg.h2_chr_count$Chromosome)


###FS
dds.fs.h1_chr = na.omit(cbind(dds.fs.h1[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.fs.h1_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fs_h1_inht_chr.csv")

dds.fs.h1_chr_count = as.data.frame(table(dds.fs.h1_chr$inheritance.fs, dds.fs.h1_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.fs.h1_chr_count$Sample = "H1"
dds.fs.h1_chr_count$Tissue = "S"
dds.fs.h1_chr_count$Sex = "Female"
dds.fs.h1_chr_count$Name = paste(dds.fs.h1_chr_count$Sample, dds.fs.h1_chr_count$Sex, dds.fs.h1_chr_count$Tissue, dds.fs.h1_chr_count$Chromosome)

dds.fs.h2_chr = na.omit(cbind(dds.fs.h2[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.fs.h2_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fs_h2_inht_chr.csv")

dds.fs.h2_chr_count = as.data.frame(table(dds.fs.h2_chr$inheritance.fs, dds.fs.h2_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.fs.h2_chr_count$Sample = "H2"
dds.fs.h2_chr_count$Tissue = "S"
dds.fs.h2_chr_count$Sex = "Female"
dds.fs.h2_chr_count$Name = paste(dds.fs.h2_chr_count$Sample, dds.fs.h2_chr_count$Sex, dds.fs.h2_chr_count$Tissue, dds.fs.h2_chr_count$Chromosome)

###MG

dds.mg.h2_chr = na.omit(cbind(dds.mg.h2[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.mg.h2_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/mg_h2_inht_chr.csv")

dds.mg.h2_chr_count = as.data.frame(table(dds.mg.h2_chr$inheritance.mg, dds.mg.h2_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.mg.h2_chr_count$Sample = "H2"
dds.mg.h2_chr_count$Tissue = "G"
dds.mg.h2_chr_count$Sex = "Male"
dds.mg.h2_chr_count$Name = paste(dds.mg.h2_chr_count$Sample, dds.mg.h2_chr_count$Sex, dds.mg.h2_chr_count$Tissue, dds.mg.h2_chr_count$Chromosome)


###MS
dds.ms.h1_chr = na.omit(cbind(dds.ms.h1[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.ms.h1_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h1_inht_chr.csv")

dds.ms.h1_chr_count = as.data.frame(table(dds.ms.h1_chr$inheritance.ms, dds.ms.h1_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.ms.h1_chr_count$Sample = "H1"
dds.ms.h1_chr_count$Tissue = "S"
dds.ms.h1_chr_count$Sex = "Male"
dds.ms.h1_chr_count$Name = paste(dds.ms.h1_chr_count$Sample, dds.ms.h1_chr_count$Sex, dds.ms.h1_chr_count$Tissue, dds.ms.h1_chr_count$Chromosome)


dds.ms.h2_chr = na.omit(cbind(dds.ms.h2[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.ms.h2_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h2_inht_chr.csv")

dds.ms.h2_chr_count = as.data.frame(table(dds.ms.h2_chr$inheritance.ms, dds.ms.h2_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.ms.h2_chr_count$Sample = "H2"
dds.ms.h2_chr_count$Tissue = "S"
dds.ms.h2_chr_count$Sex = "Male"
dds.ms.h2_chr_count$Name = paste(dds.ms.h2_chr_count$Sample, dds.ms.h2_chr_count$Sex, dds.ms.h2_chr_count$Tissue, dds.ms.h2_chr_count$Chromosome)


##WM


dds.wm.h1_chr = na.omit(cbind(dds.wm.h1[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.wm.h1_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h1_inht_chr.csv")

dds.wm.h1_chr_count = as.data.frame(table(dds.wm.h1_chr$inheritance.wm, dds.wm.h1_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.wm.h1_chr_count$Sample = "H1"
dds.wm.h1_chr_count$Tissue = "W"
dds.wm.h1_chr_count$Sex = "Male"
dds.wm.h1_chr_count$Name = paste(dds.wm.h1_chr_count$Sample, dds.wm.h1_chr_count$Sex, dds.wm.h1_chr_count$Tissue, dds.wm.h1_chr_count$Chromosome)


dds.wm.h2_chr = na.omit(cbind(dds.wm.h2[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(dds.wm.h2_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h2_inht_chr.csv")


dds.wm.h2_chr_count = as.data.frame(table(dds.wm.h2_chr$inheritance.wm, dds.wm.h2_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

dds.wm.h2_chr_count$Sample = "H2"
dds.wm.h2_chr_count$Tissue = "W"
dds.wm.h2_chr_count$Sex = "Male"
dds.wm.h2_chr_count$Name = paste(dds.wm.h2_chr_count$Sample, dds.wm.h2_chr_count$Sex, dds.wm.h2_chr_count$Tissue, dds.wm.h2_chr_count$Chromosome)

cumulative_chromosomal_inht = rbind(dds.fg.h1_chr_count, dds.fg.h2_chr_count, dds.fs.h1_chr_count, dds.fs.h2_chr_count,
                                    dds.ms.h1_chr_count, dds.ms.h2_chr_count, dds.mg.h2_chr_count, dds.wm.h1_chr_count, dds.wm.h2_chr_count)

cumulative_chromosomal_inht_a = cumulative_chromosomal_inht %>% filter(cumulative_chromosomal_inht$Inheritance != "ambiguous")

##making plot
library(ggstats)
ggplot(cumulative_chromosomal_inht_a) + 
  aes(cumulative_chromosomal_inht_a$Chromosome, fill = cumulative_chromosomal_inht_a$Inheritance, weight = cumulative_chromosomal_inht_a$Number, by = as.factor(Name)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Sample+Sex+Tissue)

######################################################################################################################################
####Autosome vs X enrichment of inheritance categories

##FUNCTION TO CALCULATE FISHER'S EXACT TEST FOR ENRICHMENT (FROM SANTIAGO)
# enrichment
enrichment <- function(x, odds.ratio=FALSE){
  pvals = oddratio = matrix(ncol=ncol(x), nrow=nrow(x))
  if (nrow(x) == 2 & ncol(x) > 2){
    for (i in 1:ncol(x)){
      x1 = cbind(x[,i],rowSums(x[,-i]))
      fet = fisher.test(x1)
      pvals[1,i] = fet$p.value
      oddratio[1,i] = fet$estimate
      fet = fisher.test(x1[c(2,1),])
      pvals[2,i] = fet$p.value
      oddratio[2,i] = fet$estimate
    }
  }
  if (nrow(x) > 2 & ncol(x) == 2){
    for (i in 1:nrow(x)){
      x1 = cbind(x[i,],colSums(x[-i,]))
      fet = fisher.test(x1)
      pvals[i,1] = fet$p.value
      oddratio[i,1] = fet$estimate
      fet = fisher.test(x1[c(2,1),])
      pvals[i,2] = fet$p.value
      oddratio[i,2] = fet$estimate
    }
  }
  if (nrow(x) > 2 & ncol(x) > 2){
    for (i in 1:nrow(x)){
      for (j in 1:ncol(x)){
        x1 = rbind(x[i,], colSums(x[-i,]))
        x1 = cbind(x1[,j], rowSums(x1[,-j]))
        fet = fisher.test(x1)
        pvals[i,j] = fet$p.value
        oddratio[i,j] = fet$estimate
      }
    }
  }
  if (odds.ratio)
    return(oddratio)
  else
    return(pvals)
}
##Comparing autosoems and X-chromosome across FG h1 & H2

#FG H1
dds.fg.h1_chr$A_X = "Autosomes"

dds.fg.h1_chr[(dds.fg.h1_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fg_a_x = dds.fg.h1_chr %>% group_by(A_X, inheritance.fg) %>% dplyr::count()
fg_a_x = fg_a_x %>% dplyr::rename(inheritance = inheritance.fg)
fg_a_x$enrichement = as.numeric(enrichment(matrix(fg_a_x$n, nrow = 7), odds.ratio = T))
fg_a_x$pval = as.numeric(enrichment(matrix(fg_a_x$n, nrow = 7), odds.ratio = F))
#fg_a_x$pval = log10(fg_a_x$pval)

fg_a_x$sig = ""
fg_a_x$sig[ fg_a_x$pval < 0.05 & abs(log2(fg_a_x$enrichement)) > 0.5 ] = "*"

fg_a_x$Sample = "FG_H1"

# FG H2
dds.fg.h2_chr$A_X = "Autosomes"


dds.fg.h2_chr[(dds.fg.h2_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fg_a_x2 = dds.fg.h2_chr %>% group_by(A_X, inheritance.fg) %>% dplyr::count()
fg_a_x2 = fg_a_x2 %>% dplyr::rename(inheritance = inheritance.fg)

fg_a_x2$enrichement = as.numeric(enrichment(matrix(fg_a_x2$n, nrow = 7), odds.ratio = T))
fg_a_x2$pval = as.numeric(enrichment(matrix(fg_a_x2$n, nrow = 7), odds.ratio = F))
#fg_a_x2$pval = log10(fg_a_x2$pval)
fg_a_x2$sig = ""
fg_a_x2$sig[ fg_a_x2$pval < 0.05 & abs(log2(fg_a_x2$enrichement)) > 0.5 ] = "*"

fg_a_x2$Sample = "FG_H2"


#FS H1
dds.fs.h1_chr
dds.fs.h1_chr$A_X = "Autosomes"

dds.fs.h1_chr[(dds.fs.h1_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fs_a_x = dds.fs.h1_chr %>% group_by(A_X, inheritance.fs) %>% dplyr::count()
fs_a_x = fs_a_x %>% dplyr::rename(inheritance = inheritance.fs)
fs_a_x$enrichement = as.numeric(enrichment(matrix(fs_a_x$n, nrow = 7), odds.ratio = T))
fs_a_x$pval = as.numeric(enrichment(matrix(fs_a_x$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)
fs_a_x$sig = ""
fs_a_x$sig[ fs_a_x$pval < 0.05 & abs(log2(fs_a_x$enrichement)) > 0.5 ] = "*"

fs_a_x$Sample = "FS_H1"


##FS H2
dds.fs.h2_chr

dds.fs.h2_chr$A_X = "Autosomes"

dds.fs.h2_chr[(dds.fs.h2_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fs_a_x2 = dds.fs.h2_chr %>% group_by(A_X, inheritance.fs) %>% dplyr::count()
fs_a_x2 = fs_a_x2 %>% dplyr::rename(inheritance = inheritance.fs)
fs_a_x2$enrichement = as.numeric(enrichment(matrix(fs_a_x2$n, nrow = 7), odds.ratio = T))
fs_a_x2$pval = as.numeric(enrichment(matrix(fs_a_x2$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)
fs_a_x2$sig = ""
fs_a_x2$sig[ fs_a_x2$pval < 0.05 & abs(log2(fs_a_x2$enrichement)) > 0.5 ] = "*"

fs_a_x2$Sample = "FS_h2"




####MALES

dds.ms.h1_chr

dds.ms.h1_chr$A_X = "Autosomes"

dds.ms.h1_chr[(dds.ms.h1_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
dds.ms.h1_chr
ms_a_x = dds.ms.h1_chr %>% group_by(A_X, inheritance.ms) %>% dplyr::count()
ms_a_x = ms_a_x %>% dplyr::rename(inheritance = inheritance.ms)
ms_a_x$enrichement = as.numeric(enrichment(matrix(ms_a_x$n, nrow = 7), odds.ratio = T))
ms_a_x$pval = as.numeric(enrichment(matrix(ms_a_x$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)
ms_a_x$sig = ""
ms_a_x$sig[ ms_a_x$pval < 0.05 & abs(log2(ms_a_x$enrichement)) > 0.5 ] = "*"

ms_a_x$Sample = "MS_H1"


dds.ms.h2_chr

dds.ms.h2_chr$A_X = "Autosomes"

dds.ms.h2_chr[(dds.ms.h2_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
dds.ms.h2_chr
ms_a_x2 = dds.ms.h2_chr %>% group_by(A_X, inheritance.ms) %>% dplyr::count()
ms_a_x2 = ms_a_x2 %>% dplyr::rename(inheritance = inheritance.ms)
ms_a_x2$enrichement = as.numeric(enrichment(matrix(ms_a_x2$n, nrow = 7), odds.ratio = T))
ms_a_x2$pval = as.numeric(enrichment(matrix(ms_a_x2$n, nrow = 7), odds.ratio = F))

ms_a_x2$sig = ""
ms_a_x2$sig[ ms_a_x2$pval < 0.05 & abs(log2(ms_a_x2$enrichement)) > 0.5 ] = "*"

ms_a_x2$Sample = "MS_h2"

###MG H2
dds.mg.h2_chr
dds.mg.h2_chr$A_X = "Autosomes"

dds.mg.h2_chr[(dds.mg.h2_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
dds.mg.h2_chr
mg_a_x2 = dds.mg.h2_chr %>% group_by(A_X, inheritance.mg) %>% dplyr::count()
mg_a_x2 = mg_a_x2 %>% dplyr::rename(inheritance = inheritance.mg)
mg_a_x2$enrichement = as.numeric(enrichment(matrix(mg_a_x2$n, nrow = 7), odds.ratio = T))
mg_a_x2$pval = as.numeric(enrichment(matrix(mg_a_x2$n, nrow = 7), odds.ratio = F))

mg_a_x2$sig = ""
mg_a_x2$sig[ mg_a_x2$pval < 0.05 & abs(log2(mg_a_x2$enrichement)) > 0.5 ] = "*"

mg_a_x2$Sample = "MG_h2"
 

###WM H1
dds.wm.h1_chr
dds.wm.h1_chr$A_X = "Autosomes"

dds.wm.h1_chr[(dds.wm.h1_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
dds.wm.h1_chr
wm_a_x1 = dds.wm.h1_chr %>% group_by(A_X, inheritance.wm) %>% dplyr::count()
wm_a_x1 = wm_a_x1 %>% dplyr::rename(inheritance = inheritance.wm)
wm_a_x1$enrichement = as.numeric(enrichment(matrix(wm_a_x1$n, nrow = 7), odds.ratio = T))
wm_a_x1$pval = as.numeric(enrichment(matrix(wm_a_x1$n, nrow = 7), odds.ratio = F))

wm_a_x1$sig = ""
wm_a_x1$sig[ wm_a_x1$pval < 0.05 & abs(log2(wm_a_x1$enrichement)) > 0.5 ] = "*"

wm_a_x1$Sample = "WM_h1"

##WM H2
dds.wm.h2_chr
dds.wm.h2_chr$A_X = "Autosomes"

dds.wm.h2_chr[(dds.wm.h2_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
dds.wm.h2_chr
wm_a_x2 = dds.wm.h2_chr %>% group_by(A_X, inheritance.wm) %>% dplyr::count()
wm_a_x2 = wm_a_x2 %>% dplyr::rename(inheritance = inheritance.wm)
wm_a_x2$enrichement = as.numeric(enrichment(matrix(wm_a_x2$n, nrow = 7), odds.ratio = T))
wm_a_x2$pval = as.numeric(enrichment(matrix(wm_a_x2$n, nrow = 7), odds.ratio = F))

wm_a_x2$sig = ""
wm_a_x2$sig[ wm_a_x2$pval < 0.05 & abs(log2(wm_a_x2$enrichement)) > 0.5 ] = "*"

wm_a_x2$Sample = "WM_h2"


##joining all the inheritance enrichment data

inheritance_a_x = rbind(fg_a_x, fg_a_x2, fs_a_x, fs_a_x2, ms_a_x, ms_a_x2, mg_a_x2, wm_a_x1, wm_a_x2)


ggplot(inheritance_a_x %>% filter(A_X == "X"), aes(x=log2(enrichement), y=inheritance)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=inheritance), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="autosomes", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="X", y=7.5, x=.75, fill="black", color="white") +
  facet_wrap(~Sample)+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_AxX_enrichment.pdf", dpi = 300, width = 12, height = 12 )






######################################################################################################################################
######################################################################################################################################
 
#regulatory divergence

##################################################################################3
##FG
fg_h1_ase_chr = rownames_to_column(fg_h1_ase)
rownames(fg_h1_ase_chr) = fg_h1_ase_chr$rowname

fg_h1_ase_chr = na.omit(cbind(fg_h1_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                                orthologs_chr$chromosome))
write.csv(fg_h1_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fg_h1_ase_chr.csv")


fg_h1_ase_chr_count = as.data.frame(table(fg_h1_ase_chr$type, fg_h1_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fg_h1_ase_chr_count$Sample = "h1"
fg_h1_ase_chr_count$Tissue = "G"
fg_h1_ase_chr_count$Sex = "Female"
fg_h1_ase_chr_count$Name = paste(fg_h1_ase_chr_count$Sample, fg_h1_ase_chr_count$Sex, fg_h1_ase_chr_count$Tissue, fg_h1_ase_chr_count$Chromosome)

##FG H2
fg_h2_ase_chr = rownames_to_column(fg_h2_ase)
rownames(fg_h2_ase_chr) = fg_h2_ase_chr$rowname

fg_h2_ase_chr = na.omit(cbind(fg_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(fg_h2_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fg_h2_ase_chr.csv")

fg_h2_ase_chr_count = as.data.frame(table(fg_h2_ase_chr$type, fg_h2_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fg_h2_ase_chr_count$Sample = "h2"
fg_h2_ase_chr_count$Tissue = "G"
fg_h2_ase_chr_count$Sex = "Female"
fg_h2_ase_chr_count$Name = paste(fg_h2_ase_chr_count$Sample, fg_h2_ase_chr_count$Sex, fg_h2_ase_chr_count$Tissue, fg_h2_ase_chr_count$Chromosome)


##FS

fs_h1_ase_chr = rownames_to_column(fs_h1_ase)
rownames(fs_h1_ase_chr) = fs_h1_ase_chr$rowname

fs_h1_ase_chr = na.omit(cbind(fs_h1_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))
write.csv(fs_h1_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fs_h1_ase_chr.csv")


fs_h1_ase_chr_count = as.data.frame(table(fs_h1_ase_chr$type, fs_h1_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fs_h1_ase_chr_count$Sample = "h1"
fs_h1_ase_chr_count$Tissue = "S"
fs_h1_ase_chr_count$Sex = "Female"
fs_h1_ase_chr_count$Name = paste(fs_h1_ase_chr_count$Sample, fs_h1_ase_chr_count$Sex, fs_h1_ase_chr_count$Tissue, fs_h1_ase_chr_count$Chromosome)

#FS H2

fs_h2_ase_chr = rownames_to_column(fs_h2_ase)
rownames(fs_h2_ase_chr) = fs_h2_ase_chr$rowname

fs_h2_ase_chr = na.omit(cbind(fs_h2_ase_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                              orthologs_chr$chromosome))

write.csv(fs_h2_ase_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/fs_h2_ase_chr.csv")


fs_h2_ase_chr_count = as.data.frame(table(fs_h2_ase_chr$type, fs_h2_ase_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

fs_h2_ase_chr_count$Sample = "h2"
fs_h2_ase_chr_count$Tissue = "S"
fs_h2_ase_chr_count$Sex = "Female"
fs_h2_ase_chr_count$Name = paste(fs_h2_ase_chr_count$Sample, fs_h2_ase_chr_count$Sex, fs_h2_ase_chr_count$Tissue, fs_h2_ase_chr_count$Chromosome)


##MS

ms_h1_ase_3_chr = ms_h1_ase_3
rownames(ms_h1_ase_3_chr) = ms_h1_ase_3$rowname

ms_h1_ase_3_chr = na.omit(cbind(ms_h1_ase_3_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                                orthologs_chr$chromosome))

write.csv(ms_h1_ase_3_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h1_ase_chr.csv")


ms_h1_ase_3_chr_count = as.data.frame(table(ms_h1_ase_3_chr$reg.div, ms_h1_ase_3_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)
ms_h1_ase_3_chr_count$Sample = "h1"
ms_h1_ase_3_chr_count$Tissue = "S"
ms_h1_ase_3_chr_count$Sex = "male"
ms_h1_ase_3_chr_count$Name = paste(ms_h1_ase_3_chr_count$Sample, ms_h1_ase_3_chr_count$Sex, ms_h1_ase_3_chr_count$Tissue, ms_h1_ase_3_chr_count$Chromosome)

#MS H2

ms_h2_ase_3_chr = ms_h2_ase_3
rownames(ms_h2_ase_3_chr) = ms_h2_ase_3$rowname

ms_h2_ase_3_chr = na.omit(cbind(ms_h2_ase_3_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                                orthologs_chr$chromosome))
write.csv(ms_h2_ase_3_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/ms_h2_ase_chr.csv")

ms_h2_ase_3_chr_count = as.data.frame(table(ms_h2_ase_3_chr$reg.div, ms_h2_ase_3_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)
ms_h2_ase_3_chr_count$Sample = "h2"
ms_h2_ase_3_chr_count$Tissue = "S"
ms_h2_ase_3_chr_count$Sex = "male"
ms_h2_ase_3_chr_count$Name = paste(ms_h2_ase_3_chr_count$Sample, ms_h2_ase_3_chr_count$Sex, ms_h2_ase_3_chr_count$Tissue, ms_h2_ase_3_chr_count$Chromosome)


##MG
mg_h2_ase_3_chr = mg_h2_ase_3
rownames(mg_h2_ase_3_chr) = mg_h2_ase_3$rowname

mg_h2_ase_3_chr = na.omit(cbind(mg_h2_ase_3_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                                orthologs_chr$chromosome))
write.csv(mg_h2_ase_3_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/mg_h2_ase_chr.csv")

mg_h2_ase_3_chr_count = as.data.frame(table(mg_h2_ase_3_chr$reg.div, mg_h2_ase_3_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

mg_h2_ase_3_chr_count$Sample = "h2"
mg_h2_ase_3_chr_count$Tissue = "G"
mg_h2_ase_3_chr_count$Sex = "male"
mg_h2_ase_3_chr_count$Name = paste(mg_h2_ase_3_chr_count$Sample, mg_h2_ase_3_chr_count$Sex, mg_h2_ase_3_chr_count$Tissue,mg_h2_ase_3_chr_count$Chromosome)


#WM

wm_h1_ase_3_chr = wm_h1_ase_3
rownames(wm_h1_ase_3_chr) = wm_h1_ase_3$rowname

wm_h1_ase_3_chr = na.omit(cbind(wm_h1_ase_3_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                                orthologs_chr$chromosome))

write.csv(wm_h1_ase_3_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h1_ase_chr.csv")


wm_h1_ase_3_chr_count = as.data.frame(table(wm_h1_ase_3_chr$reg.div, wm_h1_ase_3_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

wm_h1_ase_3_chr_count$Sample = "h1"
wm_h1_ase_3_chr_count$Tissue = "W"
wm_h1_ase_3_chr_count$Sex = "Male"
wm_h1_ase_3_chr_count$Name = paste(wm_h1_ase_3_chr_count$Sample, wm_h1_ase_3_chr_count$Sex, wm_h1_ase_3_chr_count$Tissue, wm_h1_ase_3_chr_count$Chromosome)


##H2
wm_h2_ase_3_chr = wm_h2_ase_3
rownames(wm_h2_ase_3_chr) = wm_h2_ase_3$rowname

wm_h2_ase_3_chr = na.omit(cbind(wm_h2_ase_3_chr[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ], orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name,
                                orthologs_chr$chromosome))

write.csv(wm_h2_ase_3_chr, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/wm_h2_ase_chr.csv")

wm_h2_ase_3_chr_count = as.data.frame(table(wm_h2_ase_3_chr$reg.div, wm_h2_ase_3_chr$`orthologs_chr$chromosome`)) %>%
  rename(Inheritance = Var1, Chromosome = Var2, Number = Freq)

wm_h2_ase_3_chr_count$Sample = "h2"
wm_h2_ase_3_chr_count$Tissue = "W"
wm_h2_ase_3_chr_count$Sex = "Male"
wm_h2_ase_3_chr_count$Name = paste(wm_h2_ase_3_chr_count$Sample, wm_h2_ase_3_chr_count$Sex, wm_h2_ase_3_chr_count$Tissue, wm_h2_ase_3_chr_count$Chromosome)


cumulative_chromosomal_regdiv = rbind(fg_h1_ase_chr_count, fg_h2_ase_chr_count, fs_h1_ase_chr_count, fs_h2_ase_chr_count,
                                      ms_h1_ase_3_chr_count, ms_h2_ase_3_chr_count, mg_h2_ase_3_chr_count, wm_h1_ase_3_chr_count, wm_h2_ase_3_chr_count)

cumulative_chromosomal_regdiv_a = cumulative_chromosomal_regdiv %>% filter(cumulative_chromosomal_regdiv$Inheritance != "ambiguous")

##making plot
library(ggstats)
ggplot(cumulative_chromosomal_regdiv_a) + 
  aes(cumulative_chromosomal_regdiv_a$Chromosome, fill = cumulative_chromosomal_regdiv_a$Inheritance, weight = cumulative_chromosomal_regdiv_a$Number, by = as.factor(Name)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Sample+Sex+Tissue)


######################################################################################################################################
####
##Comparing autosomes and X-chromosome across FG h1 & H2

#FG H1
fg_h1_ase_chr$A_X = "Autosomes"

fg_h1_ase_chr[(fg_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fg_a_x_ase = fg_h1_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

fg_a_x_ase$enrichement = as.numeric(enrichment(matrix(fg_a_x_ase$n, nrow = 7), odds.ratio = T))
fg_a_x_ase$pval = as.numeric(enrichment(matrix(fg_a_x_ase$n, nrow = 7), odds.ratio = F))
#fg_a_x$pval = log10(fg_a_x$pval)

fg_a_x_ase$sig = ""
fg_a_x_ase$sig[ fg_a_x_ase$pval < 0.05 & abs(log2(fg_a_x_ase$enrichement)) > 0.5 ] = "*"

fg_a_x_ase$Sample = "FG_H1"


#FG H2
fg_h2_ase_chr$A_X = "Autosomes"

fg_h2_ase_chr[(fg_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fg_a_x2_ase = fg_h2_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

fg_a_x2_ase$enrichement = as.numeric(enrichment(matrix(fg_a_x2_ase$n, nrow = 7), odds.ratio = T))
fg_a_x2_ase$pval = as.numeric(enrichment(matrix(fg_a_x2_ase$n, nrow = 7), odds.ratio = F))
#fg_a_x$pval = log10(fg_a_x$pval)

fg_a_x2_ase$sig = ""
fg_a_x2_ase$sig[ fg_a_x2_ase$pval < 0.05 & abs(log2(fg_a_x2_ase$enrichement)) > 0.5 ] = "*"

fg_a_x2_ase$Sample = "FG_h2"



#fs H1
fs_h1_ase_chr$A_X = "Autosomes"

fs_h1_ase_chr[(fs_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fs_a_x_ase = fs_h1_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

fs_a_x_ase$enrichement = as.numeric(enrichment(matrix(fs_a_x_ase$n, nrow = 7), odds.ratio = T))
fs_a_x_ase$pval = as.numeric(enrichment(matrix(fs_a_x_ase$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

fs_a_x_ase$sig = ""
fs_a_x_ase$sig[ fs_a_x_ase$pval < 0.05 & abs(log2(fs_a_x_ase$enrichement)) > 0.5 ] = "*"

fs_a_x_ase$Sample = "fs_H1"


#fs H2
fs_h2_ase_chr$A_X = "Autosomes"

fs_h2_ase_chr[(fs_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
fs_a_x2_ase = fs_h2_ase_chr %>% group_by(A_X, type) %>% dplyr::count()

fs_a_x2_ase$enrichement = as.numeric(enrichment(matrix(fs_a_x2_ase$n, nrow = 7), odds.ratio = T))
fs_a_x2_ase$pval = as.numeric(enrichment(matrix(fs_a_x2_ase$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

fs_a_x2_ase$sig = ""
fs_a_x2_ase$sig[ fs_a_x2_ase$pval < 0.05 & abs(log2(fs_a_x2_ase$enrichement)) > 0.5 ] = "*"

fs_a_x2_ase$Sample = "FS_h2"

#joining all the inheritance enrichment data

female_regdiv_a_x = rbind(fg_a_x_ase, fg_a_x2_ase, fs_a_x_ase, fs_a_x2_ase)


ggplot(female_regdiv_a_x %>% filter(A_X == "X"), aes(x=log2(enrichement), y=type)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=type), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="autosomes", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="X", y=7.5, x=.75, fill="black", color="white") +
  facet_wrap(~Sample)+
  scale_x_continuous(breaks=c(-3, -2, -1, 0, 3, 2, 1),labels=c(-3, -2, -1, 0, 3, 2, 1))+
  # xlim(-4, 4)+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_femalesAxX_enrichment.pdf", dpi = 300, width = 12, height = 12 )


##MALES
#IN MALES, THE x-LINKED GENE WERE CATEGORISED SLIGHTLY DIFFERENTLY, HENCE WE NEED TO CHANGE AUTOSOMAL GENES TO THAT SO ITS A DIRECT COMPARISON

##ms h1
ms_h1_ase_3_chr$A_X = "Autosomes"
ms_h1_ase_3_chr[(ms_h1_ase_3_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

ms_h1_ase_3_chr$inht_changed = ms_h1_ase_3_chr$reg.div
ms_h1_ase_3_chr$inht_changed[ ms_h1_ase_3_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
ms_h1_ase_3_chr$inht_changed[ ms_h1_ase_3_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
ms_h1_ase_3_chr$inht_changed[ ms_h1_ase_3_chr$inht_changed == "ambiguous" ] = "Other"
table(ms_h1_ase_3_chr$inht_changed)

ms_a_x_ase = ms_h1_ase_3_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

ms_a_x_ase$enrichement = as.numeric(enrichment(matrix(ms_a_x_ase$n, nrow = 5), odds.ratio = T))
ms_a_x_ase$pval = as.numeric(enrichment(matrix(ms_a_x_ase$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

ms_a_x_ase$sig = ""
ms_a_x_ase$sig[ ms_a_x_ase$pval < 0.05 & abs(log2(ms_a_x_ase$enrichement)) > 0.5 ] = "*"

ms_a_x_ase$Sample = "MS_H1"


##MS H2
ms_h2_ase_3_chr$A_X = "Autosomes"
ms_h2_ase_3_chr[(ms_h2_ase_3_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

ms_h2_ase_3_chr$inht_changed = ms_h2_ase_3_chr$reg.div
ms_h2_ase_3_chr$inht_changed[ ms_h2_ase_3_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
ms_h2_ase_3_chr$inht_changed[ ms_h2_ase_3_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
ms_h2_ase_3_chr$inht_changed[ ms_h2_ase_3_chr$inht_changed == "ambiguous" ] = "Other"
table(ms_h2_ase_3_chr$inht_changed)

ms_a_x2_ase = ms_h2_ase_3_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

ms_a_x2_ase$enrichement = as.numeric(enrichment(matrix(ms_a_x2_ase$n, nrow = 5), odds.ratio = T))
ms_a_x2_ase$pval = as.numeric(enrichment(matrix(ms_a_x2_ase$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

ms_a_x2_ase$sig = ""
ms_a_x2_ase$sig[ ms_a_x2_ase$pval < 0.05 & abs(log2(ms_a_x2_ase$enrichement)) > 0.5 ] = "*"

ms_a_x2_ase$Sample = "MS_H2"


##MG H2
mg_h2_ase_3_chr$A_X = "Autosomes"
mg_h2_ase_3_chr[(mg_h2_ase_3_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

mg_h2_ase_3_chr$inht_changed = mg_h2_ase_3_chr$reg.div
mg_h2_ase_3_chr$inht_changed[ mg_h2_ase_3_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
mg_h2_ase_3_chr$inht_changed[ mg_h2_ase_3_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
mg_h2_ase_3_chr$inht_changed[ mg_h2_ase_3_chr$inht_changed == "ambiguous" ] = "Other"
table(mg_h2_ase_3_chr$inht_changed)

mg_a_x2_ase = mg_h2_ase_3_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

mg_a_x2_ase$enrichement = as.numeric(enrichment(matrix(mg_a_x2_ase$n, nrow = 5), odds.ratio = T))
mg_a_x2_ase$pval = as.numeric(enrichment(matrix(mg_a_x2_ase$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

mg_a_x2_ase$sig = ""
mg_a_x2_ase$sig[ mg_a_x2_ase$pval < 0.05 & abs(log2(mg_a_x2_ase$enrichement)) > 0.5 ] = "*"

mg_a_x2_ase$Sample = "MG_H2"


##WM H1

wm_h1_ase_3_chr$A_X = "Autosomes"
wm_h1_ase_3_chr[(wm_h1_ase_3_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

wm_h1_ase_3_chr$inht_changed = wm_h1_ase_3_chr$reg.div
wm_h1_ase_3_chr$inht_changed[ wm_h1_ase_3_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
wm_h1_ase_3_chr$inht_changed[ wm_h1_ase_3_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
wm_h1_ase_3_chr$inht_changed[ wm_h1_ase_3_chr$inht_changed == "ambiguous" ] = "Other"
table(wm_h1_ase_3_chr$inht_changed)

wm_a_x_ase = wm_h1_ase_3_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

wm_a_x_ase$enrichement = as.numeric(enrichment(matrix(wm_a_x_ase$n, nrow = 5), odds.ratio = T))
wm_a_x_ase$pval = as.numeric(enrichment(matrix(wm_a_x_ase$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

wm_a_x_ase$sig = ""
wm_a_x_ase$sig[ wm_a_x_ase$pval < 0.05 & abs(log2(wm_a_x_ase$enrichement)) > 0.5 ] = "*"

wm_a_x_ase$Sample = "WM_H1"

##WM H2

wm_h2_ase_3_chr$A_X = "Autosomes"
wm_h2_ase_3_chr[(wm_h2_ase_3_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

wm_h2_ase_3_chr$inht_changed = wm_h2_ase_3_chr$reg.div
wm_h2_ase_3_chr$inht_changed[ wm_h2_ase_3_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
wm_h2_ase_3_chr$inht_changed[ wm_h2_ase_3_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
wm_h2_ase_3_chr$inht_changed[ wm_h2_ase_3_chr$inht_changed == "ambiguous" ] = "Other"
table(wm_h2_ase_3_chr$inht_changed)

wm_a_x2_ase = wm_h2_ase_3_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

wm_a_x2_ase$enrichement = as.numeric(enrichment(matrix(wm_a_x2_ase$n, nrow = 5), odds.ratio = T))
wm_a_x2_ase$pval = as.numeric(enrichment(matrix(wm_a_x2_ase$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

wm_a_x2_ase$sig = ""
wm_a_x2_ase$sig[ wm_a_x2_ase$pval < 0.05 & abs(log2(wm_a_x2_ase$enrichement)) > 0.5 ] = "*"

wm_a_x2_ase$Sample = "WM_H2"


#joining all the inheritance enrichment data

male_regdiv_a_x = rbind(mg_a_x2_ase, ms_a_x_ase, ms_a_x2_ase, wm_a_x_ase, wm_a_x2_ase)


ggplot(male_regdiv_a_x %>% filter(A_X == "X"), aes(x=log2(enrichement), y=inht_changed)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=inht_changed), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="autosomes", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="X", y=7.5, x=.75, fill="black", color="white") +
  facet_wrap(~Sample)+
  xlim(-4, 4)+
  scale_x_continuous(breaks=c(-3, -2, -1, 0, 3, 2, 1),labels=c(-3, -2, -1, 0, 3, 2, 1))+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))

#Positive values of the log2 odds ratios for females and males denote enrichment for the X-chromosome 
#and negative values indicate enrichment for autosomes.

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_malesAxX_enrichment.pdf", dpi = 300, width = 12, height = 12 )


######################################################################################################################
######COMPARING MALES AND FEMALES BY JOINING CATEGORIES TO MAKE TOEHR CAETGORY IN FEMALES 

####FIG 7#########################################            
##fs h1
fs_h1_ase_chr$A_X = "Autosomes"
fs_h1_ase_chr[(fs_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

fs_h1_ase_chr$inht_changed = fs_h1_ase_chr$type
fs_h1_ase_chr$inht_changed[ fs_h1_ase_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
fs_h1_ase_chr$inht_changed[ fs_h1_ase_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
fs_h1_ase_chr$inht_changed[ fs_h1_ase_chr$inht_changed == "ambiguous" ] = "Other"
table(fs_h1_ase_chr$inht_changed)

fs_a_x_ase_1 = fs_h1_ase_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

fs_a_x_ase_1$enrichement = as.numeric(enrichment(matrix(fs_a_x_ase_1$n, nrow = 5), odds.ratio = T))
fs_a_x_ase_1$pval = as.numeric(enrichment(matrix(fs_a_x_ase_1$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

fs_a_x_ase_1$sig = ""
fs_a_x_ase_1$sig[ fs_a_x_ase_1$pval < 0.05 & abs(log2(fs_a_x_ase_1$enrichement)) > 0.5 ] = "*"

fs_a_x_ase_1$Sample = "fs_H1"


##fs H2
fs_h2_ase_chr$A_X = "Autosomes"
fs_h2_ase_chr[(fs_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

fs_h2_ase_chr$inht_changed = fs_h2_ase_chr$type
fs_h2_ase_chr$inht_changed[ fs_h2_ase_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
fs_h2_ase_chr$inht_changed[ fs_h2_ase_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
fs_h2_ase_chr$inht_changed[ fs_h2_ase_chr$inht_changed == "ambiguous" ] = "Other"
table(fs_h2_ase_chr$inht_changed)

fs_a_x2_ase_2 = fs_h2_ase_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

fs_a_x2_ase_2$enrichement = as.numeric(enrichment(matrix(fs_a_x2_ase_2$n, nrow = 5), odds.ratio = T))
fs_a_x2_ase_2$pval = as.numeric(enrichment(matrix(fs_a_x2_ase_2$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

fs_a_x2_ase_2$sig = ""
fs_a_x2_ase_2$sig[ fs_a_x2_ase_2$pval < 0.05 & abs(log2(fs_a_x2_ase_2$enrichement)) > 0.5 ] = "*"

fs_a_x2_ase_2$Sample = "fs_H2"


##fg h1
fg_h1_ase_chr$A_X = "Autosomes"
fg_h1_ase_chr[(fg_h1_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

fg_h1_ase_chr$inht_changed = fg_h1_ase_chr$type
fg_h1_ase_chr$inht_changed[ fg_h1_ase_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
fg_h1_ase_chr$inht_changed[ fg_h1_ase_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
fg_h1_ase_chr$inht_changed[ fg_h1_ase_chr$inht_changed == "ambiguous" ] = "Other"
table(fg_h1_ase_chr$inht_changed)

fg_a_x_ase_1 = fg_h1_ase_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

fg_a_x_ase_1$enrichement = as.numeric(enrichment(matrix(fg_a_x_ase_1$n, nrow = 5), odds.ratio = T))
fg_a_x_ase_1$pval = as.numeric(enrichment(matrix(fg_a_x_ase_1$n, nrow = 5), odds.ratio = F))
#fg_a_x$pval = log10(fg_a_x$pval)

fg_a_x_ase_1$sig = ""
fg_a_x_ase_1$sig[ fg_a_x_ase_1$pval < 0.05 & abs(log2(fg_a_x_ase_1$enrichement)) > 0.5 ] = "*"

fg_a_x_ase_1$Sample = "fg_H1"


##fg H2
fg_h2_ase_chr$A_X = "Autosomes"
fg_h2_ase_chr[(fg_h2_ase_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"

fg_h2_ase_chr$inht_changed = fg_h2_ase_chr$type
fg_h2_ase_chr$inht_changed[ fg_h2_ase_chr$inht_changed == "cis x trans (compensatory)" ] = "Other"
fg_h2_ase_chr$inht_changed[ fg_h2_ase_chr$inht_changed == "cis + trans (enhancing)" ] = "Other"
fg_h2_ase_chr$inht_changed[ fg_h2_ase_chr$inht_changed == "ambiguous" ] = "Other"
table(fg_h2_ase_chr$inht_changed)

fg_a_x2_ase_2 = fg_h2_ase_chr %>% group_by(A_X, inht_changed) %>% dplyr::count()

fg_a_x2_ase_2$enrichement = as.numeric(enrichment(matrix(fg_a_x2_ase_2$n, nrow = 5), odds.ratio = T))
fg_a_x2_ase_2$pval = as.numeric(enrichment(matrix(fg_a_x2_ase_2$n, nrow = 5), odds.ratio = F))
#fg_a_x$pval = log10(fg_a_x$pval)

fg_a_x2_ase_2$sig = ""
fg_a_x2_ase_2$sig[ fg_a_x2_ase_2$pval < 0.05 & abs(log2(fg_a_x2_ase_2$enrichement)) > 0.5 ] = "*"

fg_a_x2_ase_2$Sample = "fg_H2"


malevsfemale_ax = rbind(mg_a_x2_ase, ms_a_x_ase, ms_a_x2_ase, wm_a_x_ase, wm_a_x2_ase, 
                        fg_a_x2_ase_2, fg_a_x_ase_1, fs_a_x2_ase_2, fs_a_x_ase_1)

ggplot(malevsfemale_ax %>% filter(A_X == "X"), aes(x=log2(enrichement), y=inht_changed)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=inht_changed), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="autosomes", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="X", y=7.5, x=.75, fill="black", color="white") +
  facet_wrap(~Sample)+
  xlim(-4, 4)+
  scale_x_continuous(breaks=c(-3, -2, -1, 0, 3, 2, 1),labels=c(-3, -2, -1, 0, 3, 2, 1))+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_malesvsfemalesAxX_enrichment.pdf", dpi = 300, width = 12, height = 12 )

################################################################################
##Comparing X-linked expression between H1 and H2 in MS AND WM

##MS 
ms_h1_ase_3_chr_X = ms_h1_ase_3_chr %>% filter(A_X == "X") #filtering data for X-linked genes
ms_h1_ase_3_chr_X$Sample = "MS_H1"

ms_h2_ase_3_chr_X = ms_h2_ase_3_chr %>% filter(A_X == "X") #filtering data for X-linked genes
ms_h2_ase_3_chr_X$Sample = "MS_H2"

ms_h1h2_X = rbind(ms_h1_ase_3_chr_X, ms_h2_ase_3_chr_X)
ms_h1h2_X = ms_h1h2_X %>% group_by(Sample, reg.div) %>% dplyr::count()

ms_h1h2_X$enrichement = as.numeric(enrichment(matrix(ms_h1h2_X$n, nrow = 5), odds.ratio = T))
ms_h1h2_X$pval = as.numeric(enrichment(matrix(ms_h1h2_X$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

ms_h1h2_X$sig = ""
ms_h1h2_X$sig[ ms_h1h2_X$pval < 0.05 & abs(log2(ms_h1h2_X$enrichement)) > 0.5 ] = "*"

ggplot(ms_h1h2_X %>% filter(Sample == "MS_H1"), aes(x=log2(enrichement), y=reg.div)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=reg.div), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="H2", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="H1", y=7.5, x=.75, fill="black", color="white") +
  #facet_wrap(~Sample)+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))

#Positive values of the log2 odds ratios for females and males denote enrichment for the H1 
#and negative values indicate enrichment for H2.

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_ms_h1h2_X_enrichment.pdf", dpi = 300, width = 12, height = 12 )


##WM 
wm_h1_ase_3_chr_X = wm_h1_ase_3_chr %>% filter(A_X == "X") #filtering data for X-linked genes
wm_h1_ase_3_chr_X$Sample = "WM_H1"

wm_h2_ase_3_chr_X = wm_h2_ase_3_chr %>% filter(A_X == "X") #filtering data for X-linked genes
wm_h2_ase_3_chr_X$Sample = "WM_H2"

wm_h1h2_X = rbind(wm_h1_ase_3_chr_X, wm_h2_ase_3_chr_X)
wm_h1h2_X = wm_h1h2_X %>% group_by(Sample, reg.div) %>% dplyr::count()

wm_h1h2_X$enrichement = as.numeric(enrichment(matrix(wm_h1h2_X$n, nrow = 5), odds.ratio = T))
wm_h1h2_X$pval = as.numeric(enrichment(matrix(wm_h1h2_X$n, nrow = 5), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

wm_h1h2_X$sig = ""
wm_h1h2_X$sig[ wm_h1h2_X$pval < 0.05 & abs(log2(wm_h1h2_X$enrichement)) > 0.5 ] = "*"

ggplot(wm_h1h2_X %>% filter(Sample == "WM_H1"), aes(x=log2(enrichement), y=reg.div)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=reg.div), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="H2", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="H1", y=7.5, x=.75, fill="black", color="white") +
  #facet_wrap(~Sample)+
  scale_x_continuous(breaks=c(-0.8, -0.4, 0, 0.4, 0.8),labels=c(-0.8, -0.4, 0, 0.4, 0.8))+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))

#Positive values of the log2 odds ratios for females and males denote enrichment for the H1 
#and negative values indicate enrichment for H2.

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_wm_h1h2_X_enrichment.pdf", dpi = 300, width = 12, height = 12 )


################################################################################
##Comparing Autosomal-linked expression between H1 and H2 in MS AND WM

##MS 
ms_h1_ase_3_chr_A = ms_h1_ase_3_chr %>% filter(A_X == "Autosomes") #filtering data for X-linked genes
ms_h1_ase_3_chr_A$Sample = "MS_H1"

ms_h2_ase_3_chr_A = ms_h2_ase_3_chr %>% filter(A_X == "Autosomes") #filtering data for X-linked genes
ms_h2_ase_3_chr_A$Sample = "MS_H2"

ms_h1h2_A = rbind(ms_h1_ase_3_chr_A, ms_h2_ase_3_chr_A)
ms_h1h2_A = ms_h1h2_A %>% group_by(Sample, reg.div) %>% dplyr::count()

ms_h1h2_A$enrichement = as.numeric(enrichment(matrix(ms_h1h2_A$n, nrow = 7), odds.ratio = T))
ms_h1h2_A$pval = as.numeric(enrichment(matrix(ms_h1h2_A$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

ms_h1h2_A$sig = ""
ms_h1h2_A$sig[ ms_h1h2_A$pval < 0.05 & abs(log2(ms_h1h2_A$enrichement)) > 0.5 ] = "*"

ggplot(ms_h1h2_A %>% filter(Sample == "MS_H1"), aes(x=log2(enrichement), y=reg.div)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=reg.div), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="H2", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="H1", y=7.5, x=.75, fill="black", color="white") +
  #facet_wrap(~Sample)+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))

#Positive values of the log2 odds ratios for females and males denote enrichment for the H1 
#and negative values indicate enrichment for H2.

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_ms_h1h2_Autosomal_enrichment.pdf", dpi = 300, width = 12, height = 12 )


##WM 
wm_h1_ase_3_chr_A = wm_h1_ase_3_chr %>% filter(A_X == "Autosomes") #filtering data for X-linked genes
wm_h1_ase_3_chr_A$Sample = "WM_H1"

wm_h2_ase_3_chr_A = wm_h2_ase_3_chr %>% filter(A_X == "Autosomes") #filtering data for X-linked genes
wm_h2_ase_3_chr_A$Sample = "WM_H2"

wm_h1h2_A = rbind(wm_h1_ase_3_chr_A, wm_h2_ase_3_chr_A)
wm_h1h2_A = wm_h1h2_A %>% group_by(Sample, reg.div) %>% dplyr::count()

wm_h1h2_A$enrichement = as.numeric(enrichment(matrix(wm_h1h2_A$n, nrow = 7), odds.ratio = T))
wm_h1h2_A$pval = as.numeric(enrichment(matrix(wm_h1h2_A$n, nrow = 7), odds.ratio = F))
#fs_a_x$pval = log10(fs_a_x$pval)

wm_h1h2_A$sig = ""
wm_h1h2_A$sig[ wm_h1h2_A$pval < 0.05 & abs(log2(wm_h1h2_A$enrichement)) > 0.5 ] = "*"

ggplot(wm_h1h2_A %>% filter(Sample == "WM_H1"), aes(x=log2(enrichement), y=reg.div)) +
  #geom_rect(data=rectX, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
  geom_vline(xintercept=0, size=1) +
  geom_bar(aes(fill=reg.div), stat="identity", position=position_dodge(width=0.9), alpha=0.7, show.legend=F)+
  annotate(geom="label", label="H2", y=7.5, x=-.75, fill="black", color="white") +
  annotate(geom="label", label="H1", y=7.5, x=.75, fill="black", color="white") +
  #facet_wrap(~Sample)+
  labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  geom_text(aes(x=log2(enrichement)/2, label=sig), size=5, show.legend=F, color="white") +
  theme(legend.text=element_text(hjust=0),
        strip.background = element_rect(fill="grey90", color=NA))

#Positive values of the log2 odds ratios for females and males denote enrichment for the H1 
#and negative values indicate enrichment for H2.

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_wm_h1h2_Autosomal_enrichment.pdf", dpi = 300, width = 12, height = 12 )




