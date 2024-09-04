library(dplyr)


####VISUALIZING PROPORTION OF GENES IN EACH SAMPLE FOR REGULATORY DIVERGENCE AND INHERITANCE CATEGORY

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

regdiv_all_counts = rbind(fg_h1_table.prop, fg_h2_table.prop, fs_h1_table.prop, fs_h2_table.prop,
                          ms_h1_table.prop, ms_h2_table.prop, mg_h2_table.prop, wm_h1_table.prop, wm_h2_table.prop)

#write.csv(regdiv_all_counts, file = "Allele Specific Expression/Hybrid Analysis/Basic Results/RegDiv_all_counts.csv")

#######Figures looking at regulatory divergence patterns

#############
##Fig3A
###Barplot of regualtory divergence frequencies across samples and hybrids

##plotting proportions across chromosomes
library(ggstats)
library(dplyr)
library(ggplot2)

regdiv_all_counts = read.csv("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/RegDiv_all_counts.csv", sep = ",")
View(regdiv_all_counts)
###subsetting to get rid of ambisuous counts
regdiv_all_counts$sample = paste(regdiv_all_counts$Sex, regdiv_all_counts$Tissue)
regdiv_all_counts_a = regdiv_all_counts %>% filter(regdiv_all_counts$Var1 != "ambiguous" & regdiv_all_counts$Var1 != "Other") %>% select(-Proportion)
View(regdiv_all_counts_a)

###Plotting proportion of genes in each category in each sample in both hybrids without amiguous category
ggplot(regdiv_all_counts_a) +
  aes(x = regdiv_all_counts_a$sample, fill = regdiv_all_counts_a$Var1, weight = regdiv_all_counts_a$Freq, by = as.factor(regdiv_all_counts_a$sample)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Regulatory Divergence Pattern")

#ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/fig3_regdivpattern.pdf", dpi = 300, width = 18, height = 8 )

###Plotting proportion of genes in each category in each sample in both hybrids with amiguous category
ggplot(regdiv_all_counts) +
  aes(x = regdiv_all_counts$sample, fill = regdiv_all_counts$Var1, weight = regdiv_all_counts$Freq, by = as.factor(regdiv_all_counts$sample)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Regulatory Divergence Pattern")

ggsave(filename = "Figures_and_Results/Part 3/fig3_2_regdivpattern_ambiguous.pdf", dpi = 300, width = 18, height = 8 )


#plotting regulatory divergence plots for only WM samples
regdiv_all_counts_WM = regdiv_all_counts_a %>% filter(regdiv_all_counts_a$Tissue == "W")
View(regdiv_all_counts_WM)
ggplot(regdiv_all_counts_WM) +
  aes(x = regdiv_all_counts_WM$sample, fill = regdiv_all_counts_WM$Var1, weight = regdiv_all_counts_WM$Freq, by = as.factor(regdiv_all_counts_WM$sample)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Regulatory Divergence Pattern")

ggplot(regdiv_all_counts_WM) +
  aes(x = regdiv_all_counts_WM$Var1, y = regdiv_all_counts_WM$Freq, fill = regdiv_all_counts_WM$Var1) +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Regulatory Divergence Pattern")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_barplot_WM.pdf", dpi = 300, width = 12, height = 8 )


##########################################################################################################3
######FIG 3B########################################################################################

####VISUALIZING INHERITANCE BETWEEN MALES AND FEMALES OF H2

###Proportion of genes between males and females 
regdiv_mgxfg_h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/RegDiv_mgxfg_h2.txt", sep="\t", head=T, comment.char="#")
View(regdiv_mgxfg_h2)

ggplot(regdiv_mgxfg_h2, aes(x=regdiv_mgxfg_h2$Female_percent, y=regdiv_mgxfg_h2$Male_percentage, color = regdiv_mgxfg_h2$Reg.Div.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 17) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Male Gonad")))) +
  xlab(expression(paste(bold("Female Gonad")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H2 MGXFG Regulatory Divergence")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/rd_fgxmg_h2.pdf", dpi = 300, width = 10, height = 8 )



##########################################################################################################3
######FIG 3C########################################################################################

###VISUALIZING INHERITANCE BETWEEN MS X FS FOR H1 and H2

regdiv_fsxms_h1_h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_msxfs_h1h2.txt", sep="\t", head=T, comment.char="#")
View(regdiv_fsxms_h1_h2)

ggplot(regdiv_fsxms_h1_h2, aes(x=regdiv_fsxms_h1_h2$Female_percent, y=regdiv_fsxms_h1_h2$Male_percentage, color = regdiv_fsxms_h1_h2$Reg.Div.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Male Soma")))) +
  xlab(expression(paste(bold("Female Soma")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H2xH1 MSXFS Regulatory Divergence")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/rd_fsxms_h1h2.pdf", dpi = 300, width = 10, height = 8 )

###VISUALIZING INHERITANCE ACROSS MALE GONAD AND MALE SOMA IN H2

regdiv_mgxms_h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/RegDiv_mgxms_h2.txt", sep="\t", head=T, comment.char="#")
View(regdiv_mgxms_h2)


ggplot(regdiv_mgxms_h2, aes(x=regdiv_mgxms_h2$Gonad_percentage, y=regdiv_mgxms_h2$Soma_percent, color = regdiv_mgxms_h2$Reg.Div.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 17) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Male Soma")))) +
  xlab(expression(paste(bold("Male Gonad")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H2 MGXMS Regulatory Divergence")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/rd_mgxms_h2.pdf", dpi = 300, width = 10, height = 8 )

###VISUALIZING INHERITANCE IN FG X FS IN H1 AND H2

reg_div_fgxfs_h1h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/regdiv_fgxfs_h1h2.txt", sep="\t", head=T, comment.char="#")

View(reg_div_fgxfs_h1h2)
ggplot(reg_div_fgxfs_h1h2, aes(x=reg_div_fgxfs_h1h2$Gonad.Percentage, y=reg_div_fgxfs_h1h2$Soma.Percentage, color = reg_div_fgxfs_h1h2$Reg.Div.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Female Soma")))) +
  xlab(expression(paste(bold("Female Gonad")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H1XH2 FGXFS Regulatory Divergence")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/rd_fgxfs_h1h2.pdf", dpi = 300, width = 10, height = 8 )

##################################################################################################################

#####FIG 3C
##identifying if sex-biased genes overlap with reg div genes 
###For example - maybe male-biased genes are upregulated in females and show up as compensatory

####1. Run 1.WT_general DGEresult to get sex-biased genes 

####Importing file with sex-biased genes
sex_biased_genes = read.csv(file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 1/mxf_res_genes.csv", header = TRUE,
                            comment.char = "#")
rownames(sex_biased_genes) = sex_biased_genes$X
View(sex_biased_genes)

###################################################################################################
####Regulatory Divergence for FG,FS,MS,MG, WM from 3.Allele Specific expression.R 
###################################################################################################

##1. Run 2.Classifying INheritance and get X-linked regulatory divergence
##2. Run AlleleSpecfic Expression.R and get regulatory divergence for Autosomal and X-linked genes
##3. Now run this entire script , it is similar to 2.2.Inheritance Plots.R

#########
###H1_FG#
#########
head(fg_h1_ase)
nrow(fg_h1_ase)

#filtering larger sex-biased genes to get fg genes
sex_biased_genes_fg_h1_ase = sex_biased_genes[fg_h1_ase$genes, ]
head(sex_biased_genes_fg_h1_ase)

sex_biased_genes_fg_h1_ase = na.omit(sex_biased_genes_fg_h1_ase)
dds.fg.h1_2_ase = fg_h1_ase[sex_biased_genes_fg_h1_ase$X, ]
ase_fg_h1_sb = cbind(sex_biased_genes_fg_h1_ase, dds.fg.h1_2_ase) %>% select("X", "mxf", "genes", "type")
ase_fg_h1_sb_table = as.data.frame(table(ase_fg_h1_sb$mxf, ase_fg_h1_sb$type))
ase_fg_h1_sb_table$Sample = "FG_H1"



#########
###H2_FG#
#########
head(fg_h2_ase)
nrow(fg_h2_ase)

#filtering larger sex-biased genes to get fg genes
sex_biased_genes_fg_h2_ase = sex_biased_genes[fg_h2_ase$genes, ]
head(sex_biased_genes_fg_h2_ase)

sex_biased_genes_fg_h2_ase = na.omit(sex_biased_genes_fg_h2_ase)
dds.fg.h2_2_ase = fg_h2_ase[sex_biased_genes_fg_h2_ase$X, ]
ase_fg_h2_sb = cbind(sex_biased_genes_fg_h2_ase, dds.fg.h2_2_ase) %>% select("X", "mxf", "genes", "type")
ase_fg_h2_sb_table = as.data.frame(table(ase_fg_h2_sb$mxf, ase_fg_h2_sb$type))
ase_fg_h2_sb_table$Sample = "FG_H2"

#########
###H1_FS#
#########
head(fs_h1_ase)


#filtering larger sex-biased genes to get fs genes
sex_biased_genes_fs_h1_ase = sex_biased_genes[fs_h1_ase$genes, ]
head(sex_biased_genes_fs_h1_ase)

sex_biased_genes_fs_h1_ase = na.omit(sex_biased_genes_fs_h1_ase)
dds.fs.h1_2_ase = fs_h1_ase[sex_biased_genes_fs_h1_ase$X, ]
ase_fs_h1_sb = cbind(sex_biased_genes_fs_h1_ase, dds.fs.h1_2_ase) %>% select("X", "mxf", "genes", "type")
ase_fs_h1_sb_table = as.data.frame(table(ase_fs_h1_sb$mxf, ase_fs_h1_sb$type))
ase_fs_h1_sb_table$Sample = "FS_H1"


#########
###H2_FS#
#########
head(fs_h2_ase)


#filtering larger sex-biased genes to get fs genes
sex_biased_genes_fs_h2_ase = sex_biased_genes[fs_h2_ase$genes, ]
head(sex_biased_genes_fs_h2_ase)

sex_biased_genes_fs_h2_ase = na.omit(sex_biased_genes_fs_h2_ase)
dds.fs.h2_2_ase = fs_h2_ase[sex_biased_genes_fs_h2_ase$X, ]
ase_fs_h2_sb = cbind(sex_biased_genes_fs_h2_ase, dds.fs.h2_2_ase) %>% select("X", "mxf", "genes", "type")
ase_fs_h2_sb_table = as.data.frame(table(ase_fs_h2_sb$mxf, ase_fs_h2_sb$type))
ase_fs_h2_sb_table$Sample = "FS_H2"


#########
###H1_MS#
#########
head(ms_h1_ase_3)

row.names(ms_h1_ase_3) = ms_h1_ase_3$rowname
#filtering larger sex-biased genes to get fs genes
sex_biased_genes_ms_h1_ase = sex_biased_genes[ms_h1_ase_3$rowname, ]
head(sex_biased_genes_ms_h1_ase)

sex_biased_genes_ms_h1_ase = na.omit(sex_biased_genes_ms_h1_ase)
dds.ms.h1_2_ase = ms_h1_ase_3[sex_biased_genes_ms_h1_ase$X, ]
ase_ms_h1_sb = cbind(sex_biased_genes_ms_h1_ase, dds.ms.h1_2_ase) %>% select("X", "mxf", "rowname", "reg.div")
ase_ms_h1_sb_table = as.data.frame(table(ase_ms_h1_sb$mxf, ase_ms_h1_sb$reg.div))
ase_ms_h1_sb_table$Sample = "MS_H1"

#########
###H2_MS#
#########
head(ms_h2_ase_3)
row.names(ms_h2_ase_3) = ms_h2_ase_3$rowname

#filtering larger sex-biased genes to get fs genes
sex_biased_genes_ms_h2_ase = sex_biased_genes[ms_h2_ase_3$rowname, ]
head(sex_biased_genes_ms_h2_ase)

sex_biased_genes_ms_h2_ase = na.omit(sex_biased_genes_ms_h2_ase)
dds.ms.h2_2_ase = ms_h2_ase_3[sex_biased_genes_ms_h2_ase$X, ]
ase_ms_h2_sb = cbind(sex_biased_genes_ms_h2_ase, dds.ms.h2_2_ase) %>% select("X", "mxf", "rowname", "reg.div")
ase_ms_h2_sb_table = as.data.frame(table(ase_ms_h2_sb$mxf, ase_ms_h2_sb$reg.div))
ase_ms_h2_sb_table$Sample = "MS_H2"


#########
###H2_MG#
#########
###MAKE SURE TO TAKE THE ASE DATA AFTER YOU REPLACED DIVERGENCE FOR X-LINKED GENES
head(mg_h2_ase_3)
row.names(mg_h2_ase_3) = mg_h2_ase_3$rowname

#filtering larger sex-biased genes to get fs genes
sex_biased_genes_mg_h2_ase = sex_biased_genes[mg_h2_ase_3$rowname, ]
head(sex_biased_genes_mg_h2_ase)
nrow(sex_biased_genes_mg_h2_ase)

sex_biased_genes_mg_h2_ase = na.omit(sex_biased_genes_mg_h2_ase)
dds.mg.h2_2_ase = mg_h2_ase_3[sex_biased_genes_mg_h2_ase$X, ]
head(mg_h2_ase_3)


ase_mg_h2_sb = cbind(sex_biased_genes_mg_h2_ase, dds.mg.h2_2_ase) %>% select("X", "mxf", "rowname", "reg.div")
ase_mg_h2_sb_table = as.data.frame(table(ase_mg_h2_sb$mxf, ase_mg_h2_sb$reg.div))
ase_mg_h2_sb_table$Sample = "MG_H2"

####Merging all frequency tables
head(ase_ms_h2_sb_table)

sex_bised_regdiv = rbind(ase_fg_h2_sb_table, ase_fg_h1_sb_table, ase_fs_h2_sb_table, ase_fs_h1_sb_table,
                         ase_mg_h2_sb_table, ase_ms_h1_sb_table, ase_ms_h2_sb_table)
View(sex_bised_regdiv)

write.csv(sex_bised_regdiv, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/sex_biasedandregdiv.csv")


sex_bised_regdiv_a = sex_bised_regdiv %>% filter(Var2 != "ambiguous" & Var2 != "Other" & Var2 != "conserved")

###Heatmap for sex-biased genes and regualtroy divergence

ggplot(sex_bised_regdiv_a %>% filter(Var2 != "conserved"), aes(x= Var2, y=Var1, fill=Freq)) +
  geom_tile() +
  xlab("") + ylab("") +
  #scale_fill_distiller(palette = "RdPu", direction=-1)
  #scale_fill_gradient(low="#ebebff", high="#0d2b95")+ 
  scale_fill_gradient(low="#e8ecf0", high="#00315c")+
  #viridis::scale_fill_viridis(option = "B", direction=-1, name="number of \ngenes")+
  facet_rep_grid(~Sample) +
  geom_text(aes(label=Freq)) +
  #scale_y_discrete(labels=type.legend[3:7], limits=type.levels[3:7]) +
  #scale_x_discrete(labels=class.legend[3:7][c(4:5,1:3)], limits=class.levels[3:7][c(4:5,1:3)]) +
  theme(strip.background = element_rect(fill="grey90", color=NA),
        strip.text.y = element_text(size=10),
        legend.title=element_text(size=10),
        legend.text =element_text(size=10),
        axis.text.y = element_text(size=10),
        #axis.text.x = element_blank())
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))

#ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/sex_biasedvsregdiv_heatmap.pdf", dpi = 300, width = 20, height = 8 )

###Barplots
##Barplot

ggplot(sex_bised_regdiv_a, aes(x = sex_bised_regdiv_a$Var1, y = sex_bised_regdiv_a$Freq, fill = sex_bised_regdiv_a$Var2))+
  geom_bar(stat = "identity")+
  facet_rep_wrap(~Sample)

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/sex_bsvsregdiv_barplot1.pdf", dpi = 300, width = 20, height = 8 )

ggplot(sex_bised_regdiv_a, aes(x = sex_bised_regdiv_a$Var2, y = sex_bised_regdiv_a$Freq, fill = sex_bised_regdiv_a$Var1))+
  geom_bar(stat = "identity")+
  facet_rep_wrap(~Sample)

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/sex_bsvsregdiv_barplot2.pdf", dpi = 300, width = 20, height = 8 )

###Line plot 
sex_bised_regdiv = read.csv("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/sex_biasedandregdiv.csv", sep = ",")

sex_bised_regdiv$Var1 = as.factor(sex_bised_regdiv$Var1)
sex_bised_regdiv$Var2 = as.factor(sex_bised_regdiv$Var2)

sex_bised_regdiv %>% filter(Var2 != "ambiguous" & Var2 != "Other") %>%
  ggplot(aes(x = factor(Var2, level = c("cis-trans (compensatory)", "cis x trans (compensatory)", "cis + trans (enhancing)", "cis-only", "trans-only", "conserved")), Freq, col = Var1, group = Var1)) +
  #geom_hline(yintercept=0, color = "black", linetype="dashed")+
  #scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  geom_text(aes(label = Freq), nudge_x = 0.25, nudge_y = 0.5) +
  geom_line() +
  theme_bw() +
  facet_wrap(~Sample) +
  #stat_smooth(method = "lm") +
  ggtitle("Sex biased vs regdiv")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/sex_bsvsregdiv_lineplot.pdf", dpi = 300, width = 20, height = 8 )

levels(sex_bised_regdiv$Var2)


#########################################################################################################3
#######VISUALIZING REGULATORY DIVERGENCE BETWEEN H1 AND H2
#########################################################################################################3

regdiv_h1xh2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/RegDiv_all_h1xh2.txt", sep="\t", head=T, comment.char="#")
View(regdiv_h1xh2)
ggplot(regdiv_h1xh2, aes(x=regdiv_h1xh2$H1_percentage, y=regdiv_h1xh2$H2_percent, color = regdiv_h1xh2$Reg.Div.Category, shape = Tissue)) +
  geom_point(alpha=0.8, show.legend=T, size=8) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("H2")))) +
  xlab(expression(paste(bold("H1")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  facet_rep_wrap(~Sex) +
  theme(strip.background = element_rect(fill="grey90", color=NA))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/regdiv_h1xh2_1.pdf", dpi = 300, width = 15, height = 8 )

ggplot(regdiv_h1xh2, aes(x=regdiv_h1xh2$H1_percentage, y=regdiv_h1xh2$H2_percent, color = regdiv_h1xh2$Reg.Div.Category, shape = Sex)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 15) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("H2")))) +
  xlab(expression(paste(bold("H1")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  facet_rep_wrap(~Tissue+Sex) +
  theme(strip.background = element_rect(fill="grey90", color=NA))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/regdiv_h1xh2_2.pdf", dpi = 300, width = 18, height = 8 )


###
##SPECIFCALLY FOR MALE SOMA, I HAD TO SPLIT THE PERCENTAGE PLOTS (H1 X H2) INTO AUTOSOMAL AND X-LINKED GENES
#INCLUDED AUTOSOMAL DATA IN THE MAIN FIRGURE FOR MALE SOMA

###creating the file with A vs X info
###check chromosomal distribution.R file
# ms_h1_ase_3_chr$A_X = "Autosomes"
# ms_h1_ase_3_chr[(ms_h1_ase_3_chr$`orthologs_chr$chromosome` == "X"), "A_X"] = "X"
# ms_h1_ax_count = as.data.frame(table(ms_h1_ase_3_chr$reg.div, ms_h1_ase_3_chr$A_X)) 
# ms_h1_ax_count$Hybrid = "H1"
# ms_h1_ax_count = ms_h1_ax_count %>% dplyr::rename(regdiv_1 = Var1, Chr_1 = Var2, H1_count = Freq)
# 
# ms_h2_ax_count = as.data.frame(table(ms_h2_ase_3_chr$reg.div, ms_h2_ase_3_chr$A_X))
# ms_h2_ax_count$Hybrid = "H2"
# ms_h2_ax_count = ms_h2_ax_count %>% dplyr::rename(regdiv_2 = Var1, Chr_2 = Var2, H2_count = Freq)
# 
# ms_h1h2_ax_counts = cbind(ms_h1_ax_count, ms_h2_ax_count)
# 
# write.csv(ms_h1h2_ax_counts, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/RegDiv_MS_A_X_counts.csv")

ms_axx_count = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/regdiv_MS_H1H1_AxX_counts.txt", sep="\t", head=T, comment.char="#")
#View(ms_axx_count)

ggplot(ms_axx_count, aes(x=ms_axx_count$H1_percentage, y=ms_axx_count$H2_percentage, color = ms_axx_count$regdiv_1)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 15) +
  scale_color_brewer(palette = "Dark2") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("H2")))) +
  xlab(expression(paste(bold("H1")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  facet_rep_wrap(~Chr_1) +
  theme(strip.background = element_rect(fill="grey90", color=NA))

#ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/regdiv_MSh1xh2_2.pdf", dpi = 300, width = 18, height = 8 )

#######################################################################################################
###COMAPRING H1 AND H2 
#######################################################################################################

###REGULATORY DIVERGENCE BETWEEN H1 AND H2

###########
###FG######
###########

####NEED TO RUN 2.CLASSIFYING INHERITANCE.R FIRST
####THEN RUN 3.ALLELE SPECIFCI EXPRESSION.R
##GET REGULATORY DIVERGENC FOR AUTOSOMAL AND X-LINKED GENES
###VARIABLES SAVED FROM 3.ALLELE SPECIFCI EXPRESSION.R

head(fg_h1_ase)

`%notin%` <- Negate(`%in%`)

##check to see if dataframes have gene names 
head(fg_h1_ase)
head(fg_h2_ase)

##If they already have rownames as column names
dds.fg.h1.1_ase = fg_h1_ase
dds.fg.h2.1_ase = fg_h2_ase

##ordering each dataframe adn filtering to make sure both datafrmae have the same number of rows
dds.fg.h1.1_ase = dds.fg.h1.1_ase[order(row.names(dds.fg.h1.1_ase)), ]
dds.fg.h1.1_ase = subset(dds.fg.h1.1_ase, dds.fg.h1.1_ase$genes %in% dds.fg.h2.1_ase$genes)

dim(dds.fg.h1.1_ase)
head(dds.fg.h2.1_ase)


dds.fg.h2.1_ase = dds.fg.h2.1_ase[order(row.names(dds.fg.h2.1_ase)), ]
dds.fg.h2.1_ase = dplyr::rename(dds.fg.h2.1_ase, logFC.sp_2 = logFC.sp, logFC.ase.h2 = logFC.ase, DE.h2 = DE, p.value.sp.h2 = p.value.sp,
         p.value.ase.h2 = p.value.ase, genes.h2 = genes, trans_effect.h2 = trans_effect, type.h2 = type, 
         species.h2 = species, sex.h2 = sex)

names(dds.fg.h2.1_ase)
nrow(dds.fg.h1.1_ase)
nrow(dds.fg.h2.1_ase)

###joining both dataframes

dds.fg.h1h2_ase = cbind(dds.fg.h1.1_ase, dds.fg.h2.1_ase) 
nrow(dds.fg.h1h2_ase)

#####Visualizing using scatterplot
ggplot(dds.fg.h1h2_ase %>% filter(type != "ambiguous" & type.h2 != "ambiguous"), aes(x = type, type.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

#####VISUALIZING USING BUBBLE PLOTS

head(as.data.frame(table(dds.fg.h1h2_ase$type, dds.fg.h1h2_ase$type.h2)))

##creating intersection in gene categories
ase.fg.h1h2.intersection = as.data.frame(table(dds.fg.h1h2_ase$type, dds.fg.h1h2_ase$type.h2)) %>%
  dplyr::rename(H1 = Var1, H2=Var2, Intersection = Freq)


ase.fg.h1h2.intersection = ase.fg.h1h2.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ase.fg.h1h2.intersection$Sample = "FG"

ggplot(data = ase.fg.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/fg.h1h2.regdiv.pdf", dpi = 300, width = 18, height = 8 )



###########
###FS######
###########

####NEED TO RUN 2.CLASSIFYING INHERITANCE.R FIRST
####THEN RUN 3.ALLELE SPECIFCI EXPRESSION.R
##GET REGULATORY DIVERGENC FOR AUTOSOMAL AND X-LINKED GENES
###VARIABLES SAVED FROM 3.ALLELE SPECIFCI EXPRESSION.R

head(fs_h1_ase)

`%notin%` <- Negate(`%in%`)

##check to see if dataframes have gene names 
head(fs_h1_ase)
head(fs_h2_ase)

##If they already have rownames as column names
dds.fs.h1.1_ase = fs_h1_ase
dds.fs.h2.1_ase = fs_h2_ase

##ordering each dataframe adn filtering to make sure both datafrmae have the same number of rows
dds.fs.h1.1_ase = dds.fs.h1.1_ase[order(row.names(dds.fs.h1.1_ase)), ]
dds.fs.h1.1_ase = subset(dds.fs.h1.1_ase, dds.fs.h1.1_ase$genes %in% dds.fs.h2.1_ase$genes)

dim(dds.fs.h1.1_ase)
head(dds.fs.h2.1_ase)


dds.fs.h2.1_ase = dds.fs.h2.1_ase[order(row.names(dds.fs.h2.1_ase)), ] %>%
  dplyr::rename(logFC.sp_2 = logFC.sp, logFC.ase.h2 = logFC.ase, DE.h2 = DE, p.value.sp.h2 = p.value.sp,
         p.value.ase.h2 = p.value.ase, genes.h2 = genes, trans_effect.h2 = trans_effect, type.h2 = type, 
         species.h2 = species, sex.h2 = sex)


nrow(dds.fs.h1.1_ase)
nrow(dds.fs.h2.1_ase)

###joining both dataframes

dds.fs.h1h2_ase = cbind(dds.fs.h1.1_ase, dds.fs.h2.1_ase) 
dim(dds.fs.h1h2_ase)

#####Visualizing using scatterplot
ggplot(dds.fs.h1h2_ase %>% filter(type != "ambiguous" & type.h2 != "ambiguous"), aes(x = type, type.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

#####VISUALIZING USING BUBBLE PLOTS

head(as.data.frame(table(dds.fs.h1h2_ase$type, dds.fs.h1h2_ase$type.h2)))

##creating intersection in gene categories
ase.fs.h1h2.intersection = as.data.frame(table(dds.fs.h1h2_ase$type, dds.fs.h1h2_ase$type.h2)) %>%
  dplyr::rename(H1 = Var1, H2=Var2, Intersection = Freq)

sum(ase.fs.h1h2.intersection$Intersection) #7755 genes in bubble plots after removing ambiguous genes
ase.fs.h1h2.intersection = ase.fs.h1h2.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ase.fs.h1h2.intersection$Sample = "FS"

ggplot(data = ase.fs.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/fs.h1h2.regdiv.pdf", dpi = 300, width = 18, height = 8 )




###########
###MS######
###########

####NEED TO RUN 2.CLASSIFYING INHERITANCE.R FIRST
####THEN RUN 3.ALLELE SPECIFCI EXPRESSION.R
##GET REGULATORY DIVERGENC FOR AUTOSOMAL AND X-LINKED GENES
###VARIABLES SAVED FROM 3.ALLELE SPECIFCI EXPRESSION.R
ms_h1_ase_3 = rbind(ms_h1_ase_2, dds.ms.h1_X_1) ### data with updated X-linked regulatory categories
ms_h2_ase_3 = rbind(ms_h2_ase_2, dds.ms.h2_X_1)

head(ms_h1_ase)
head(ms_h1_ase_3)

table(ms_h1_ase_3$reg.div)
table(ms_h1_ase$type)


`%notin%` <- Negate(`%in%`)

##check to see if dataframes have gene names 
head(ms_h1_ase_3)
head(ms_h2_ase_3)
nrow(ms_h1_ase_3)
head(ms_h2_ase_3)

##If they already have rownames as column names
dds.ms.h1.1_ase = ms_h1_ase_3 %>% dplyr::rename(genes = rowname)
dds.ms.h2.1_ase = ms_h2_ase_3 %>% dplyr::rename(genes = rowname)
nrow(dds.ms.h1.1_ase)
##ordering each dataframe adn filtering to make sure both datafrmae have the same number of rows
dds.ms.h1.1_ase = dds.ms.h1.1_ase[order((dds.ms.h1.1_ase$genes)), ]
dds.ms.h1.1_ase = subset(dds.ms.h1.1_ase, dds.ms.h1.1_ase$genes %in% dds.ms.h2.1_ase$genes)

dim(dds.ms.h1.1_ase)
dim(dds.ms.h2.1_ase)

dds.ms.h2.1_ase = dds.ms.h2.1_ase[order((dds.ms.h2.1_ase$genes)), ] %>%
  dplyr::rename(genes.h2 = genes, reg.div.h2 = reg.div, chr_info.h2 = chr_info)

names(dds.ms.h2.1_ase)

head(dds.ms.h2.1_ase)
head(dds.ms.h1.1_ase)
# nrow(dds.ms.h1.1_ase)
# nrow(dds.ms.h2.1_ase)

###joining both dataframes

dds.ms.h1h2_ase = cbind(dds.ms.h1.1_ase, dds.ms.h2.1_ase) 
names(dds.ms.h1h2_ase)
nrow(dds.ms.h1h2_ase)

#####Visualizing using scatterplot
ggplot(dds.ms.h1h2_ase %>% filter(reg.div != "ambiguous" & reg.div.h2 != "ambiguous"), aes(x = reg.div, reg.div.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  facet_wrap(~chr_info)+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

ggplot(dds.ms.h1h2_ase, aes(x = reg.div, reg.div.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  facet_wrap(~chr_info)+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

#####VISUALIZING USING BUBBLE PLOTS

head(as.data.frame(table(dds.ms.h1h2_ase$reg.div, dds.ms.h1h2_ase$reg.div.h2)))

##creating intersection in gene categories
##for autosomal genes
dds.ms.h1h2_ase_A = dds.ms.h1h2_ase %>% filter(chr_info == "A" & chr_info.h2 == "A") ##11120 genes

nrow(dds.ms.h1h2_ase %>% filter(chr_info.h2 == "A"))

nrow(dds.ms.h1h2_ase_A)


ase.ms.h1h2.intersection_A = as.data.frame(table(dds.ms.h1h2_ase_A$reg.div, dds.ms.h1h2_ase_A$reg.div.h2)) %>%
  dplyr::rename(H1 = Var1, H2=Var2, Intersection = Freq)

sum(ase.ms.h1h2.intersection_A$Intersection) #11120 before ambiguous REMOVED, 6515 after removing ambiguous

ase.ms.h1h2.intersection_A = ase.ms.h1h2.intersection_A %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ase.ms.h1h2.intersection_A$Sample = "MS"

ggplot(data = ase.ms.h1h2.intersection_A, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/ms.h1h2_Autosomal.regdiv.pdf", dpi = 300, width = 18, height = 8 )

##For X-linked genes

dds.ms.h1h2_ase_X = dds.ms.h1h2_ase %>% filter(chr_info == "X" & chr_info.h2 == "X")
nrow(dds.ms.h1h2_ase_X)

ase.ms.h1h2.intersection = as.data.frame(table(dds.ms.h1h2_ase_X$reg.div, dds.ms.h1h2_ase_X$reg.div.h2)) %>%
  dplyr::rename(H1 = Var1, H2=Var2, Intersection = Freq)
ase.ms.h1h2.intersection = ase.ms.h1h2.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ase.ms.h1h2.intersection$Sample = "MS"

ggplot(data = ase.ms.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/ms.h1h2_X_linked.regdiv.pdf", dpi = 300, width = 18, height = 8 )


##########################################################################################33
###MAKING CUMULATIVE BUBBLEPLOT TO MAKE ALL UNITS/SCALE FOR LEGEND SAME

##COMMBINING ALL INTERSECTION DATA (except MS X-linked and WM)

head(ase.fg.h1h2.intersection)
head(ase.fs.h1h2.intersection)
head(ase.ms.h1h2.intersection_A)

cumulative_bubbleplot = rbind(ase.fg.h1h2.intersection, ase.fs.h1h2.intersection, ase.ms.h1h2.intersection_A)


ggplot(data = cumulative_bubbleplot, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  facet_wrap(~Sample)+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size_continuous("Number of genes", breaks=c(50, 100, 500,1000,1500,2000,2500,3000, 4000,5000),labels=c(50, 100, 500,1000,1500,2000,2500,3000, 4000,5000),range = c(1,15))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/all.h1h2.regdiv.pdf", dpi = 300, width = 22, height = 8 )


###########
###WM######
###########

####NEED TO RUN 2.CLASSIFYING INHERITANCE.R FIRST
####THEN RUN 3.ALLELE SPECIFCI EXPRESSION.R
##GET REGULATORY DIVERGENC FOR AUTOSOMAL AND X-LINKED GENES
###VARIABLES SAVED FROM 3.ALLELE SPECIFCI EXPRESSION.R
wm_h1_ase_3 = rbind(wm_h1_ase_2, dds.wm.h1_X_1) ### data with updated X-linked regulatory categories
wm_h2_ase_3 = rbind(wm_h2_ase_2, dds.wm.h2_X_1)

head(wm_h1_ase)
head(wm_h1_ase_3)

table(wm_h1_ase_3$reg.div)
table(wm_h1_ase$type)


`%notin%` <- Negate(`%in%`)

##check to see if dataframes have gene names 
head(wm_h1_ase_3)
head(wm_h2_ase_3)
nrow(wm_h1_ase_3)
nrow(wm_h2_ase_3)

##If they already have rownames as column names
dds.wm.h1.1_ase = wm_h1_ase_3 %>% dplyr::rename(genes = rowname)
dds.wm.h2.1_ase = wm_h2_ase_3 %>% dplyr::rename(genes = rowname)
nrow(dds.wm.h1.1_ase)
##ordering each dataframe adn filtering to make sure both datafrmae have the same number of rows
dds.wm.h1.1_ase = dds.wm.h1.1_ase[order((dds.wm.h1.1_ase$genes)), ]
dds.wm.h1.1_ase = subset(dds.wm.h1.1_ase, dds.wm.h1.1_ase$genes %in% dds.wm.h2.1_ase$genes)

dim(dds.wm.h1.1_ase)
dim(dds.wm.h2.1_ase)


dds.wm.h2.1_ase = dds.wm.h2.1_ase[order((dds.wm.h2.1_ase$genes)), ] %>%
  dplyr::rename(genes.h2 = genes, reg.div.h2 = reg.div, chr_info.h2 = chr_info)

names(dds.wm.h2.1_ase)

head(dds.wm.h2.1_ase)
head(dds.wm.h1.1_ase)
# nrow(dds.wm.h1.1_ase)
# nrow(dds.wm.h2.1_ase)

###joining both dataframes

dds.wm.h1h2_ase = cbind(dds.wm.h1.1_ase, dds.wm.h2.1_ase) 
names(dds.wm.h1h2_ase)
nrow(dds.wm.h1h2_ase)

#####Visualizing using scatterplot
ggplot(dds.wm.h1h2_ase %>% filter(reg.div != "ambiguous" & reg.div.h2 != "ambiguous"), aes(x = reg.div, reg.div.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  facet_wrap(~chr_info)+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

ggplot(dds.wm.h1h2_ase, aes(x = reg.div, reg.div.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  facet_wrap(~chr_info)+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

#####VISUALIZING USING BUBBLE PLOTS

head(as.data.frame(table(dds.wm.h1h2_ase$reg.div, dds.wm.h1h2_ase$reg.div.h2)))

##creating intersection in gene categories
##for autosomal genes
dds.wm.h1h2_ase_A = dds.wm.h1h2_ase %>% filter(chr_info == "A" & chr_info.h2 == "A")

nrow(dds.wm.h1h2_ase %>% filter(chr_info.h2 == "A"))

nrow(dds.wm.h1h2_ase_A)


ase.wm.h1h2.intersection = as.data.frame(table(dds.wm.h1h2_ase_A$reg.div, dds.wm.h1h2_ase_A$reg.div.h2)) %>%
  dplyr::rename(H1 = Var1, H2=Var2, Intersection = Freq)
ase.wm.h1h2.intersection = ase.wm.h1h2.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ase.wm.h1h2.intersection$Sample = "wm"

ggplot(data = ase.wm.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/wm.h1h2_Autosomal.regdiv.pdf", dpi = 300, width = 18, height = 8 )

##For X-linked genes

dds.wm.h1h2_ase_X = dds.wm.h1h2_ase %>% filter(chr_info == "X" & chr_info.h2 == "X")
nrow(dds.wm.h1h2_ase_X)

ase.wm.h1h2.intersection = as.data.frame(table(dds.wm.h1h2_ase_X$reg.div, dds.wm.h1h2_ase_X$reg.div.h2)) %>%
  dplyr::rename(H1 = Var1, H2=Var2, Intersection = Freq)
ase.wm.h1h2.intersection = ase.wm.h1h2.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ase.wm.h1h2.intersection$Sample = "wm"

ggplot(data = ase.wm.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/wm.h1h2_X_linked.regdiv.pdf", dpi = 300, width = 18, height = 8 )



##########################################################################################################
#####Looking at regulatory divergence underlying transgressive genes
############################################################################################3

######
##H1_FG##
######

###VARIABLES SAVED FROM 2.Classifying Inheritance.R

#inheritance
View(dds.fg.h1)
#allele specific expression
View(fg_h1_ase)

##getting only gene names and inheritance for each gene
dds.fg.h1_a = dds.fg.h1
fg_h1_ase_a = fg_h1_ase[dds.fg.h1_a$rowname, ] ##filtering genes in inheritance file
fg_h1_rd_inht = cbind(dds.fg.h1_a, fg_h1_ase_a) %>% select("inheritance.fg", "type", "genes")
fg_h1_rd_inht_frq = as.data.frame(table(fg_h1_rd_inht$inheritance.fg, fg_h1_rd_inht$type))
fg_h1_rd_inht_frq$Sample = "FG_H1"

######
##H2_FG##
######

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R

#inheritance
View(dds.fg.h2)
#allele specific expression
View(fg_h2_ase)

##getting only gene names and inheritance for each gene
dds.fg.h2_a = dds.fg.h2
fg_h2_ase_a = fg_h2_ase[dds.fg.h2_a$rowname, ] ##filtering genes in inheritance file
fg_h2_rd_inht = cbind(dds.fg.h2_a, fg_h2_ase_a) %>% select("inheritance.fg", "type", "genes")
fg_h2_rd_inht_frq = as.data.frame(table(fg_h2_rd_inht$inheritance.fg, fg_h2_rd_inht$type))
fg_h2_rd_inht_frq$Sample = "FG_H2"


######
##H1_FS##
######

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R

#inheritance
View(dds.fs.h1)
#allele specific expression
View(fs_h1_ase)

##getting only gene names and inheritance for each gene
dds.fs.h1_a = dds.fs.h1
fs_h1_ase_a = fs_h1_ase[dds.fs.h1_a$rowname, ] ##filtering genes in inheritance file
fs_h1_rd_inht = cbind(dds.fs.h1_a, fs_h1_ase_a) %>% select("inheritance.fs", "type", "genes")
fs_h1_rd_inht_frq = as.data.frame(table(fs_h1_rd_inht$inheritance.fs, fs_h1_rd_inht$type))
fs_h1_rd_inht_frq$Sample = "FS_H1"

######
##H2_FS##
######

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R

#inheritance
View(dds.fs.h2)
#allele specific expression
View(fs_h2_ase)

##getting only gene names and inheritance for each gene
dds.fs.h2_a = dds.fs.h2
fs_h2_ase_a = fs_h2_ase[dds.fs.h2_a$rowname, ] ##filtering genes in inheritance file
fs_h2_rd_inht = cbind(dds.fs.h2_a, fs_h2_ase_a) %>% select("inheritance.fs", "type", "genes")
fs_h2_rd_inht_frq = as.data.frame(table(fs_h2_rd_inht$inheritance.fs, fs_h2_rd_inht$type))
fs_h2_rd_inht_frq$Sample = "FS_H2"


######
##H1_MS##
######

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R

#inheritance
View(dds.ms.h1)
#allele specific expression
View(ms_h1_ase)

##getting only gene names and inheritance for each gene
dds.ms.h1_a = dds.ms.h1
ms_h1_ase_a = ms_h1_ase[dds.ms.h1_a$rowname, ] ##filtering genes in inheritance file
ms_h1_rd_inht = cbind(dds.ms.h1_a, ms_h1_ase_a) %>% select("inheritance.ms", "type", "genes")
ms_h1_rd_inht_frq = as.data.frame(table(ms_h1_rd_inht$inheritance.ms, ms_h1_rd_inht$type))
ms_h1_rd_inht_frq$Sample = "MS_H1"

######
##H2_MS##
######

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R

#inheritance
View(dds.ms.h2)
#allele specific expression
View(ms_h2_ase)

##getting only gene names and inheritance for each gene
dds.ms.h2_a = dds.ms.h2
ms_h2_ase_a = ms_h2_ase[dds.ms.h2_a$rowname, ] ##filtering genes in inheritance file
ms_h2_rd_inht = cbind(dds.ms.h2_a, ms_h2_ase_a) %>% select("inheritance.ms", "type", "genes")
ms_h2_rd_inht_frq = as.data.frame(table(ms_h2_rd_inht$inheritance.ms, ms_h2_rd_inht$type))
ms_h2_rd_inht_frq$Sample = "MS_H2"


######
##H2_MG##
######

###VARIABLES SAVED FROM 2.Classifying Inheritance.R

#inheritance
View(dds.mg.h2)
#allele specific expression
View(mg_h2_ase)

##getting only gene names and inheritance for each gene
dds.mg.h2_a = dds.mg.h2
mg_h2_ase_a = mg_h2_ase[dds.mg.h2_a$rowname, ] ##filtering genes in inheritance file
mg_h2_rd_inht = cbind(dds.mg.h2_a, mg_h2_ase_a) %>% select("inheritance.mg", "type", "genes")
mg_h2_rd_inht_frq = as.data.frame(table(mg_h2_rd_inht$inheritance.mg, mg_h2_rd_inht$type))
mg_h2_rd_inht_frq$Sample = "MG_H2"


###COMBING ALL THE DATA TABLES

inhtxregdiv_all = rbind(fg_h1_rd_inht_frq, fg_h2_rd_inht_frq, fs_h1_rd_inht_frq, fs_h2_rd_inht_frq,
                        ms_h1_rd_inht_frq, ms_h2_rd_inht_frq, mg_h2_rd_inht_frq)

inhtxregdiv_all_a = inhtxregdiv_all %>% filter(Var2 != "ambiguous" & Var2 != "conserved" & Var1 != "ambiguous" & Var1 != "no change")
inhtxregdiv_all_b = inhtxregdiv_all %>% filter(Var2 != "ambiguous" & Var1 != "ambiguous" & Var1 != "no change")



####BARPLOTS 
ggplot(inhtxregdiv_all_a, aes(x = inhtxregdiv_all_a$Var1, y = inhtxregdiv_all_a$Freq, fill = inhtxregdiv_all_a$Var2))+
  geom_bar(stat = "identity")+
  facet_rep_wrap(~Sample)

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/inhtvsregdiv_barplots.pdf", dpi = 300, width = 20, height = 8 )


#######HEATMAP 
ggplot(inhtxregdiv_all_b, aes(x= Var2, y=Var1, fill=Freq)) +
  geom_tile() +
  xlab("") + ylab("") +
  #scale_fill_distiller(palette = "RdPu", direction=-1)
  scale_fill_gradient(low="white", high="black")+ 
  #viridis::scale_fill_viridis(option = "B", direction=-1, name="number of \ngenes")+
  facet_rep_grid(~Sample) +
  geom_text(aes(label=Freq)) +
  #scale_y_discrete(labels=type.legend[3:7], limits=type.levels[3:7]) +
  #scale_x_discrete(labels=class.legend[3:7][c(4:5,1:3)], limits=class.levels[3:7][c(4:5,1:3)]) +
  theme(strip.background = element_rect(fill="grey90", color=NA),
        strip.text.y = element_text(size=10),
        legend.title=element_text(size=10),
        legend.text =element_text(size=10),
        axis.text.y = element_text(size=10),
        #axis.text.x = element_blank())
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/inhtvsregdiv_heatmap_b.pdf", dpi = 300, width = 35, height = 8 )



##############################################################################################################################################
#FIG 11
##############################################################################################################################################


##############################################################################################################################################
###CREATING BUBBLE PLOTS BETWEEN MS AND WM
###################################################################################################################################
##H1
wm_h1_ase_3 = rbind(wm_h1_ase_2, dds.wm.h1_X_1) ### data with updated X-linked regulatory categories
ms_h1_ase_3 = rbind(ms_h1_ase_2, dds.ms.h1_X_1) ### data with updated X-linked regulatory categories

head(wm_h1_ase)
table(wm_h1_ase_3$reg.div)
head(ms_h1_ase_3)
table(ms_h1_ase_3$reg.div)

`%notin%` <- Negate(`%in%`)

##check to see if dataframes have gene names 
head(ms_h1_ase_3)
head(wm_h1_ase_3)
nrow(ms_h1_ase_3)
nrow(wm_h1_ase_3)

##If they already have rownames as column names
dds.ms.h1.1_ase = ms_h1_ase_3 %>% dplyr::rename(genes = rowname)
dds.wm.h1.1_ase = wm_h1_ase_3 %>% dplyr::rename(genes = rowname)

nrow(dds.ms.h1.1_ase)

###filtering genes to retain same genes in boht ms and wm in H1
dds.ms.h1.1_ase = dds.ms.h1.1_ase[order((dds.ms.h1.1_ase$genes)), ]
dds.ms.h1.1_ase = subset(dds.ms.h1.1_ase, dds.ms.h1.1_ase$genes %in% dds.wm.h1.1_ase$genes)

dds.wm.h1.1_ase = dds.wm.h1.1_ase[order((dds.wm.h1.1_ase$genes)), ]
dds.wm.h1.1_ase = subset(dds.wm.h1.1_ase, dds.wm.h1.1_ase$genes %in% dds.ms.h1.1_ase$genes)%>%
  dplyr::rename(genes.wm = genes, reg.div.wm = reg.div, chr_info.wm = chr_info)

dim(dds.ms.h1.1_ase)
dim(dds.wm.h1.1_ase)

head(dds.wm.h1.1_ase)
###joining both dataframes

dds.mswm.h1_ase = cbind(dds.ms.h1.1_ase, dds.wm.h1.1_ase) 


##creating intersection in gene categories

ase.wmms.h1.intersection = as.data.frame(table(dds.mswm.h1_ase$reg.div, dds.mswm.h1_ase$reg.div.wm)) %>%
  dplyr::rename(MS = Var1, WM =Var2, Intersection = Freq)

write.csv(ase.wmms.h1.intersection, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/wm_ms_h1_WITHAMBIGUOUS.csv")

sum(ase.wmms.h1.intersection$Intersection) #11120 before ambiguous REMOVED, 6515 after removing ambiguous

ase.wmms.h1.intersection = ase.wmms.h1.intersection %>% filter(MS != "ambiguous" & WM != "ambiguous" & Intersection > 0)
ase.wmms.h1.intersection$Sample = "H1"

ggplot(data = ase.wmms.h1.intersection, aes(x=MS, y=WM, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  #scale_size(range = c(2, 20), name="Number of Genes")+
  scale_size_continuous("Number of genes", breaks=c(50, 100, 500,1000,1500,2000,2500,3000, 4000,5000),labels=c(50, 100, 500,1000,1500,2000,2500,3000, 4000,5000),range = c(1,15))



ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/ms.wm.h1.regdiv.pdf", dpi = 300, width = 18, height = 8 )



#############################################################################
##H2
##########################################################################
wm_h2_ase_3 = rbind(wm_h2_ase_2, dds.wm.h2_X_1)
ms_h2_ase_3 = rbind(ms_h2_ase_2, dds.ms.h2_X_1)

head(wm_h2_ase_3)
table(wm_h2_ase$type)

`%notin%` <- Negate(`%in%`)

##check to see if dataframes have gene names 

head(wm_h2_ase_3)
head(ms_h2_ase_3)
nrow(wm_h2_ase_3)
nrow(ms_h2_ase_3)

##If they already have rownames as column names
dds.ms.h2.1_ase = ms_h2_ase_3 %>% dplyr::rename(genes = rowname)
dds.wm.h2.1_ase = wm_h2_ase_3 %>% dplyr::rename(genes = rowname)
nrow(dds.wm.h2.1_ase)

###filtering genes to retain same genes in boht ms and wm in h2
dds.ms.h2.1_ase = dds.ms.h2.1_ase[order((dds.ms.h2.1_ase$genes)), ]
dds.ms.h2.1_ase = subset(dds.ms.h2.1_ase, dds.ms.h2.1_ase$genes %in% dds.wm.h2.1_ase$genes)

dds.wm.h2.1_ase = dds.wm.h2.1_ase[order((dds.wm.h2.1_ase$genes)), ]
dds.wm.h2.1_ase = subset(dds.wm.h2.1_ase, dds.wm.h2.1_ase$genes %in% dds.ms.h2.1_ase$genes)%>%
  dplyr::rename(genes.wm = genes, reg.div.wm = reg.div, chr_info.wm = chr_info)

dim(dds.ms.h2.1_ase)
dim(dds.wm.h2.1_ase)

head(dds.wm.h2.1_ase)
###joining both dataframes

dds.mswm.h2_ase = cbind(dds.ms.h2.1_ase, dds.wm.h2.1_ase) 


##creating intersection in gene categories

ase.wmms.h2.intersection = as.data.frame(table(dds.mswm.h2_ase$reg.div, dds.mswm.h2_ase$reg.div.wm)) %>%
  dplyr::rename(MS = Var1, WM =Var2, Intersection = Freq)

write.csv(ase.wmms.h2.intersection, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/wm_ms_h2_WITHAMBIGUOUS.csv")

sum(ase.wmms.h2.intersection$Intersection) #13734 before ambiguous REMOVED, 6515 after removing ambiguous

ase.wmms.h2.intersection = ase.wmms.h2.intersection %>% filter(MS != "ambiguous" & WM != "ambiguous" & Intersection > 0)
ase.wmms.h2.intersection$Sample = "h2"

ggplot(data = ase.wmms.h2.intersection, aes(x=MS, y=WM, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  #scale_size(range = c(2, 20), name="Number of Genes")+
  scale_size_continuous("Number of genes", breaks=c(50, 100, 500,1000,1500,2000,2500,3000, 4000,5000),labels=c(50, 100, 500,1000,1500,2000,2500,3000, 4000,5000),range = c(1,15))



ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 5/ms.wmh2.regdiv.pdf", dpi = 300, width = 18, height = 8 )







