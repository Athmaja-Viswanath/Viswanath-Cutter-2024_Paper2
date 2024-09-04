

#######################################################################################################
#####FIGURE 4.9 AND SUPPLEMENTAL
############################################################################################3


########################################################################################################
#####Looking at regulatory divergence underlying transgressive genes
############################################################################################3

######
##H1_FG##
######

###VARIABLES SAVED FROM 2.Classifying Inheritance.R

#inheritance
View(dds.fg.h1)


##autosomal
#allele specific expression
View(fg_h1_ase_autosome)

##getting only gene names and inheritance for each gene
dds.fg.h1_a = rownames_to_column(dds.fg.h1)
fg_h1_ase_a = fg_h1_ase_autosome[dds.fg.h1_a$rowname, ] ##filtering genes in inheritance file
fg_h1_rd_inht_A = cbind(dds.fg.h1_a, fg_h1_ase_a) %>% select("inheritance.fg", "type", "genes")
fg_h1_rd_inht_frq_A = as.data.frame(table(fg_h1_rd_inht_A$inheritance.fg, fg_h1_rd_inht_A$type))
fg_h1_rd_inht_frq_A$Sample = "FG_H1"

head(fg_h1_ase_a)
##x-linked

colnames(fg_h1_ase_X)

##getting only gene names and inheritance for each gene
dds.fg.h1_a = rownames_to_column(dds.fg.h1)
fg_h1_ase_ax = subset(fg_h1_ase_X, fg_h1_ase_X$genes %in% dds.fg.h1_a$rowname) ##filtering genes in inheritance file
dds.fg.h1_ax = subset(dds.fg.h1_a, dds.fg.h1_a$rowname %in% fg_h1_ase_ax$genes)
fg_h1_rd_inht_X = cbind(dds.fg.h1_ax, fg_h1_ase_ax) %>% select("inheritance.fg", "combined_regdiv", "genes")
fg_h1_rd_inht_frq_X = as.data.frame(table(fg_h1_rd_inht_X$inheritance.fg, fg_h1_rd_inht_X$combined_regdiv))
fg_h1_rd_inht_frq_X$Sample = "FG_H1"

colnames(cbind(dds.fg.h1_ax, fg_h1_ase_ax))
nrow(dds.fg.h1_a)

########
##H2 FG
########

#inheritance
dds.fg.h2
###AUTOSOME
fg_h2_ase_autosome

##getting only gene names and inheritance for each gene
dds.fg.h2_a = rownames_to_column(dds.fg.h2)
fg_h2_ase_a = fg_h2_ase_autosome[dds.fg.h2_a$rowname, ] ##filtering genes in inheritance file
fg_h2_rd_inht_A = cbind(dds.fg.h2_a, fg_h2_ase_a) %>% select("inheritance.fg", "type", "genes")
fg_h2_rd_inht_frq_A = as.data.frame(table(fg_h2_rd_inht_A$inheritance.fg, fg_h2_rd_inht_A$type))
fg_h2_rd_inht_frq_A$Sample = "FG_H2"

##x-linked
colnames(fg_h2_ase_X)

##getting only gene names and inheritance for each gene
dds.fg.h2_a = rownames_to_column(dds.fg.h2)
fg_h2_ase_ax = subset(fg_h2_ase_X, fg_h2_ase_X$genes %in% dds.fg.h2_a$rowname) ##filtering genes in inheritance file
dds.fg.h2_ax = subset(dds.fg.h2_a, dds.fg.h2_a$rowname %in% fg_h2_ase_ax$genes)
fg_h2_rd_inht_X = cbind(dds.fg.h2_ax, fg_h2_ase_ax) %>% select("inheritance.fg", "combined_regdiv", "genes")
fg_h2_rd_inht_frq_X = as.data.frame(table(fg_h2_rd_inht_X$inheritance.fg, fg_h2_rd_inht_X$combined_regdiv))
fg_h2_rd_inht_frq_X$Sample = "FG_H2"


########
##H1 FS
########

###AUTOSOME

#inheritance
View(dds.fs.h1)

##autosomal
#allele specific expression
View(fs_h1_ase_autosome)

##getting only gene names and inheritance for each gene
dds.fs.h1_a = rownames_to_column(dds.fs.h1)
fs_h1_ase_a = fs_h1_ase_autosome[dds.fs.h1_a$rowname, ] ##filtering genes in inheritance file
fs_h1_rd_inht_A = cbind(dds.fs.h1_a, fs_h1_ase_a) %>% select("inheritance.fs", "type", "genes")
fs_h1_rd_inht_frq_A = as.data.frame(table(fs_h1_rd_inht_A$inheritance.fs, fs_h1_rd_inht_A$type))
fs_h1_rd_inht_frq_A$Sample = "FS_H1"

##x-linked
colnames(fs_h1_ase_X)

##getting only gene names and inheritance for each gene
dds.fs.h1_a = rownames_to_column(dds.fs.h1)
fs_h1_ase_ax = subset(fs_h1_ase_X, fs_h1_ase_X$genes %in% dds.fs.h1_a$rowname) ##filtering genes in inheritance file
dds.fs.h1_ax = subset(dds.fs.h1_a, dds.fs.h1_a$rowname %in% fs_h1_ase_ax$genes)
fs_h1_rd_inht_X = cbind(dds.fs.h1_ax, fs_h1_ase_ax) %>% select("inheritance.fs", "combined_regdiv", "genes")
fs_h1_rd_inht_frq_X = as.data.frame(table(fs_h1_rd_inht_X$inheritance.fs, fs_h1_rd_inht_X$combined_regdiv))
fs_h1_rd_inht_frq_X$Sample = "FS_H1"

########
##H2 FS
########

###AUTOSOME
#inheritance
View(dds.fs.h2)

##autosomal
#allele specific expression
View(fs_h2_ase_autosome)

##getting only gene names and inheritance for each gene
dds.fs.h2_a = rownames_to_column(dds.fs.h2)
fs_h2_ase_a = fs_h2_ase_autosome[dds.fs.h2_a$rowname, ] ##filtering genes in inheritance file
fs_h2_rd_inht_A = cbind(dds.fs.h2_a, fs_h2_ase_a) %>% select("inheritance.fs", "type", "genes")
fs_h2_rd_inht_frq_A = as.data.frame(table(fs_h2_rd_inht_A$inheritance.fs, fs_h2_rd_inht_A$type))
fs_h2_rd_inht_frq_A$Sample = "FS_H2"

##x-linked
colnames(fs_h2_ase_X)

##getting only gene names and inheritance for each gene
dds.fs.h2_a = rownames_to_column(dds.fs.h2)
fs_h2_ase_ax = subset(fs_h2_ase_X, fs_h2_ase_X$genes %in% dds.fs.h2_a$rowname) ##filtering genes in inheritance file
dds.fs.h2_ax = subset(dds.fs.h2_a, dds.fs.h2_a$rowname %in% fs_h2_ase_ax$genes)
fs_h2_rd_inht_X = cbind(dds.fs.h2_ax, fs_h2_ase_ax) %>% select("inheritance.fs", "combined_regdiv", "genes")
fs_h2_rd_inht_frq_X = as.data.frame(table(fs_h2_rd_inht_X$inheritance.fs, fs_h2_rd_inht_X$combined_regdiv))
fs_h2_rd_inht_frq_X$Sample = "FS_H2"

######
#H1 MS
#####
#AUTOSOME

#inheritance
View(dds.ms.h1)

##autosomal
#allele specific expression
View(ms_h1_ase_autosome)

##getting only gene names and inheritance for each gene
dds.ms.h1_a = rownames_to_column(dds.ms.h1)
ms_h1_ase_a = ms_h1_ase_autosome[dds.ms.h1_a$rowname, ] ##filtering genes in inheritance file
ms_h1_rd_inht_A = cbind(dds.ms.h1_a, ms_h1_ase_a) %>% select("inheritance.ms", "type", "genes")
ms_h1_rd_inht_frq_A = as.data.frame(table(ms_h1_rd_inht_A$inheritance.ms, ms_h1_rd_inht_A$type))
ms_h1_rd_inht_frq_A$Sample = "MS_H1"

##X-LINKED
colnames(dds.ms.h1_X)

##getting only gene names and inheritance for each gene
dds.ms.h1_a = rownames_to_column(dds.ms.h1)
ms_h1_ase_ax = subset(dds.ms.h1_X, dds.ms.h1_X$rowname %in% dds.ms.h1_a$rowname) ##filtering genes in inheritance file
dds.ms.h1_ax = subset(dds.ms.h1_a, dds.ms.h1_a$rowname %in% ms_h1_ase_ax$rowname)
ms_h1_rd_inht_X = cbind(dds.ms.h1_ax, ms_h1_ase_ax) %>% select("inheritance.ms", "reg.div", "rowname")
ms_h1_rd_inht_frq_X = as.data.frame(table(ms_h1_rd_inht_X$inheritance.ms, ms_h1_rd_inht_X$reg.div))
ms_h1_rd_inht_frq_X$Sample = "MS_H1"

######
#H2 MS
#####
dds.ms.h2
#AUTOSOME
ms_h2_ase_autosome

##getting only gene names and inheritance for each gene
dds.ms.h2_a = rownames_to_column(dds.ms.h2)
ms_h2_ase_a = ms_h2_ase_autosome[dds.ms.h2_a$rowname, ] ##filtering genes in inheritance file
ms_h2_rd_inht_A = cbind(dds.ms.h2_a, ms_h2_ase_a) %>% select("inheritance.ms", "type", "genes")
ms_h2_rd_inht_frq_A = as.data.frame(table(ms_h2_rd_inht_A$inheritance.ms, ms_h2_rd_inht_A$type))
ms_h2_rd_inht_frq_A$Sample = "MS_H2"

#X-LINKED
colnames(dds.ms.h2_X)

##getting only gene names and inheritance for each gene
dds.ms.h2_a = rownames_to_column(dds.ms.h2)
ms_h2_ase_ax = subset(dds.ms.h2_X, dds.ms.h2_X$rowname %in% dds.ms.h2_a$rowname) ##filtering genes in inheritance file
dds.ms.h2_ax = subset(dds.ms.h2_a, dds.ms.h2_a$rowname %in% ms_h2_ase_ax$rowname)
ms_h2_rd_inht_X = cbind(dds.ms.h2_ax, ms_h2_ase_ax) %>% select("inheritance.ms", "reg.div", "rowname")
ms_h2_rd_inht_frq_X = as.data.frame(table(ms_h2_rd_inht_X$inheritance.ms, ms_h2_rd_inht_X$reg.div))
ms_h2_rd_inht_frq_X$Sample = "MS_H2"

#########
####H2 MG
#########
dds.mg.h2
#Autosomal
mg_h2_ase_autosome

##getting only gene names and inheritance for each gene
dds.mg.h2_a = rownames_to_column(dds.mg.h2)
mg_h2_ase_a = mg_h2_ase_autosome[dds.mg.h2_a$rowname, ] ##filtering genes in inheritance file
mg_h2_rd_inht_A = cbind(dds.mg.h2_a, mg_h2_ase_a) %>% select("inheritance.mg", "type", "genes")
mg_h2_rd_inht_frq_A = as.data.frame(table(mg_h2_rd_inht_A$inheritance.mg, mg_h2_rd_inht_A$type))
mg_h2_rd_inht_frq_A$Sample = "MG_H2"

#X-linked
dds.mg.h2_X

##getting only gene names and inheritance for each gene
dds.mg.h2_a = rownames_to_column(dds.mg.h2)
mg_h2_ase_ax = subset(dds.mg.h2_X, dds.mg.h2_X$rowname %in% dds.mg.h2_a$rowname) ##filtering genes in inheritance file
dds.mg.h2_ax = subset(dds.mg.h2_a, dds.mg.h2_a$rowname %in% mg_h2_ase_ax$rowname)
mg_h2_rd_inht_X = cbind(dds.mg.h2_ax, mg_h2_ase_ax) %>% select("inheritance.mg", "reg.div", "rowname")
mg_h2_rd_inht_frq_X = as.data.frame(table(mg_h2_rd_inht_X$inheritance.mg, mg_h2_rd_inht_X$reg.div))
mg_h2_rd_inht_frq_X$Sample = "MG_H2"

########
#H1 WM##
########

dds.wm.h1
colnames(dds.wm.h2)
###AUTOSOMAL
wm_h1_ase_autosome

##getting only gene names and inheritance for each gene
#dds.wm.h1_a = rownames_to_column(dds.wm.h1)
dds.wm.h1_a = dds.wm.h1
wm_h1_ase_a = wm_h1_ase_autosome[dds.wm.h1_a$rowname, ] ##filtering genes in inheritance file
wm_h1_rd_inht_A = cbind(dds.wm.h1_a, wm_h1_ase_a) %>% select("inheritance.wm", "type", "genes")
wm_h1_rd_inht_frq_A = as.data.frame(table(wm_h1_rd_inht_A$inheritance.wm, wm_h1_rd_inht_A$type))
wm_h1_rd_inht_frq_A$Sample = "WM_H1"

###X-LINKED
colnames(dds.wm.h1_X)

##getting only gene names and inheritance for each gene
dds.wm.h1_a = dds.wm.h1
wm_h1_ase_ax = subset(dds.wm.h1_X, dds.wm.h1_X$rowname %in% dds.wm.h1_a$rowname) ##filtering genes in inheritance file
dds.wm.h1_ax = subset(dds.wm.h1_a, dds.wm.h1_a$rowname %in% wm_h1_ase_ax$rowname)
wm_h1_rd_inht_X = cbind(dds.wm.h1_ax, wm_h1_ase_ax) %>% select("inheritance.wm", "reg.div", "rowname")
wm_h1_rd_inht_frq_X = as.data.frame(table(wm_h1_rd_inht_X$inheritance.wm, wm_h1_rd_inht_X$reg.div))
wm_h1_rd_inht_frq_X$Sample = "WM_H1"

########
#H2 WM##
########

dds.wm.h2
###AUTOSOMAL

wm_h2_ase_autosome

##getting only gene names and inheritance for each gene
#dds.wm.h1_a = rownames_to_column(dds.wm.h1)
dds.wm.h2_a = dds.wm.h2
wm_h2_ase_a = wm_h2_ase_autosome[dds.wm.h2_a$rowname, ] ##filtering genes in inheritance file
wm_h2_rd_inht_A = cbind(dds.wm.h2_a, wm_h2_ase_a) %>% select("inheritance.wm", "type", "genes")
wm_h2_rd_inht_frq_A = as.data.frame(table(wm_h2_rd_inht_A$inheritance.wm, wm_h2_rd_inht_A$type))
wm_h2_rd_inht_frq_A$Sample = "WM_H2"

###X-LINKED

colnames(dds.wm.h2_X)

##getting only gene names and inheritance for each gene
dds.wm.h2_a = dds.wm.h2
wm_h2_ase_ax = subset(dds.wm.h2_X, dds.wm.h2_X$rowname %in% dds.wm.h2_a$rowname) ##filtering genes in inheritance file
dds.wm.h2_ax = subset(dds.wm.h2_a, dds.wm.h2_a$rowname %in% wm_h2_ase_ax$rowname)
wm_h2_rd_inht_X = cbind(dds.wm.h2_ax, wm_h2_ase_ax) %>% select("inheritance.wm", "reg.div", "rowname")
wm_h2_rd_inht_frq_X = as.data.frame(table(wm_h2_rd_inht_X$inheritance.wm, wm_h2_rd_inht_X$reg.div))
wm_h2_rd_inht_frq_X$Sample = "WM_H2"

###COMBING ALL THE DATA TABLES

###AUTOSOMAL DATATABLES

inhtxregdiv_Autosomes = rbind(fg_h1_rd_inht_frq_A, fg_h2_rd_inht_frq_A, fs_h1_rd_inht_frq_A, fs_h2_rd_inht_frq_A,
                        ms_h1_rd_inht_frq_A, ms_h2_rd_inht_frq_A, mg_h2_rd_inht_frq_A, wm_h1_rd_inht_frq_A, wm_h2_rd_inht_frq_A)

inhtxregdiv_Autosomes_a = inhtxregdiv_Autosomes %>% filter(Var2 != "ambiguous" & Var2 != "conserved" & Var1 != "ambiguous" & Var1 != "no change")
inhtxregdiv_all_b = inhtxregdiv_all %>% filter(Var2 != "ambiguous" & Var1 != "ambiguous" & Var1 != "no change")

#######HEATMAP 
ggplot(inhtxregdiv_Autosomes_a, aes(x= factor(Var2, level = c("cis-trans (compensatory)", "cis x trans (compensatory)", "cis + trans (enhancing)", "cis-only", "trans-only")), 
                                    y=factor(Var1, level = c("Additive", "C. latens dominant", "C. remanei dominant", "Overdominant", "Underdominant")), fill=Freq)) +
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

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/regulatory divergence figs for thesis/inhtvsregdiv_heatmap_Autosome.pdf", dpi = 300, width = 35, height = 8 )

#####X-LINKED

inhtxregdiv_X = rbind(fg_h1_rd_inht_frq_X, fg_h2_rd_inht_frq_X, fs_h1_rd_inht_frq_X, fs_h2_rd_inht_frq_X,
                              ms_h1_rd_inht_frq_X, ms_h2_rd_inht_frq_X, mg_h2_rd_inht_frq_X, wm_h1_rd_inht_frq_X, wm_h2_rd_inht_frq_X)

inhtxregdiv_X_a = inhtxregdiv_X %>% filter(Var2 != "ambiguous" & Var2 != "conserved" & Var1 != "ambiguous" & Var1 != "no change")


#######HEATMAP 
ggplot(inhtxregdiv_X_a, aes(x= factor(Var2, level = c("cis-trans (compensatory)", "cis-only", "trans-only", "Other")), 
                                    y=factor(Var1, level = c("Additive", "C. latens dominant", "C. remanei dominant", "Overdominant", "Underdominant")), fill=Freq)) +
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

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/regulatory divergence figs for thesis/inhtvsregdiv_heatmap_Xlinked.pdf", dpi = 300, width = 35, height = 8 )

