#######Figures looking at expression dominance/inheritance patterns


####SCATTERPLOTS FOR INHERITANCE CATEGORIES PER SAMPLE

#H1

#ms

ms_h1_xx = (cbind(dds.ms.cre.clat.res["ms.cre.clat"], dds.ms.cre.H1.res["ms.cre.H1"], dds.ms.clat.H1.res["ms.clat.H1"], 
           dds.ms.cre.H1.res["log2FoldChange"])) %>% rename(Fc_h1vscre = "log2FoldChange")
ms_h1_xx = cbind(ms_h1_xx,dds.ms.clat.H1.res["log2FoldChange"]) %>% rename(Fc_h1vsclat = "log2FoldChange")

ms_h1_xx = na.omit(ms_h1_xx)
ms_h1_xx$inheritance.ms = classify_inheritance(ms_h1_xx)
table(ms_h1_xx$inheritance.ms)
ms_h1_xx = ms_h1_xx %>% filter(ms_h1_xx$inheritance.ms != "ambiguous")

ggplot(ms_h1_xx, aes(ms_h1_xx$Fc_h1vscre, ms_h1_xx$Fc_h1vsclat, color = ms_h1_xx$inheritance.ms)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept=0, linetype=2) +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H1"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H1"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_ms_h1_scatterplot.pdf", dpi = 300, width = 12, height = 8 )



###H2 MS

ms_h2_xx = (cbind(dds.ms.cre.clat.res["ms.cre.clat"], dds.ms.cre.H2.res["ms.cre.H2"], dds.ms.clat.H2.res["ms.clat.H2"], 
                  dds.ms.cre.H2.res["log2FoldChange"])) %>% rename(Fc_h2vscre = "log2FoldChange")
ms_h2_xx = cbind(ms_h2_xx,dds.ms.clat.H2.res["log2FoldChange"]) %>% rename(Fc_h2vsclat = "log2FoldChange")

ms_h2_xx = na.omit(ms_h2_xx)
ms_h2_xx$inheritance.ms = classify_inheritance(ms_h2_xx)
table(ms_h2_xx$inheritance.ms)
table(dds.ms.h2$inheritance.ms) #verifying
ms_h2_xx = ms_h2_xx %>% filter(ms_h2_xx$inheritance.ms != "ambiguous")

ggplot(ms_h2_xx, aes(ms_h2_xx$Fc_h2vscre, ms_h2_xx$Fc_h2vsclat, color = ms_h2_xx$inheritance.ms)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H2"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H2"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_ms_h2_scatterplot.pdf", dpi = 300, width = 12, height = 8 )


##MG H2
mg_h2_xx = (cbind(dds.mg.cre.clat.res["mg.cre.clat"], dds.mg.cre.H2.res["mg.cre.H2"], dds.mg.clat.H2.res["mg.clat.H2"], 
                  dds.mg.cre.H2.res["log2FoldChange"])) %>% rename(Fc_h2vscre = "log2FoldChange")
mg_h2_xx = cbind(mg_h2_xx,dds.mg.clat.H2.res["log2FoldChange"]) %>% rename(Fc_h2vsclat = "log2FoldChange")

mg_h2_xx = na.omit(mg_h2_xx)
mg_h2_xx$inheritance.mg = classify_inheritance(mg_h2_xx)
table(mg_h2_xx$inheritance.mg)
table(dds.mg.h2$inheritance.mg) #verifying
mg_h2_xx = mg_h2_xx %>% filter(mg_h2_xx$inheritance.mg != "ambiguous")

ggplot(mg_h2_xx, aes(mg_h2_xx$Fc_h2vscre, mg_h2_xx$Fc_h2vsclat, color = mg_h2_xx$inheritance.mg)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H2"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H2"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_mg_h2_scatterplot.pdf", dpi = 300, width = 12, height = 8 )

####FG H1

fg_h1_xx = (cbind(dds.fg.cre.clat.res["fg.cre.clat"], dds.fg.cre.H1.res["fg.cre.H1"], dds.fg.clat.H1.res["fg.clat.H1"], 
                  dds.fg.cre.H1.res["log2FoldChange"])) %>% rename(Fc_h1vscre = "log2FoldChange")
fg_h1_xx = cbind(fg_h1_xx,dds.fg.clat.H1.res["log2FoldChange"]) %>% rename(Fc_h1vsclat = "log2FoldChange")

fg_h1_xx = na.omit(fg_h1_xx)
fg_h1_xx$inheritance.fg = classify_inheritance(fg_h1_xx)
table(fg_h1_xx$inheritance.fg)
fg_h1_xx = fg_h1_xx %>% filter(fg_h1_xx$inheritance.fg != "ambiguous")

ggplot(fg_h1_xx, aes(fg_h1_xx$Fc_h1vscre, fg_h1_xx$Fc_h1vsclat, color = fg_h1_xx$inheritance.fg)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept=0, linetype=2) +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H1"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H1"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_fg_h1_scatterplot.pdf", dpi = 300, width = 12, height = 8 )

####FG H

fg_h2_xx = (cbind(dds.fg.cre.clat.res["fg.cre.clat"], dds.fg.cre.H2.res["fg.cre.H2"], dds.fg.clat.H2.res["fg.clat.H2"], 
                  dds.fg.cre.H2.res["log2FoldChange"])) %>% rename(Fc_h2vscre = "log2FoldChange")
fg_h2_xx = cbind(fg_h2_xx,dds.fg.clat.H2.res["log2FoldChange"]) %>% rename(Fc_h2vsclat = "log2FoldChange")

fg_h2_xx = na.omit(fg_h2_xx)
fg_h2_xx$inheritance.fg = classify_inheritance(fg_h2_xx)
table(fg_h2_xx$inheritance.fg)
fg_h2_xx = fg_h2_xx %>% filter(fg_h2_xx$inheritance.fg != "ambiguous")

ggplot(fg_h2_xx, aes(fg_h2_xx$Fc_h2vscre, fg_h2_xx$Fc_h2vsclat, color = fg_h2_xx$inheritance.fg)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept=0, linetype=2) +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H2"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H2"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_fg_h2_scatterplot.pdf", dpi = 300, width = 12, height = 8 )


##FS H1


fs_h1_xx = (cbind(dds.fs.cre.clat.res["fs.cre.clat"], dds.fs.cre.H1.res["fs.cre.H1"], dds.fs.clat.H1.res["fs.clat.H1"], 
                  dds.fs.cre.H1.res["log2FoldChange"])) %>% rename(Fc_h1vscre = "log2FoldChange")
fs_h1_xx = cbind(fs_h1_xx,dds.fs.clat.H1.res["log2FoldChange"]) %>% rename(Fc_h1vsclat = "log2FoldChange")

fs_h1_xx = na.omit(fs_h1_xx)
fs_h1_xx$inheritance.fs = classify_inheritance(fs_h1_xx)
table(fs_h1_xx$inheritance.fs)
fs_h1_xx = fs_h1_xx %>% filter(fs_h1_xx$inheritance.fs != "ambiguous")

ggplot(fs_h1_xx, aes(fs_h1_xx$Fc_h1vscre, fs_h1_xx$Fc_h1vsclat, color = fs_h1_xx$inheritance.fs)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept=0, linetype=2) +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H1"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H1"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_fs_h1_scatterplot.pdf", dpi = 300, width = 12, height = 8 )

####fs H2

fs_h2_xx = (cbind(dds.fs.cre.clat.res["fs.cre.clat"], dds.fs.cre.H2.res["fs.cre.H2"], dds.fs.clat.H2.res["fs.clat.H2"], 
                  dds.fs.cre.H2.res["log2FoldChange"])) %>% rename(Fc_h2vscre = "log2FoldChange")
fs_h2_xx = cbind(fs_h2_xx,dds.fs.clat.H2.res["log2FoldChange"]) %>% rename(Fc_h2vsclat = "log2FoldChange")

fs_h2_xx = na.omit(fs_h2_xx)
fs_h2_xx$inheritance.fs = classify_inheritance(fs_h2_xx)
table(fs_h2_xx$inheritance.fs)
fs_h2_xx = fs_h2_xx %>% filter(fs_h2_xx$inheritance.fs != "ambiguous")

ggplot(fs_h2_xx, aes(fs_h2_xx$Fc_h2vscre, fs_h2_xx$Fc_h2vsclat, color = fs_h2_xx$inheritance.fs)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept=0, linetype=2) +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H2"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H2"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_fs_h2_scatterplot.pdf", dpi = 300, width = 12, height = 8 )




##H1 WM

wm_h1_xx = (cbind(dds.wm.cre.clat.res["wm.cre.clat"], dds.wm.cre.H1.res["wm.cre.H1"], dds.wm.clat.H1.res["wm.clat.H1"], 
                  dds.wm.cre.H1.res["log2FoldChange"])) %>% rename(Fc_h1vscre = "log2FoldChange")
wm_h1_xx = cbind(wm_h1_xx,dds.wm.clat.H1.res["log2FoldChange"]) %>% rename(Fc_h1vsclat = "log2FoldChange")

wm_h1_xx = na.omit(wm_h1_xx)
wm_h1_xx$inheritance.wm = classify_inheritance(wm_h1_xx)
table(wm_h1_xx$inheritance.wm)
table(dds.wm.h1$inheritance.wm)
wm_h1_xx = wm_h1_xx %>% filter(wm_h1_xx$inheritance.wm != "ambiguous")

ggplot(wm_h1_xx, aes(wm_h1_xx$Fc_h1vscre, wm_h1_xx$Fc_h1vsclat, color = wm_h1_xx$inheritance.wm)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H1"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H1"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_wm_h1_scatterplot.pdf", dpi = 300, width = 12, height = 8 )

##H2 WM

wm_h2_xx = (cbind(dds.wm.cre.clat.res["wm.cre.clat"], dds.wm.cre.H2.res["wm.cre.H2"], dds.wm.clat.H2.res["wm.clat.H2"], 
                  dds.wm.cre.H2.res["log2FoldChange"])) %>% rename(Fc_h2vscre = "log2FoldChange")
wm_h2_xx = cbind(wm_h2_xx,dds.wm.clat.H2.res["log2FoldChange"]) %>% rename(Fc_h2vsclat = "log2FoldChange")

wm_h2_xx = na.omit(wm_h2_xx)
wm_h2_xx$inheritance.wm = classify_inheritance(wm_h2_xx)
table(wm_h2_xx$inheritance.wm)
table(dds.wm.h2$inheritance.wm)
wm_h2_xx = wm_h2_xx %>% filter(wm_h2_xx$inheritance.wm != "ambiguous")

ggplot(wm_h2_xx, aes(wm_h2_xx$Fc_h2vscre, wm_h2_xx$Fc_h2vsclat, color = wm_h2_xx$inheritance.wm)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Set2") +
  xlim(-10, 20) +
  ylim(-10, 20) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste( " [log"[2],"(",italic("H2"),"/",italic("Cre"),")]", sep=''))) +
  xlab(expression(paste("[log"[2],"(",italic("H2"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_wm_h2_scatterplot.pdf", dpi = 300, width = 12, height = 8 )


#############
##Fig2A
###Barplot of inheritance divergence frequencies across samples and hybrids

##plotting proportions across chromosomes
library(ggstats)

###subsetting to get rid of ambisuous counts
inheritance_all_counts$sample = paste(inheritance_all_counts$Sex, inheritance_all_counts$Tissue)
inheritance_all_counts_a = inheritance_all_counts %>% filter(inheritance_all_counts$Var1 != "ambiguous") %>% select(-Proportion)
View(inheritance_all_counts)

###Plotting proportion of genes in each category in each sample in both hybrids
ggplot(inheritance_all_counts_a) +
  aes(x = inheritance_all_counts_a$sample, fill = inheritance_all_counts_a$Var1, weight = inheritance_all_counts_a$Freq, by = as.factor(inheritance_all_counts_a$sample)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Inheritance Pattern (Actual")

#ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/fig2_inheritancepattern.pdf", dpi = 300, width = 18, height = 8 )

##inheritance plot for WM 

inheritance_all_counts_wm = inheritance_all_counts_a %>% filter(inheritance_all_counts_a$Tissue == "W")
View(inheritance_all_counts_wm)

ggplot(inheritance_all_counts_wm) +
  aes(x = inheritance_all_counts_wm$sample, fill = inheritance_all_counts_wm$Var1, weight = inheritance_all_counts_wm$Freq, by = as.factor(inheritance_all_counts_wm$sample)) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Inheritance Pattern WM")

ggplot(inheritance_all_counts_wm) +
  aes(x = inheritance_all_counts_wm$Var1, y = inheritance_all_counts_wm$Freq, fill = inheritance_all_counts_wm$Var1) +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~Hybrid)+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))+
  ggtitle("Inheritance Pattern for WM")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inht_barplot_WM.pdf", dpi = 300, width = 12, height = 8 )


##########################################################################################################3
######FIG 2B########################################################################################

####VISUALIZING INHERITANCE BETWEEN MALES AND FEMALES OF H2

###Proportion of genes between males and females 
inht_mgxfg_h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inheritance_mgxfg_h2.txt", sep="\t", head=T, comment.char="#")
View(inht_mgxfg_h2)

ggplot(inht_mgxfg_h2, aes(x=inht_mgxfg_h2$Female_Percentage, y=inht_mgxfg_h2$Male_Percentage, color = inht_mgxfg_h2$Inheritance.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 17) +
  scale_color_brewer(palette = "Set2") +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Male Gonad")))) +
  xlab(expression(paste(bold("Female Gonad")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H2 MGXFG")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/fgxmg_h2.pdf", dpi = 300, width = 10, height = 8 )

##########################################################################################################3
######FIG 2C########################################################################################

###VISUALIZING INHERITANCE BETWEEN MS X FS FOR H1 and H2

fsxms_h1_h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inheritance_msxfs_h1h2.txt", sep="\t", head=T, comment.char="#")
View(fsxms_h1_h2)

ggplot(fsxms_h1_h2, aes(x=fsxms_h1_h2$Female_Percentage, y=fsxms_h1_h2$Male_Percentage, color = fsxms_h1_h2$Inheritance.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8) +
  scale_color_brewer(palette = "Set2") +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Male Soma")))) +
  xlab(expression(paste(bold("Female Soma")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H2xH1 MSXFS")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/fsxms_h1h2.pdf", dpi = 300, width = 10, height = 8 )



###VISUALIZING INHERITANCE ACROSS MALE GONAD AND MALE SOMA IN H2

inht_mgxms_h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inheritance_mgxms_h2.txt", sep="\t", head=T, comment.char="#")
View(inht_mgxms_h2)


ggplot(inht_mgxms_h2, aes(x=inht_mgxms_h2$Gonad_Percentage, y=inht_mgxms_h2$Soma_Percentage, color = inht_mgxms_h2$Inheritance.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 17) +
  scale_color_brewer(palette = "Set2") +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Male Soma")))) +
  xlab(expression(paste(bold("Male Gonad")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H2 MGXMS")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/mgxms_h2.pdf", dpi = 300, width = 10, height = 8 )


###VISUALIZING INHERITANCE IN FG X FS IN H1 AND H2

fgxfs_h1h2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/inheritance_fgxfs_h1H2.txt", sep="\t", head=T, comment.char="#")

View(fgxfs_h1h2)
ggplot(fgxfs_h1h2, aes(x=fgxfs_h1h2$Gonad_Percentage, y=fgxfs_h1h2$Soma_Percentage, color = fgxfs_h1h2$Inheritance.Category, shape = Hybrid)) +
  geom_point(alpha=0.8, show.legend=T, size=8) +
  scale_color_brewer(palette = "Set2") +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("Female Soma")))) +
  xlab(expression(paste(bold("Female Gonad")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  ggtitle("H1XH2 FGXFS")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/fgxfs_h1h2.pdf", dpi = 300, width = 10, height = 8 )



##################################################################################################################

#####FIG 2C
##identifying if sex-biased genes overlap with transgressive genes (underdominant and overdominant genes)
###For example - maybe male-biased genes are upregulated in females and shwo up as overdominant

####1. Run 1.WT_general DGEresult to get sex-biased genes 

####Importing file with sex-biased genes
sex_biased_genes = read.csv(file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 1/mxf_res_genes.csv", header = TRUE,
                             comment.char = "#")
rownames(sex_biased_genes) = sex_biased_genes$X
View(sex_biased_genes)

###################################################################################################
####Inheritance for FG,FS,MG,MS H2 from 2.Classifying inheritance.R
###################################################################################################

#########
###H2_FG#
#########
dds.fg.h2 = rownames_to_column(dds.fg.h2)
rownames(dds.fg.h2) = dds.fg.h2$rowname

#filtering larger sex-biased genes to get fg genes
sex_biased_genes_fg_h2 = sex_biased_genes[dds.fg.h2$rowname, ]
sex_biased_genes_fg_h2 = na.omit(sex_biased_genes_fg_h2)
dds.fg.h2_2 = dds.fg.h2[sex_biased_genes_fg_h2$X, ]
fg_h2_sb = cbind(sex_biased_genes_fg_h2, dds.fg.h2_2) %>% select("X", "mxf", "rowname", "inheritance.fg")
fg_h2_sb_table = as.data.frame(table(fg_h2_sb$mxf, fg_h2_sb$inheritance.fg))
fg_h2_sb_table$Sample = "FG_H2"

#########
###H2_MG#
#########
dds.mg.h2 = rownames_to_column(dds.mg.h2)
rownames(dds.mg.h2) = dds.mg.h2$rowname

#filtering larger sex-biased genes to get mg genes
sex_biased_genes_mg_h2 = sex_biased_genes[dds.mg.h2$rowname, ]
sex_biased_genes_mg_h2 = na.omit(sex_biased_genes_mg_h2)
dds.mg.h2_2 = dds.mg.h2[sex_biased_genes_mg_h2$X, ]
mg_h2_sb = cbind(sex_biased_genes_mg_h2, dds.mg.h2_2) %>% select("X", "mxf", "rowname", "inheritance.mg")
mg_h2_sb_table = as.data.frame(table(mg_h2_sb$mxf, mg_h2_sb$inheritance.mg))
mg_h2_sb_table$Sample = "MG_H2"

View(mg_h2_sb_table)

#########
###H2_FS#
#########
dds.fs.h2 = rownames_to_column(dds.fs.h2)
rownames(dds.fs.h2) = dds.fs.h2$rowname

#filtering larger sex-biased genes to get fs genes
sex_biased_genes_fs_h2 = sex_biased_genes[dds.fs.h2$rowname, ]
sex_biased_genes_fs_h2 = na.omit(sex_biased_genes_fs_h2)
dds.fs.h2_2 = dds.fs.h2[sex_biased_genes_fs_h2$X, ]
fs_h2_sb = cbind(sex_biased_genes_fs_h2, dds.fs.h2_2) %>% select("X", "mxf", "rowname", "inheritance.fs")
fs_h2_sb_table = as.data.frame(table(fs_h2_sb$mxf, fs_h2_sb$inheritance.fs))
fs_h2_sb_table$Sample = "FS_H2"

View(fs_h2_sb_table)


#########
###H2_MS#
#########
dds.ms.h2 = rownames_to_column(dds.ms.h2)
rownames(dds.ms.h2) = dds.ms.h2$rowname

#filtering larger sex-biased genes to get MS genes
sex_biased_genes_ms_h2 = sex_biased_genes[dds.ms.h2$rowname, ]
sex_biased_genes_ms_h2 = na.omit(sex_biased_genes_ms_h2)
dds.ms.h2_2 = dds.ms.h2[sex_biased_genes_ms_h2$X, ]
ms_h2_sb = cbind(sex_biased_genes_ms_h2, dds.ms.h2_2) %>% select("X", "mxf", "rowname", "inheritance.ms")
ms_h2_sb_table = as.data.frame(table(ms_h2_sb$mxf, ms_h2_sb$inheritance.ms))
ms_h2_sb_table$Sample = "ms_H2"

View(ms_h2_sb_table)

###################################################################################################33
#####InheriTance for GF, FS, MS IN H1
####################################################################################################3

#########
###H1_FG#
#########
dds.fg.h1 = rownames_to_column(dds.fg.h1)
rownames(dds.fg.h1) = dds.fg.h1$rowname
View(dds.fg.h1)

#filtering larger sex-biased genes to get fg genes
sex_biased_genes_fg_h1 = sex_biased_genes[dds.fg.h1$rowname, ]
sex_biased_genes_fg_h1 = na.omit(sex_biased_genes_fg_h1)
dds.fg.h1_2 = dds.fg.h1[sex_biased_genes_fg_h1$X, ]
fg_h1_sb = cbind(sex_biased_genes_fg_h1, dds.fg.h1_2) %>% select("X", "mxf", "rowname", "inheritance.fg")
fg_h1_sb_table = as.data.frame(table(fg_h1_sb$mxf, fg_h1_sb$inheritance.fg))
fg_h1_sb_table$Sample = "FG_H1"

View(fg_h1_sb_table)

#########
###H1_FS#
#########
#FS
dds.fs.h1 = rownames_to_column(dds.fs.h1)
rownames(dds.fs.h1) = dds.fs.h1$rowname

#filtering larger sex-biased genes to get fs genes
sex_biased_genes_fs_h1 = sex_biased_genes[dds.fs.h1$rowname, ]
sex_biased_genes_fs_h1 = na.omit(sex_biased_genes_fs_h1)
dds.fs.h1_2 = dds.fs.h1[sex_biased_genes_fs_h1$X, ]
fs_h1_sb = cbind(sex_biased_genes_fs_h1, dds.fs.h1_2) %>% select("X", "mxf", "rowname", "inheritance.fs")
fs_h1_sb_table = as.data.frame(table(fs_h1_sb$mxf, fs_h1_sb$inheritance.fs))
fs_h1_sb_table$Sample = "FS_h1"

View(fs_h1_sb_table)

#########
###H1_MS#
#########
dds.ms.h1 = rownames_to_column(dds.ms.h1)
rownames(dds.ms.h1) = dds.ms.h1$rowname


#filtering larger sex-biased genes to get MS genes
sex_biased_genes_ms_h1 = sex_biased_genes[dds.ms.h1$rowname, ]
sex_biased_genes_ms_h1 = na.omit(sex_biased_genes_ms_h1)
dds.ms.h1_2 = dds.ms.h1[sex_biased_genes_ms_h1$X, ]
ms_h1_sb = cbind(sex_biased_genes_ms_h1, dds.ms.h1_2) %>% select("X", "mxf", "rowname", "inheritance.ms")
ms_h1_sb_table = as.data.frame(table(ms_h1_sb$mxf, ms_h1_sb$inheritance.ms))
ms_h1_sb_table$Sample = "ms_h1"

View(ms_h1_sb_table)


####Merging all frequency tables

sex_bised_inht = rbind(fg_h2_sb_table, fg_h1_sb_table, fs_h2_sb_table, fs_h1_sb_table,
                       mg_h2_sb_table, ms_h1_sb_table, ms_h2_sb_table)
View(sex_bised_inht)
write.csv(sex_bised_inht, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/sex_biasedandinht.csv")



###Heatmap for sex-biased genes and inheritance

ggplot(sex_bised_inht %>% filter(Var2 != "ambiguous"), aes(x= Var2, y=Var1, fill=Freq)) +
  geom_tile() +
  xlab("") + ylab("") +
  #scale_fill_distiller(palette = "RdPu", direction=-1)
  scale_fill_gradient(low="#fff4f4", high="#3f0a0f")+ 
  #viridis::scale_fill_viridis(option = "B", direction=-1, name="number of \ngenes")+
  facet_grid(~ Sample) +
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

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/sex_biasedvsinht_heatmap.pdf", dpi = 300, width = 20, height = 8 )

##Barplot
sex_bised_inht_a = sex_bised_inht %>% filter(Var2 != "ambiguous")
ggplot(sex_bised_inht_a, aes(x = sex_bised_inht_a$Var1, y = sex_bised_inht_a$Freq, fill = sex_bised_inht_a$Var2))+
  geom_bar(stat = "identity")+
  facet_wrap(~Sample)

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/sex_bsvsinht_barplot1.pdf", dpi = 300, width = 20, height = 8 )

ggplot(sex_bised_inht_a, aes(x = sex_bised_inht_a$Var2, y = sex_bised_inht_a$Freq, fill = sex_bised_inht_a$Var1))+
  geom_bar(stat = "identity")+
  facet_wrap(~Sample)

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/sex_bsvsinht_barplot2.pdf", dpi = 300, width = 20, height = 8 )


###Line plot 
sex_bised_inht = read.csv("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/sex_biasedandinht.csv", sep = ",")
str(sex_bised_inht)
sex_bised_inht$Var1 = as.factor(sex_bised_inht$Var1)
sex_bised_inht$Var2 = as.factor(sex_bised_inht$Var2)

sex_bised_inht %>% filter(Var2 != "ambiguous" & Var2 != "no change") %>%
  ggplot(aes(x = factor(Var2, level = c("C. latens dominant", "C. remanei dominant", "Additive", "Overdominant", "Underdominant")), Freq, col = Var1, group = Var1)) +
  #geom_hline(yintercept=0, color = "black", linetype="dashed")+
  #scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  geom_text(aes(label = Freq)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~Sample) +
  #stat_smooth(method = "lm") +
  ggtitle("Sex biased vs inht")
ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 2/sex_bsvsinht_lineplot.pdf", dpi = 300, width = 20, height = 8 )

#########################################################################################################3
#######VISUALIZING INHERITANCE BETWEEN H1 AND H2
#########################################################################################################3

inht_h1xh2 = read.table("Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/inheritance_all_h1xh2.txt", sep="\t", head=T, comment.char="#")
View(inht_h1xh2)
ggplot(inht_h1xh2, aes(x=inht_h1xh2$H1_percentage, y=inht_h1xh2$H2_percentage, color = inht_h1xh2$Inheritance_Category, shape = Tissue)) +
  geom_point(alpha=0.8, show.legend=T, size=8) +
  scale_color_brewer(palette = "Set2") +
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

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/inht_h1xh2_1.pdf", dpi = 300, width = 15, height = 8 )

ggplot(inht_h1xh2, aes(x=inht_h1xh2$H1_percentage, y=inht_h1xh2$H2_percentage, color = inht_h1xh2$Inheritance_Category, shape = Sex)) +
  geom_point(alpha=0.8, show.legend=T, size=8, shape = 15) +
  #scale_shape_manual(values=c(15, 23))+
  scale_color_brewer(palette = "Set2") +
  xlim(0, 50) +
  ylim(0, 50) +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_abline(a=0, b=1)+
  ylab(expression(paste(bold("H2")))) +
  xlab(expression(paste(bold("H1")))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  facet_rep_wrap(~Tissue+Sex) +
  theme(strip.background = element_rect(fill="grey90", color=NA))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/inht_h1xh2_2.pdf", dpi = 300, width = 18, height = 8 )



#######################################################################################################
###COMAPRING H1 AND H2 

###INHERITANCE CATEGORY BETWEEN H1 AND H2

###########
###FG######
###########

####NEED TO RUN 2.CLASSIFYING INHERITANCE.R FIRST
###VARIABLES SAVED FROM 2.Classifying Inheritance.R

dim(dds.fg.h1)

`%notin%` <- Negate(`%in%`)

##chanign rownames as column
dds.fg.h1.1 = rownames_to_column(dds.fg.h1)
dds.fg.h2.1 = rownames_to_column(dds.fg.h2)

##If they already have rownames as column names
dds.fg.h1.1 = dds.fg.h1
dds.fg.h2.1 = dds.fg.h2

##ordering each dataframe adn filtering to make sure both datafrmae have the same number of rows
dds.fg.h1.1 = dds.fg.h1.1[order(dds.fg.h1.1$rowname), ]
dds.fg.h1.1 = subset(dds.fg.h1.1, dds.fg.h1.1$rowname %in% dds.fg.h2.1$rowname)

dds.fg.h2.1 = dds.fg.h2.1[order(dds.fg.h2.1$rowname), ] %>%
  rename(fg.cre.clat_2 = fg.cre.clat, rowname.h2 = rowname, inheritance.fg.h2 = inheritance.fg)


nrow(dds.fg.h1.1)
nrow(dds.fg.h2.1)

###joining both dataframes

dds.fg.h1h2 = cbind(dds.fg.h1.1, dds.fg.h2.1) 
head(dds.fg.h1h2)

#####Visualizing using scatterplot
ggplot(dds.fg.h1h2 %>% filter(inheritance.fg != "ambiguous" & inheritance.fg.h2 != "ambiguous"), aes(x = inheritance.fg, inheritance.fg.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

#####VISUALIZING USING BUBBLE PLOTS

head(as.data.frame(table(dds.fg.h1h2$inheritance.fg, dds.fg.h1h2$inheritance.fg.h2)))

##creating intersection in gene categories
fg.h1h2.intersection = as.data.frame(table(dds.fg.h1h2$inheritance.fg, dds.fg.h1h2$inheritance.fg.h2)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq)
fg.h1h2.intersection = fg.h1h2.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
fg.h1h2.intersection$Sample = "FG"
View(fg.h1h2.intersection)

ggplot(data = fg.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/fg.h1h2.inht.pdf", dpi = 300, width = 20, height = 18 )

###########
###FS######
###########

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R
head(dds.fs.h2)
nrow(dds.fs.h2)

head(dds.fs.h1)
nrow(dds.fs.h1)

##chanign rownames as column
dds.fs.h1.1 = rownames_to_column(dds.fs.h1) 
dds.fs.h2.1 = rownames_to_column(dds.fs.h2)

head(dds.fs.h2.1)
##if rownames are already changed 
dds.fs.h1.1 = (dds.fs.h1)
dds.fs.h2.1 = (dds.fs.h2)

##ordering each dataframe
dds.fs.h1.1 = dds.fs.h1.1[order(dds.fs.h1.1$rowname), ]
dds.fs.h2.1 = dds.fs.h2.1[order(dds.fs.h2.1$rowname), ] %>%
  rename(fs.cre.clat_2 = fs.cre.clat, rowname.h2 = rowname, inheritance.fs.h2 = inheritance.fs)
head(dds.fs.h2.1)

###joining both dataframes
dds.fs.h1h2 = cbind(dds.fs.h1.1, dds.fs.h2.1) 
head(dds.fs.h1h2)

##scatterplot
ggplot(dds.fs.h1h2 %>% filter(inheritance.fs != "ambiguous" & inheritance.fs.h2 != "ambiguous"), aes(x = inheritance.fs, inheritance.fs.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

###BUBBLE PLOT

#creating intersection in gene categories
fs.h1h2.intersection = as.data.frame(table(dds.fs.h1h2$inheritance.fs, dds.fs.h1h2$inheritance.fs.h2)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq) %>% 
  filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
fs.h1h2.intersection$Sample = "FS"

ggplot(data = fs.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/fs.h1h2.inht.pdf", dpi = 300, width = 20, height = 18 )

###########
###MS######
###########

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R

head(dds.ms.h1)
nrow(dds.ms.h1)

head(dds.ms.h2)
nrow(dds.ms.h2)

#chanign rownames to columns
dds.ms.h1.1 = rownames_to_column(dds.ms.h1) 

dds.ms.h2.1 = rownames_to_column(dds.ms.h2)

##if rownames are already changed
dds.ms.h1.1 = dds.ms.h1
dds.ms.h2.1 = dds.ms.h2
#ordering both dataframes
dds.ms.h1.1 = dds.ms.h1.1[order(dds.ms.h1.1$rowname), ]
dds.ms.h2.1 = dds.ms.h2.1[order(dds.ms.h2.1$rowname), ] %>%
  rename(ms.cre.clat_2 = ms.cre.clat, rowname.h2 = rowname, inheritance.ms.h2 = inheritance.ms) 
colnames(dds.ms.h2.1)

#joining dataframes
dds.ms.h1h2 = cbind(dds.ms.h1.1, dds.ms.h2.1) 
head(dds.ms.h1h2)

##scatterplot
ggplot(dds.ms.h1h2 %>% filter(inheritance.ms != "ambiguous" & inheritance.ms.h2 != "ambiguous"), aes(x = inheritance.ms, inheritance.ms.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

###BUBBLE PLOT

#creating intersection between gene categories
ms.h1h2.intersection = as.data.frame(table(dds.ms.h1h2$inheritance.ms, dds.ms.h1h2$inheritance.ms.h2)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq) %>% 
  filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
ms.h1h2.intersection$Sample = "MS"
ggplot(data = ms.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/ms.h1h2.inht.pdf", dpi = 300, width = 20, height = 18 )


###########
###WM######
###########

###VARIABLES SAVED FROM 2.Classifying Inheritance.R
head(dds.wm.h1)
nrow(dds.wm.h1)
head(dds.wm.h2)
nrow(dds.wm.h2)

#changing rownames to columnames
dds.wm.h1.1 = rownames_to_column(dds.wm.h1) 
dds.wm.h2.1 = rownames_to_column(dds.wm.h2)

#if already done
dds.wm.h1.1 = (dds.wm.h1) 
dds.wm.h2.1 = (dds.wm.h2)

nrow(dds.wm.h1.1 %>% filter(rowname %in% dds.wm.h2.1$rowname))

#ordering both dataframes
dds.wm.h1.1 = dds.wm.h1.1[order(dds.wm.h1.1$rowname), ]
dds.wm.h2.1 = dds.wm.h2.1[order(dds.wm.h2.1$rowname), ] %>%
  rename(wm.cre.clat_2 = wm.cre.clat, rowname.h2 = rowname, inheritance.wm.h2 = inheritance.wm) 
colnames(dds.wm.h2.1)

#joining dataframes
dds.wm.h1h2 = cbind(dds.wm.h1.1, dds.wm.h2.1) 
head(dds.wm.h1h2)

##scatterplot
ggplot(dds.wm.h1h2 %>% filter(inheritance.wm != "ambiguous" & inheritance.wm.h2 != "ambiguous"), aes(x = inheritance.wm, inheritance.wm.h2)) +
  geom_jitter(show.legend=T, alpha=0.5, color = "#4D4545", width = 0.2, height = 0.25)+
  xlab("H1")+
  ylab("H2")+
  #labs(title = "Downregulated")+
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = rel(0.7), face = "bold"))+
  theme(axis.text.y = element_text(size = 6, face = "bold"))

###BUBBLE PLOT

##creating intersection between gene categories
wm.h1h2.intersection = as.data.frame(table(dds.wm.h1h2$inheritance.wm, dds.wm.h1h2$inheritance.wm.h2)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq) %>% 
  filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)
wm.h1h2.intersection$Sample = "FG"
ggplot(data = wm.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/wm.h1h2.inht.pdf", dpi = 300, width = 20, height = 18 )


##########################################################################################33
###MAKING CUMULATIVE BUBBLEPLOT TO MAKE ALL UNITS/SCALE FOR LEGEND SAME

##COMMBINING ALL INTERSECTION DATA
inht.h1h2.intersection = rbind(fg.h1h2.intersection,fs.h1h2.intersection, ms.h1h2.intersection)

ggplot(data = inht.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  facet_wrap(~Sample)+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size_continuous("Number of genes", breaks=c(500,1000,1500,2000,2500,3000),labels=c(500,1000,1500,2000,2500,3000),range = c(1,15))
  
ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 4/all.h1h2.inht.pdf", dpi = 300, width = 33, height = 10 )

ggplot(data = inht.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  facet_wrap(~Sample)+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size_continuous("Number of genes", breaks=c(5, 50, 250, 500,1000,1500,2000,2500,3000),labels=c(5, 50, 250, 500,1000,1500,2000,2500,3000), range = c(1,15))
 