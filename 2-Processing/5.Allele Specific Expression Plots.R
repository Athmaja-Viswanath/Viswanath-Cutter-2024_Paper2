###PLOTS
library(ggplot2)
library(cowplot) #add on to ggplot for better themes and figure customization
library(lemon) #to work with legends and axes in ggplot2
library(dplyr)
library(gdata)
library(RColorBrewer)
library(colorBlindness)
library(colorspace)


######FIGURE 5 AND RELATED TO SUPPLEMENTARY FIGURES
###Regulatory Divergence

#FG H1
type.legend = c("conserved","ambiguous",expression(paste(italic("trans"), "-only")), expression(paste(italic("cis"),"-only")),
                expression(paste(italic("cis + trans"), " (enhancing)")), expression(paste(italic("cis x trans"), " (compensatory)")),
                expression(paste(italic("cis-trans"), " (compensatory)")))

ggplot(fg_h1_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/fg_h1_ase.pdf", dpi = 300, width = 20, height = 8 )

#FG H2
type.legend = c("conserved","ambiguous",expression(paste(italic("trans"), "-only")), expression(paste(italic("cis"),"-only")),
                expression(paste(italic("cis + trans"), " (enhancing)")), expression(paste(italic("cis x trans"), " (compensatory)")),
                expression(paste(italic("cis-trans"), " (compensatory)")))

ggplot(fg_h2_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/fg_h2_ase.pdf", dpi = 300, width = 20, height = 8 )

#For shapes 21 - 25 , you can colour outline and fill inside separately as below

# ggplot(a, aes(logFC.sp, logFC.ase)) +
#   geom_point(alpha=0.6, show.legend=T, shape=21, size=2, color = "#484848", aes(fill=factor(type))) +
#   scale_fill_brewer(palette = "Dark2") +
#   geom_hline(yintercept=0, linetype=2) +
#   geom_vline(xintercept=0, linetype=2) +
#   ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
#   xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
#   background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
#   #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
#   theme(strip.background = element_rect(fill="grey90", color=NA)) +
#   theme(strip.background = element_blank(), strip.text=element_blank(),
#         axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

#FS
ggplot(fs_h1_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=3) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  ggtitle("FS H1")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/fs_h1_ase.pdf", dpi = 300, width = 20, height = 8 )


#FS H2
ggplot(fs_h2_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("FS H2")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/fs_h2_ase.pdf", dpi = 300, width = 20, height = 8 )

#MS
ms_h1_table = table(ms_h1_ase$type)

ggplot(ms_h1_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("MS H1")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/ms_h1_ase.pdf", dpi = 300, width = 20, height = 8 )


#WM
wm_h1_table = table(wm_h1_ase$type)
ggplot(wm_h1_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("WM H1")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/wm_h1_ase.pdf", dpi = 300, width = 20, height = 8 )

#####H2

#FS
fs_h2_table = table(fs_h2_ase$type)
ggplot(fs_h2_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

##MS
ms_h2_table = table(ms_h2_ase$type) 
ggplot(ms_h2_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("MS H2")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/ms_h2_ase.pdf", dpi = 300, width = 20, height = 8 )



#WM
wm_h2_table = table(wm_h2_ase$type)

ggplot(wm_h2_ase, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("WM H2")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/wm_h2_ase.pdf", dpi = 300, width = 20, height = 8 )

##MG
mg_h2_table = table(mg_h2_ase_3$type)

ggplot(mg_h2_ase_3, aes(logFC.sp, logFC.ase, color = type)) +
  geom_point(alpha=0.8, show.legend=T, size=4) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cre"),"/",italic("Clat"),")]", sep=''))) +
  background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
  #facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("MG H2")+
  ylim(-20, 20)+
  xlim(-20, 20)+
  theme(strip.background = element_blank(), strip.text=element_blank(),
        axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/Part 3/mg_h2_ase.pdf", dpi = 300, width = 20, height = 8 )
