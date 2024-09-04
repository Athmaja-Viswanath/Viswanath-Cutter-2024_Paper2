#####Geeting data for autosomes adn X-chromsome separately

##FG H1
head(fg_h1_ase_chr)
fg_h1_ase_autosome = fg_h1_ase_chr %>% filter(fg_h1_ase_chr$A_X == "Autosomes") #info for autosome
fg_h1_ase_autosome_table = table(fg_h1_ase_autosome$type)

fg_h1_ase_X = fg_h1_ase_chr %>% filter(fg_h1_ase_chr$A_X == "X") #info for x-chrosomosomes

fg_h1_ase_X$combined_regdiv = "Other"
fg_h1_ase_X[(fg_h1_ase_X$type == "cis-only"), "combined_regdiv"] = "cis-only"
fg_h1_ase_X[(fg_h1_ase_X$type == "trans-only"), "combined_regdiv"] = "trans-only"
fg_h1_ase_X[(fg_h1_ase_X$type == "cis-trans (compensatory)"), "combined_regdiv"] = "cis-trans (compensatory)"
fg_h1_ase_X[(fg_h1_ase_X$type == "conserved"), "combined_regdiv"] = "conserved"

fg_h1_ase_X_table = table(fg_h1_ase_X$combined_regdiv)

##FG H2
head(fg_h2_ase_chr)
##separating information for autosomes and x-chromosome
fg_h2_ase_autosome = fg_h2_ase_chr %>% filter(fg_h2_ase_chr$A_X == "Autosomes")
fg_h2_ase_X = fg_h2_ase_chr %>% filter(fg_h2_ase_chr$A_X == "X")

fg_h2_ase_autosome_table = table(fg_h2_ase_autosome$type)

fg_h2_ase_X$combined_regdiv = "Other"
fg_h2_ase_X[(fg_h2_ase_X$type == "cis-only"), "combined_regdiv"] = "cis-only"
fg_h2_ase_X[(fg_h2_ase_X$type == "trans-only"), "combined_regdiv"] = "trans-only"
fg_h2_ase_X[(fg_h2_ase_X$type == "cis-trans (compensatory)"), "combined_regdiv"] = "cis-trans (compensatory)"
fg_h2_ase_X[(fg_h2_ase_X$type == "conserved"), "combined_regdiv"] = "conserved"

fg_h2_ase_X_table = table(fg_h2_ase_X$combined_regdiv)

##FS H1
##separating information for autosomes and x-chromosome
fs_h1_ase_autosome = fs_h1_ase_chr %>% filter(fs_h1_ase_chr$A_X == "Autosomes")
fs_h1_ase_X = fs_h1_ase_chr %>% filter(fs_h1_ase_chr$A_X == "X")

fs_h1_ase_autosome_table = table(fs_h1_ase_autosome$type)

fs_h1_ase_X$combined_regdiv = "Other"
fs_h1_ase_X[(fs_h1_ase_X$type == "cis-only"), "combined_regdiv"] = "cis-only"
fs_h1_ase_X[(fs_h1_ase_X$type == "trans-only"), "combined_regdiv"] = "trans-only"
fs_h1_ase_X[(fs_h1_ase_X$type == "cis-trans (compensatory)"), "combined_regdiv"] = "cis-trans (compensatory)"
fs_h1_ase_X[(fs_h1_ase_X$type == "conserved"), "combined_regdiv"] = "conserved"

fs_h1_ase_X_table = table(fs_h1_ase_X$combined_regdiv)

##FS H2
##separating information for autosomes and x-chromosome
fs_h2_ase_autosome = fs_h2_ase_chr %>% filter(fs_h2_ase_chr$A_X == "Autosomes")
fs_h2_ase_X = fs_h2_ase_chr %>% filter(fs_h2_ase_chr$A_X == "X")
nrow(fs_h2_ase_autosome) #10465 genes detected on autosomes

fs_h2_ase_autosome_table = table(fs_h2_ase_autosome$type)

fs_h2_ase_X$combined_regdiv = "Other"
fs_h2_ase_X[(fs_h2_ase_X$type == "cis-only"), "combined_regdiv"] = "cis-only"
fs_h2_ase_X[(fs_h2_ase_X$type == "trans-only"), "combined_regdiv"] = "trans-only"
fs_h2_ase_X[(fs_h2_ase_X$type == "cis-trans (compensatory)"), "combined_regdiv"] = "cis-trans (compensatory)"
fs_h2_ase_X[(fs_h2_ase_X$type == "conserved"), "combined_regdiv"] = "conserved"

fs_h2_ase_X_table = table(fs_h2_ase_X$combined_regdiv)

##MS h1
ms_h1_ase_autosome = ms_h1_ase_chr %>% filter(ms_h1_ase_chr$A_X == "Autosomes")
nrow(ms_h1_ase_autosome) #10814 detectable on autosomes out of 11088 overall on autosomes


table(dds.ms.h1_X$reg.div)
nrow(dds.ms.h1_X) #2624 detectable out of possible 2665 X-linked genes

ms_h1_ase_autosome_table = table(ms_h1_ase_autosome$type)

ms_h1_ase_X_table = table(dds.ms.h1_X$reg.div)


##MS H2
ms_h2_ase_autosome = ms_h2_ase_chr %>% filter(ms_h2_ase_chr$A_X == "Autosomes")
nrow(ms_h2_ase_autosome) #10814 detectable on autosomes out of 11088 overall on autosomes

table(dds.ms.h2_X$reg.div)
nrow(dds.ms.h2_X) #2624 detectable out of possible 2665 X-linked genes

ms_h2_ase_autosome_table = table(ms_h2_ase_autosome$type)
ms_h2_ase_X_table = table(dds.ms.h2_X$reg.div)

##WM H1
wm_h1_ase_autosome = wm_h1_ase_chr %>% filter(wm_h1_ase_chr$A_X == "Autosomes")
nrow(wm_h1_ase_autosome) #10824 detectable on autosomes out of 11088 overall on autosomes

table(wm_h1_ase_autosome$type)

table(dds.wm.h1_X$reg.div)
head(dds.wm.h1_X)
nrow(dds.wm.h1_X) #2626 detectable out of possible 2665 X-linked genes

wm_h1_ase_autosome_table = table(wm_h1_ase_autosome$type)
wm_h1_ase_X_table = table(dds.wm.h1_X$reg.div)

##WM H2
wm_h2_ase_autosome = wm_h2_ase_chr %>% filter(wm_h2_ase_chr$A_X == "Autosomes")
nrow(wm_h2_ase_autosome) #10824 detectable on autosomes out of 11088 overall on autosomes
table(wm_h2_ase_autosome$type)

table(dds.wm.h2_X$reg.div)
nrow(dds.wm.h2_X) #2626 detectable out of possible 2665 X-linked genes

wm_h2_ase_autosome_table = table(wm_h2_ase_autosome$type)
wm_h2_ase_X_table = table(dds.wm.h2_X$reg.div)

#MG H2
mg_h2_ase_autosome = mg_h2_ase_chr %>% filter(mg_h2_ase_chr$A_X == "Autosomes")
nrow(mg_h2_ase_autosome) #10480 detectable on autosomes out of 11088 overall on autosomes
table(mg_h2_ase_autosome$type)

table(dds.mg.h2_X$reg.div)
head(dds.mg.h2_X)
nrow(dds.mg.h2_X) #2326 detectable out of possible 2665 X-linked genes

mg_h2_ase_autosome_table = table(mg_h2_ase_autosome$type)
mg_h2_ase_X_table = table(dds.mg.h2_X$reg.div)


#############################################################################################################
###Making cumulative data for autosomes
fg_h1_table.prop_A = cbind(as.data.frame(fg_h1_ase_autosome_table),as.data.frame(fg_h1_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(fg_h1_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("G", 7))) %>%
                           mutate(Hybrid = c(rep("H1", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))


fg_h2_table.prop_A = cbind(as.data.frame(fg_h2_ase_autosome_table),as.data.frame(fg_h2_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(fg_h2_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("G", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))  

fs_h1_table.prop_A = cbind(as.data.frame(fs_h1_ase_autosome_table),as.data.frame(fs_h1_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(fs_h1_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("S", 7))) %>%
                           mutate(Hybrid = c(rep("H1", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))

fs_h2_table.prop_A = cbind(as.data.frame(fs_h2_ase_autosome_table),as.data.frame(fs_h2_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(fs_h2_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("F", 7))) %>%
                           mutate(Tissue = c(rep("S", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))

ms_h1_table.prop_A = cbind(as.data.frame(ms_h1_ase_autosome_table),as.data.frame(ms_h1_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(ms_h1_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("M", 7))) %>%
                           mutate(Tissue = c(rep("S", 7))) %>%
                           mutate(Hybrid = c(rep("H1", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

ms_h2_table.prop_A = cbind(as.data.frame(ms_h2_ase_autosome_table),as.data.frame(ms_h2_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(ms_h2_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("M", 7))) %>%
                           mutate(Tissue = c(rep("S", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 
wm_h1_table.prop_A = cbind(as.data.frame(wm_h1_ase_autosome_table),as.data.frame(wm_h1_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(wm_h1_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("M", 7))) %>%
                           mutate(Tissue = c(rep("W", 7))) %>%
                           mutate(Hybrid = c(rep("H1", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

wm_h2_table.prop_A = cbind(as.data.frame(wm_h2_ase_autosome_table),as.data.frame(wm_h2_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(wm_h2_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("M", 7))) %>%
                           mutate(Tissue = c(rep("W", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid))

mg_h2_table.prop_A = cbind(as.data.frame(mg_h2_ase_autosome_table),as.data.frame(mg_h2_ase_autosome_table) %>%
                           mutate(Proportion = prop.table(mg_h2_ase_autosome_table)) %>%
                           mutate(Sex = c(rep("M", 7))) %>%
                           mutate(Tissue = c(rep("G", 7))) %>%
                           mutate(Hybrid = c(rep("H2", 7))) %>%
                           dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

regdiv_autosomal_counts = rbind(fg_h1_table.prop_A, fg_h2_table.prop_A, fs_h1_table.prop_A, fs_h2_table.prop_A,
                          ms_h1_table.prop_A, ms_h2_table.prop_A, mg_h2_table.prop_A, wm_h1_table.prop_A, wm_h2_table.prop_A)

regdiv_autosomal_counts = regdiv_autosomal_counts %>% rename(RegDiv = Var1, Num_of_genes = Freq)
#write.csv(regdiv_autosomal_counts, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/regdiv_autosomal_counts.csv")

############################################################################################################
###Making cumulative data for X-linked genes
fg_h1_table.prop_X = cbind(as.data.frame(fg_h1_ase_X_table),as.data.frame(fg_h1_ase_X_table) %>%
                             mutate(Proportion = prop.table(fg_h1_ase_X_table)) %>%
                             mutate(Sex = c(rep("F", 5))) %>%
                             mutate(Tissue = c(rep("G", 5))) %>%
                             mutate(Hybrid = c(rep("H1", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid))


fg_h2_table.prop_X = cbind(as.data.frame(fg_h2_ase_X_table),as.data.frame(fg_h2_ase_X_table) %>%
                             mutate(Proportion = prop.table(fg_h2_ase_X_table)) %>%
                             mutate(Sex = c(rep("F", 5))) %>%
                             mutate(Tissue = c(rep("G", 5))) %>%
                             mutate(Hybrid = c(rep("H2", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid))  

fs_h1_table.prop_X = cbind(as.data.frame(fs_h1_ase_X_table),as.data.frame(fs_h1_ase_X_table) %>%
                             mutate(Proportion = prop.table(fs_h1_ase_X_table)) %>%
                             mutate(Sex = c(rep("F", 5))) %>%
                             mutate(Tissue = c(rep("S", 5))) %>%
                             mutate(Hybrid = c(rep("H1", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid))

fs_h2_table.prop_X = cbind(as.data.frame(fs_h2_ase_X_table),as.data.frame(fs_h2_ase_X_table) %>%
                             mutate(Proportion = prop.table(fs_h2_ase_X_table)) %>%
                             mutate(Sex = c(rep("F", 5))) %>%
                             mutate(Tissue = c(rep("S", 5))) %>%
                             mutate(Hybrid = c(rep("H2", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid))

ms_h1_table.prop_X = cbind(as.data.frame(ms_h1_ase_X_table),as.data.frame(ms_h1_ase_X_table) %>%
                             mutate(Proportion = prop.table(ms_h1_ase_X_table)) %>%
                             mutate(Sex = c(rep("M", 5))) %>%
                             mutate(Tissue = c(rep("S", 5))) %>%
                             mutate(Hybrid = c(rep("H1", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

ms_h2_table.prop_X = cbind(as.data.frame(ms_h2_ase_X_table),as.data.frame(ms_h2_ase_X_table) %>%
                             mutate(Proportion = prop.table(ms_h2_ase_X_table)) %>%
                             mutate(Sex = c(rep("M", 5))) %>%
                             mutate(Tissue = c(rep("S", 5))) %>%
                             mutate(Hybrid = c(rep("H2", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid)) 
wm_h1_table.prop_X = cbind(as.data.frame(wm_h1_ase_X_table),as.data.frame(wm_h1_ase_X_table) %>%
                             mutate(Proportion = prop.table(wm_h1_ase_X_table)) %>%
                             mutate(Sex = c(rep("M", 5))) %>%
                             mutate(Tissue = c(rep("W", 5))) %>%
                             mutate(Hybrid = c(rep("H1", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

wm_h2_table.prop_X = cbind(as.data.frame(wm_h2_ase_X_table),as.data.frame(wm_h2_ase_X_table) %>%
                             mutate(Proportion = prop.table(wm_h2_ase_X_table)) %>%
                             mutate(Sex = c(rep("M", 5))) %>%
                             mutate(Tissue = c(rep("W", 5))) %>%
                             mutate(Hybrid = c(rep("H2", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid))

mg_h2_table.prop_X = cbind(as.data.frame(mg_h2_ase_X_table),as.data.frame(mg_h2_ase_X_table) %>%
                             mutate(Proportion = prop.table(mg_h2_ase_X_table)) %>%
                             mutate(Sex = c(rep("M", 5))) %>%
                             mutate(Tissue = c(rep("G", 5))) %>%
                             mutate(Hybrid = c(rep("H2", 5))) %>%
                             dplyr::select(Proportion, Sex, Tissue, Hybrid)) 

regdiv_X_counts = rbind(fg_h1_table.prop_X, fg_h2_table.prop_X, fs_h1_table.prop_X, fs_h2_table.prop_X,
                                ms_h1_table.prop_X, ms_h2_table.prop_X, mg_h2_table.prop_X, wm_h1_table.prop_X, wm_h2_table.prop_X)

regdiv_X_counts = regdiv_X_counts %>% rename(RegDiv = Var1, Num_of_genes = Freq)
#write.csv(regdiv_X_counts, file = "Allele Specific Expression/Hybrid Analysis/Figures_and_Results/regdiv_X_counts.csv")



