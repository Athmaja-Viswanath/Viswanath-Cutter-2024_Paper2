###COMAPRING H1 AND H2 

###INHERITANCE CATEGORY BETWEEN H1 AND H2

###########
###FG######
###########

###VARIABLES SAVED FROM 2.cClassifying Inheritance.R
head(rownames_to_column(dds.fg.h1))
dds.fg.h1

`%notin%` <- Negate(`%in%`)

##chanign rownames as column
dds.fg.h1.1 = rownames_to_column(dds.fg.h1) 
dds.fg.h2.1 = rownames_to_column(dds.fg.h2)

##ordering each dataframe
dds.fg.h1.1 = dds.fg.h1.1[order(dds.fg.h1.1$rowname), ]
dds.fg.h2.1 = dds.fg.h2.1[order(dds.fg.h2.1$rowname), ] %>%
  rename(fg.cre.clat_2 = fg.cre.clat, rowname.h2 = rowname, inheritance.fg.h2 = inheritance.fg)
head(dds.fg.h2.1)

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


ggplot(data = fg.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Part3/fg.h1h2.inht.pdf", dpi = 300, width = 18, height = 8 )


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

ggplot(data = fs.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Part3/fs.h1h2.inht.pdf", dpi = 300, width = 18, height = 8 )

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

nrow(dds.ms.h1.1 %>% filter(rowname %in% dds.ms.h2.1$rowname))

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

ggplot(data = ms.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Part3/ms.h1h2.inht.pdf", dpi = 300, width = 18, height = 8 )

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

ggplot(data = wm.h1h2.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")


ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Part3/wm.h1h2.inht.pdf", dpi = 300, width = 18, height = 8 )


#######################################################################################################
####VISAULIZING BUBBLEPLOTS FOR REGULATORY DIVERGECNE BETWEEN H1 & H2##################################
#######################################################################################################

###########
###FG######
###########

###VARIABLES SAVED FROM 3.Allelespecific Expression.R
fg_h1_ase
fg_h2_ase

`%notin%` <- Negate(`%in%`)

##changing rownames to column
fg_h1_ase.1 = rownames_to_column(fg_h1_ase) 
fg_h2_ase.1 = rownames_to_column(fg_h2_ase)

nrow(fg_h1_ase.1 %>% filter(rowname %in% fg_h2_ase.1$rowname))

#ordering dataframes
fg_h1_ase.1 = fg_h1_ase.1[order(fg_h1_ase.1$rowname), ]
fg_h2_ase.1 = fg_h2_ase.1[order(fg_h2_ase.1$rowname), ] 
head(fg_h2_ase.1)

##joining dataframes and renaming 
fg.h1h2.ase = as.data.frame(cbind(fg_h1_ase.1$rowname, fg_h1_ase.1$type,fg_h2_ase.1$rowname,fg_h2_ase.1$type)) %>%
  rename(H1_name = V1, H1_regdiv = V2, H2_name = V3, H2_regdiv = V4)
head(as.data.frame(fg.h1h2.ase))
head((fg.h1h2.ase))

#getting intersection between gene categories
fg.h1h2.ase.intersection = as.data.frame(table(fg.h1h2.ase$H1_regdiv, fg.h1h2.ase$H2_regdiv)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq)
fg.h1h2.ase.intersection = fg.h1h2.ase.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)


##Bubbleplot
ggplot(data = fg.h1h2.ase.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Part3/fg.h1h2_regdiv.pdf", dpi = 300, width = 18, height = 8 )


###########
###FS######
###########

###VARIABLES SAVED FROM 3.Allelespecific Expression.R
fs_h1_ase
fs_h2_ase

#Changing rownames to column
fs_h1_ase.1 = rownames_to_column(fs_h1_ase) 
fs_h2_ase.1 = rownames_to_column(fs_h2_ase)

nrow(fs_h1_ase.1 %>% filter(rowname %in% fs_h2_ase.1$rowname))

##ordering both dataframes
fs_h1_ase.1 = fs_h1_ase.1[order(fs_h1_ase.1$rowname), ]
fs_h2_ase.1 = fs_h2_ase.1[order(fs_h2_ase.1$rowname), ] 
head(fs_h2_ase.1)

##combining both dataframes and renaming columns
fs.h1h2.ase = as.data.frame(cbind(fs_h1_ase.1$rowname, fs_h1_ase.1$type,fs_h2_ase.1$rowname,fs_h2_ase.1$type)) %>%
  rename(H1_name = V1, H1_regdiv = V2, H2_name = V3, H2_regdiv = V4)
head(as.data.frame(fs.h1h2.ase))
head((fs.h1h2.ase))

###getting intersection frequencies for gene categories
head(table(fs.h1h2.ase$H1_regdiv, fs.h1h2.ase$H2_regdiv))
fs.h1h2.ase.intersection = as.data.frame(table(fs.h1h2.ase$H1_regdiv, fs.h1h2.ase$H2_regdiv)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq)
fs.h1h2.ase.intersection = fs.h1h2.ase.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)

##Bubble plot
ggplot(data = fs.h1h2.ase.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Part3/fs.h1h2_regdiv.pdf", dpi = 300, width = 18, height = 8 )


###########
###MS######
###########

###VARIABLES SAVED FROM 3.Allelespecific Expression.R
#with updated categories for X-linked genes, got from 3.Allele specific expression.R
ms_h1_ase_3 #
ms_h2_ase_3

##ordeing the dataframe
ms_h1_ase.1 = ms_h1_ase_3[order(ms_h1_ase_3$rowname), ]
ms_h2_ase.1 = ms_h2_ase_3[order(ms_h2_ase_3$rowname), ] 
head(ms_h1_ase.1)

####Joining H1 and H2 information and renaming columns
ms.h1h2.ase = as.data.frame(cbind(ms_h1_ase.1$rowname, ms_h1_ase.1$reg.div, ms_h2_ase.1$rowname,ms_h2_ase.1$reg.div)) %>%
  rename(H1_name = V1, H1_regdiv = V2, H2_name = V3, H2_regdiv = V4)
head(as.data.frame(ms.h1h2.ase))
#View((ms.h1h2.ase))

#####Finding intersection in categories between H1 and H2
head(table(ms.h1h2.ase$H1_regdiv, ms.h1h2.ase$H2_regdiv))
ms.h1h2.ase.intersection = as.data.frame(table(ms.h1h2.ase$H1_regdiv, ms.h1h2.ase$H2_regdiv)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq)
ms.h1h2.ase.intersection = ms.h1h2.ase.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)

###Plotting bubble plot for all genes other than ambiguous genes or catgories with 0 intersection
ggplot(data = ms.h1h2.ase.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Basic Results/Part 3/ms.h1h2_regdiv.pdf", dpi = 300, width = 18, height = 8 )


###########
###WM######
###########

###VARIABLES SAVED FROM 3.Allelespecific Expression.R
##with updated categories for X-linked genes, got from 3.Allele specific expression.R
wm_h1_ase_3
wm_h2_ase_3

##ORDERING DATAFRAMES 
wm_h1_ase.1 = wm_h1_ase_3[order(wm_h1_ase_3$rowname), ]
wm_h2_ase.1 = wm_h2_ase_3[order(wm_h2_ase_3$rowname), ] 
head(wm_h2_ase.1)

##combining and renaming columns
wm.h1h2.ase = as.data.frame(cbind(wm_h1_ase.1$rowname, wm_h1_ase.1$reg.div,wm_h2_ase.1$rowname,wm_h2_ase.1$reg.div)) %>%
  rename(H1_name = V1, H1_regdiv = V2, H2_name = V3, H2_regdiv = V4)
head(as.data.frame(wm.h1h2.ase))
#View((wm.h1h2.ase))

##getting intersection between gene categories
head(table(wm.h1h2.ase$H1_regdiv, wm.h1h2.ase$H2_regdiv))
wm.h1h2.ase.intersection = as.data.frame(table(wm.h1h2.ase$H1_regdiv, wm.h1h2.ase$H2_regdiv)) %>%
  rename(H1 = Var1, H2=Var2, Intersection = Freq)
wm.h1h2.ase.intersection = wm.h1h2.ase.intersection %>% filter(H1 != "ambiguous" & H2 != "ambiguous" & Intersection > 0)

####Bubble Plot
ggplot(data = wm.h1h2.ase.intersection, aes(x=H1, y=H2, size = Intersection))+
  geom_count()+
  geom_text(aes(label = Intersection), size = 5, colour = "red", nudge_x = 0.1, nudge_y = 0.1)+
  scale_size(range = c(2, 20), name="Number of Genes")

ggsave(filename = "Allele Specific Expression/Hybrid Analysis/Basic Results/Part 3/wm.h1h2_regdiv.pdf", dpi = 300, width = 18, height = 8 )


