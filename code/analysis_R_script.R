#New analysis 2024

library(dplyr)
library(ggplot2)
library(vegan)
library(pheatmap)
library(reshape2)
library(stringr)
library(DESeq2)

metadata <- read.csv("Metadata.csv", row.names = 1)
alpha <- read.table("alpha_results.text", header=T) 
all.equal(metadata$X.SampleID, alpha$sample)

alpha$fish <- metadata$Fish.No.
alpha$group <- metadata$group
alpha$group2 <-metadata$group2
alpha$Population <- metadata$Population
alpha$treatment <- metadata$treatment
alpha$Treatment <- metadata$Treatment
alpha$elipse <-metadata$elipse

levels(alpha$Treatment) <- c("Control", "0.01 mg/L", "0.02 mg/L")

###alpha diveristy stats  
#tried alll the metrics, but sticking with shannon and chao

lm_chao <-lm(log(chao1)~Population * treatment, data = alpha)
shapiro.test(residuals(lm_chao))
anova(lm_chao)

lm_shannon <-lm(shannon_entropy~Population * treatment, data = alpha)
shapiro.test(residuals(lm_shannon))
anova(lm_shannon)


#graphs

ggplot(data=alpha,aes(x=Treatment, y=shannon_entropy, fill=Population))+ 
  geom_boxplot(col="black")+
  scale_fill_manual(values=c("white", "lightsteelblue2"))+
  theme_bw(base_size = 12)

ggplot(data=alpha,aes(x=Treatment, y=chao1, fill=Population))+ 
  geom_boxplot(col="black")+
  scale_fill_manual(values=c("white", "lightsteelblue2"))+
  theme_bw(base_size = 12)

####beta diversity

##nmds analysis of distance matrices (square format)

bc_dist=read.table("braycurtis_distance.tsv", header=T)
all.equal(rownames(metadata), rownames(bc_dist))

bc_nmds= metaMDS(bc_dist)
bc_nmds #stress=0.20 , data stored in bc_nmds$points
bc_points_NMDS <- as.data.frame(bc_nmds$points)
bc_points_NMDS$group <- metadata$group
bc_points_NMDS$population <- metadata$Population
bc_points_NMDS$treatment <- metadata$Treatment
bc_points_NMDS$elipse <- metadata$elipse

bc_points_NMDS$group <- factor(bc_points_NMDS$group, levels= c("0-naive",  "0.01-naive", "0.02-naive", "0-pre_ex", "0.01-pre_ex", "0.02-pre_ex"))
bc_points_NMDS$population <- factor(bc_points_NMDS$population, levels= c("Naive", "Pre_exposed"))

wuf_dist=read.table("W_unifrac_distance.tsv", header=T)
wuf_nmds= metaMDS(wuf_dist)
wuf_nmds #stress=0.14 , data stored in bc_nmds$points
wuf_points_NMDS <- as.data.frame(wuf_nmds$points)
wuf_points_NMDS$group <- metadata$group
wuf_points_NMDS$population <- metadata$Population
wuf_points_NMDS$treatment <- metadata$Treatment
wuf_points_NMDS$elipse <- metadata$elipse

wuf_points_NMDS$group <- factor(wuf_points_NMDS$group, levels= c("0-naive",  "0.01-naive", "0.02-naive", "0-pre_ex", "0.01-pre_ex", "0.02-pre_ex"))
wuf_points_NMDS$population <- factor(wuf_points_NMDS$population, levels= c("Naive", "Pre_exposed"))

ggplot(data=wuf_points_NMDS, aes(x=MDS1, y=MDS2, shape=population, col=group, size=population))+
  geom_point()+
  scale_colour_manual(values=c("dodgerblue3", "forestgreen", "firebrick2", "skyblue1", "palegreen", "lightpink"))+
  scale_shape_manual(values= c(16,17))+
  scale_size_manual(values= c(3,3))+
  theme_bw(base_size = 12)

ggplot(data=wuf_points_NMDS, aes(x=MDS1, y=MDS2, shape=population, col=group, size=population))+
  geom_point()+
  scale_colour_manual(values=c("dodgerblue3", "forestgreen", "firebrick2", "skyblue1", "palegreen", "lightpink"))+
  scale_shape_manual(values= c(16,17))+
  scale_size_manual(values= c(3,3))+
  stat_ellipse(data=wuf_points_NMDS, aes(group=elipse, linetype=elipse), size=0.5, type="norm", level=0.95)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  theme_bw(base_size = 12)+
  theme(legend.position = "bottom")

##BC (not as good stress, but nmds is only a visulaisation)

ggplot(data=bc_points_NMDS, aes(x=MDS1, y=MDS2, shape=population, fill=group, size=population))+
  geom_point()+
  scale_fill_manual(values=c("dodgerblue3", "forestgreen", "firebrick2", "skyblue1", "palegreen", "lightpink"))+
  scale_shape_manual(values= c(21,24))+
  scale_size_manual(values= c(4,3))+
  stat_ellipse(data=bc_points_NMDS, aes(group=elipse, linetype=elipse), size=0.5, type="norm", level=0.95)+
  scale_colour_manual(values=c("dodgerblue3", "forestgreen", "firebrick2", "skyblue3", "green", "red"))+
  scale_linetype_manual(values=c("solid", "longdash"))+
  theme_bw(base_size = 12)+
  theme(legend.position = "bottom")

ggplot(data=bc_points_NMDS, aes(x=MDS1, y=MDS2, shape=population, fill=group, size=population))+
  geom_point()+
  scale_fill_manual(values=c("dodgerblue3", "forestgreen", "firebrick2", "skyblue1", "palegreen", "lightpink"))+
  scale_shape_manual(values= c(21,24))+
  scale_size_manual(values= c(4,3))+
  scale_colour_manual(values=c("dodgerblue3", "forestgreen", "firebrick2", "skyblue3", "green", "red"))+
  scale_linetype_manual(values=c("solid", "longdash"))+
  theme_bw(base_size = 12)+
  theme(legend.position = "bottom")

#Beta diversity stats

adonis2(wuf_dist ~ metadata$Population * metadata$treatment, permutations = 99999)
adonis2(bc_dist ~ metadata$Population * metadata$treatment, permutations = 99999)

#beta disper
beta_df <- data.frame(Distance_to_centroid=beta$distances,Group=beta$group)
all.equal(rownames(metadata), rownames(beta_df))
beta_df$treatment <- metadata$Treatment
beta_df$population <- metadata$Population

ggplot(data=beta_df ,aes(x=treatment, y=Distance_to_centroid, fill= population))+ 
  geom_boxplot()+
  scale_fill_manual(values=c("white", "dodgerblue4"))+
  theme_bw(base_size = 12)+
  ylim(0,1)+
  ylab("Distance to centroid") + xlab(" ")+
  theme(legend.position = "none", legend.direction="horizontal")



#first format the ASV table:open the text file in excel, remove first line, add unique annotation, and save as a csv

####heatmaps

ASVs <- read.csv("asv_table.csv", row.names = 1, check.names = F)
metadata <- read.csv("metadata.csv", row.names=1)

levels(metadata$treatment) <- c("0.00", "0.01", "0.02")

#sort by most abundant OTUs, then remove the sum colum
ASVs <- cbind(ASVs, rowSums(ASVs))
ASVs <- ASVs[order(rowSums(-ASVs)),]
ASVs[,c("rowSums(ASVs)")] <- list(NULL)

logtop <- head(log2(ASVs+1), 100)

all.equal(colnames(logtop), row.names(metadata))
annot <- subset(metadata, select=c(Population,treatment))

my_colour = list(Population = c('Pre_exposed'= "chocolate1", 'Naive'="darkmagenta"), treatment = c('0.00' ="dodgerblue", '0.01' = "lightgreen", '0.02' = "firebrick"))
pheatmap(logtop, border_color = "grey40", 
         annotation_col =annot, 
         annotation_colors=my_colour, 
         cluster_rows=TRUE, show_rownames=TRUE, 
         show_colnames = FALSE, cluster_cols=TRUE, 
         fontsize_row = 5, fontsize_col = 6)


###Stacked barplots
#selecting the most abundant 20 ASVs
Bar_OTUs <- read.csv("asv_table.csv", row.names = 1, check.names = F)

Bar_OTUs <- cbind(Bar_OTUs, rowSums(Bar_OTUs))
Bar_OTUs <- Bar_OTUs[order(rowSums(-Bar_OTUs)),]
Bar_OTUs[,c("rowSums(Bar_OTUs)")] <- list(NULL)
Bar20 <- head(Bar_OTUs, 20)

tBar20 <- t(Bar20) #transposing data
mtop20 = melt(tBar20, id.vars = c("ASV"))

mtop20$Population <- metadata$Population
mtop20$Treatment <- metadata$treatment
mtop20$Group <- metadata$group
mtop20$Rep <- metadata$rep

colours = c("#8B0000",  "#008000", "#FF0000","#90EE90",  "#008080" , "#F6AE2D", "#DC134C", "#7B68EE", "#DDA0DD", "#F9ECCC", "#000000",  "#A54657",  "#582630", "#F7EE7F", "#191970", "#679289",  "#33658A","#F1A66A", "#86BBD8", "#4DAA57")

ggplot(mtop20, aes(x = Rep, fill = Var2, y = value)) + 
  geom_bar(stat = "identity", colour = "black")+ facet_grid(rows= vars(Population), cols=vars(Treatment), scales = "free")+
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        axis.text.y = element_text(size = 10)) + 
  labs(x = "", y = "Abundance (counts)", fill = "OTU")+
  scale_fill_manual(values = colours)

ggplot(mtop20, aes(x = Rep, y = Var2)) + 
  geom_point(aes(size= value, fill=Var2), alpha=0.75, shape=21,)+
  facet_grid(cols=vars(Treatment), rows=vars(Population), scales = "free", space= 'free')+
  scale_size_continuous(limits = c(1, 9000), range = c(1,10), breaks = c(20,100,500,2500)) +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, ), 
        axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), 
        axis.text.y = element_text(size = 8)) + 
  labs(x = "", y = "", fill = "", size= "Counts")

  scale_fill_manual(values = palette30)



#20 most abundant genera
#from bar plots, level 6, removed additional info, made new taxonomic names

Bar_gen <- read.csv("level_6_top20.csv", row.names = 1, check.names = F)

Bar_gen <- cbind(Bar_gen, rowSums(Bar_gen))
Bar_gen <- Bar_gen[order(rowSums(-Bar_gen)),]
Bar_gen[,c("rowSums(Bar_gen)")] <- list(NULL)
Bar20gen <- head(Bar_gen, 20)

tBar20gen <- t(Bar20gen) #transposing data
mtop20gen = melt(tBar20gen, id.vars = c("Taxonomy"))

mtop20gen$Population <- metadata$Population
mtop20gen$Treatment <- metadata$treatment
mtop20gen$Group <- metadata$group
mtop20gen$Rep <- metadata$rep

mtop20gen$Var2 <- reorder(mtop20gen$Var2, mtop20gen$value)

colours = c("#556B2F", "#C1FFC1", "#679289", "#008080" , "#CCE5FF", "#DC134C", "#7B68EE", "#DDA0DD", "#F9ECCC", "#000000",  "#A54657",  "#808080", "#191970", "#F7EE7F", "#8B0000",  "#FF0000", "#F6AE2D", "#86BBD8","#33658A","#4DAA57")

#facet labels
Treatment.labs <- c("0 mg/L", "0.01 mg/L", "0.02 mg/L")
names(Treatment.labs) <- c("0", "0.01", "0.02")

Pop.labs <- c("Naive population", "Pre-exposed population")
names(Pop.labs) <- c("Naive", "Pre_exposed")

ggplot(mtop20gen, aes(x = Rep, fill = Var2, y = value)) + 
  geom_bar(stat = "identity", colour = "black")+ facet_grid(Population~Treatment, scales = "free", labeller = labeller(Population=Pop.labs, Treatment=Treatment.labs))+
  theme(axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size=12, face= "bold"),
        panel.background = element_rect(fill= "white", colour = "black"),
        panel.grid.major = element_line(colour="grey")) + 
  labs(x = "", y = "Abundance (counts)", fill = "Genera")+
  scale_fill_manual(values = colours)

#20 most abundant genera
#from bar plots, level 6, removed additional info, made new taxonomic names> level_6.csv
#to make names neater, just did the top20 in excel 

Bar_gen <- read.csv("level_6_top20.csv", row.names = 1, check.names = F)

tBar20gen <- t(Bar_gen) #transposing data
mtop20gen = melt(tBar20gen, id.vars = c("Genus"))

mtop20gen$Population <- metadata$Population
mtop20gen$Treatment <- metadata$treatment
mtop20gen$Group <- metadata$group
mtop20gen$Rep <- metadata$rep

#mtop20gen$Var2 <- reorder(mtop20gen$Var2, mtop20gen$value)

colours = c("#4DAA57", "#C1FFC1", "yellow","#FF0000", "#7B68EE","#CCE5FF", "#DC134C", "#008080" ,"#679289", "#DDA0DD", "#F9ECCC", "#000000",  "#A54657", "#F7EE7F", "#86BBD8",  "#808080", "#191970","#33658A", "#8B0000","#556B2F" )

#colours = c("#F7EE7F" ,"#191970", "purple", "#FF0000", "#DC134C","#C1FFC1", "#4DAA57", "#CCE5FF", "#8B0000","#808080", "pink1" , "#556B2F", "#86BBD8", "#7B68EE","#F9ECCC", "#000000", "#DDA0DD","#F6AE2D", "limegreen", "#33658A")
#colours = c("#4DAA57","#CCE5FF", "yellow", "#FF0000", "#DC134C","#C1FFC1","#F9ECCC", "#4DAA57", "#191970", "#8B0000","#808080", "pink1" , "#556B2F", "#86BBD8", "#7B68EE", "#F7EE7F" ,"#000000", "#DDA0DD","#F6AE2D", "limegreen", "#33658A")

#facet labels
Treatment.labs <- c("0 mg/L", "0.01 mg/L", "0.02 mg/L")
names(Treatment.labs) <- c("0", "0.01", "0.02")

Pop.labs <- c("Naive population", "Pre-exposed population")
names(Pop.labs) <- c("Naive", "Pre_exposed")

ggplot(mtop20gen, aes(x = Rep, fill = Var2, y = value)) + 
  geom_bar(stat = "identity", colour = "black")+ facet_grid(Population~Treatment, scales = "free", labeller = labeller(Population=Pop.labs, Treatment=Treatment.labs))+
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size=12, face= "bold"),
        panel.background = element_rect(fill= "white", colour = "black"),
        panel.grid.major = element_line(colour="grey")) + 
  labs(x = "", y = "Abundance (counts)", fill = "Genera")+
  scale_fill_manual(values = colours)

ggplot(mtop20gen, aes(x = Rep, y = Var2)) + 
  geom_point(aes(size= value, fill=Var2), alpha=0.75, shape=21,)+
  facet_grid(cols=vars(Treatment), rows=vars(Population), scales = "free", space= 'free')+
  scale_size_continuous(limits = c(1, 9000), range = c(1,10), breaks = c(20,100,500,2500)) +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, ), 
        axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), 
        axis.text.y = element_text(size = 8)) + 
  labs(x = "", y = "", fill = "", size= "Counts")



#DESeq2

#Different this time--- analysis to match the RNA-seq and RRBS--- need to split datasets 
#1 Testing effect of pre-exposure: naive-0 v prex-0
#2 Testing the effect of copper exposure on both populations


#will need to subset the datset to do differnet comparsions- adding a unique name that can be used for subsetting

pre_countData <- t(read.csv("asv_table.csv", row.names = "ASV", check.names = F))
row.names(pre_countData)

colData <- read.csv("metadata.csv", row.names=1)
colData$uniquename <- paste(colData$group, colData$rep)

all.equal(row.names(pre_countData),row.names(colData))

row.names(colData) <- colData$uniquename
row.names(pre_countData) <- colData$uniquename
row.names(pre_countData)

#now subsetting the Countdata and Coldata, and transpose countData

colData_0cop <- colData[grep(pattern = "0-", rownames(colData)),]
countData_0cop <- t(pre_countData[grep(pattern = "0-", rownames(pre_countData)),])

colData_prex <- colData[grep(pattern = "pre_ex", rownames(colData)),]
countData_prex <- t(pre_countData[grep(pattern = "pre_ex", rownames(pre_countData)),])

colData_naive <- colData[grep(pattern = "naive", rownames(colData)),]
countData_naive <- t(pre_countData[grep(pattern = "naive", rownames(pre_countData)),])

#1.prex effect

countData_0cop <- (countData_0cop+1)
dim(countData_0cop)
colnames(countData_0cop)
row.names(countData_0cop)

dim(colData_0cop)
colnames(colData_0cop)
row.names(colData_0cop)

all(rownames(colData_0cop) %in% colnames(countData_0cop))

DES <- DESeqDataSetFromMatrix(countData = countData_0cop, colData = colData_0cop, design = ~ Population)
DES <- DES[rowSums(counts(DES) > 3) >=8 ] 

dim(DES)   #36ASVs remaining 

design(DES) <- formula(~ Population)
DES <- DESeq(DES) 
resultsNames(DES)

summary(results(DES, alpha=0.05, contrast= c("Population", "Pre_exposed", "Naive")), alpha=0.05)
pop_res <- results(DES, alpha=0.05, contrast= c("Population", "Pre_exposed", "Naive"))
sig_pop_otus <- rownames(subset(pop_res, results(DES, alpha=0.05, contrast= c("Population", "Pre_exposed", "Naive"))$padj <0.05)) #sig results
write.csv(pop_res, "Prex_v_naive_results.csv")

#2. cop effect in naive

countData_naive <- (countData_naive+1)
dim(countData_naive)
colnames(countData_naive)
row.names(countData_naive)

dim(colData_naive)
colnames(colData_naive)
row.names(colData_naive)

all(rownames(colData_naive) %in% colnames(countData_naive))

DES_naive <- DESeqDataSetFromMatrix(countData = countData_naive, colData = colData_naive, design = ~ Treatment)
DES_naive <- DES_naive[rowSums(counts(DES_naive) > 3) >=8 ] 

dim(DES_naive)   #78ASVs remaining 

design(DES_naive) <- formula(~ Treatment)
DES_naive <- DESeq(DES_naive) 
resultsNames(DES_naive)

summary(results(DES_naive, alpha=0.05, contrast= c("Treatment", "b", "a")), alpha=0.05)
summary(results(DES_naive, alpha=0.05, contrast= c("Treatment", "c", "a")), alpha=0.05)

naive_res <- results(DES_naive, alpha=0.05, contrast= c("Treatment", "b", "a"))
naive_res2 <- results(DES_naive, alpha=0.05, contrast= c("Treatment", "c", "a"))
write.csv(naive_res, "Naive_0.01v0_results.csv")
write.csv(naive_res2, "Naive_0.02v0_results.csv")

#3 cop effect in prex

countData_prex <- (countData_prex+1)
dim(countData_prex)
colnames(countData_prex)
row.names(countData_prex)

dim(colData_prex)
colnames(colData_prex)
row.names(colData_prex)

all(rownames(colData_prex) %in% colnames(countData_prex))

DES_prex <- DESeqDataSetFromMatrix(countData = countData_prex, colData = colData_prex, design = ~ Treatment)
DES_prex <- DES_prex[rowSums(counts(DES_prex) > 3) >=10 ] 

dim(DES_prex)   #95 ASVs remaining 

design(DES_prex) <- formula(~ Treatment)
DES_prex <- DESeq(DES_prex) 
resultsNames(DES_prex)

summary(results(DES_prex, alpha=0.05, contrast= c("Treatment", "b", "a")), alpha=0.05)
summary(results(DES_prex, alpha=0.05, contrast= c("Treatment", "c", "a")), alpha=0.05)

prex_res <- results(DES_prex, alpha=0.05, contrast= c("Treatment", "b", "a"))
prex_res2 <- results(DES_prex, alpha=0.05, contrast= c("Treatment", "c", "a"))
write.csv(prex_res, "Prex_0.01v0_results.csv")
write.csv(prex_res2, "Prex_0.02v0_results.csv")
