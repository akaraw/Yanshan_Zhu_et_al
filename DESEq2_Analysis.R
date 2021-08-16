####################################### ~ Difference between age groups with or without treatment effect####
####################################### ~ Loading libraries ~ #####
library("DESeq2")
library(dplyr)
library("RColorBrewer")
library("gplots")
library( "genefilter" )
library("pheatmap")
library(tools)
library("ggplot2")
library(magrittr)
library("biomaRt")
library(apeglm)
library("genefilter")
library(tximport)
library(pcaExplorer)
library(GenomicFeatures)
library(goseq)
library(GO.db)
library(fgsea)
library('EnhancedVolcano')
library(PCAtools)
library(dplyr)
####################################### ~ DECLARE YOUR VARIABLES HERE ~ #####
myspecies <- "H.sapiens"

mygenes <- c("MYD88", "TMPRSS2", "TRAIP", "LRRFIP1", "PIK3IP1", "TRAM1", "DDX58", "TLR3") #, "TLR7") 
#"TICAM1", "IRF3", "IRAK4", "IRAK2", "TRAF6", "IFIH1", "DHX58",
#"MAVS")
directory <- "D:/17.SARSC-V2_experiment/"
setwd("D:/17.SARSC-V2_experiment/R-project/SCVO2_KVSA/")
resultsdir <- paste0(directory, "results")

baselevel <-  "Adult" 
baselevel2 <- "Kid"
ageclass <- "mocks"
myresdir <- "analysis1"
testgroup <- "/salmon_res/"

dir.create(paste0(resultsdir))
resultsdir
dir.create(paste0( resultsdir, testgroup))

resultsdir <- paste0(resultsdir,testgroup)
resultsdir
dir.create(paste0(resultsdir, myresdir))
resultsdir <- paste0(resultsdir, myresdir)
resultsdir

####################################### ~ Count tables and metadata ###############################################
sample_table <- read.table(paste0("D:/17.SARSC-V2_experiment/YZ_experiment.txt"), sep = "\t", header=TRUE)
length(sample_table)
sample_table

(pull(sample_table, `ID`))

sampleFiles <- paste0("../../SCOV2_Study_19122020/salmon_quant/",pull(sample_table, `ID`),"/", "quant.sf")
sampleFiles
class(sampleFiles)

all(file.exists(sampleFiles))

names(sampleFiles) <- paste0(pull(sample_table, 'ID'), '_', 
                             pull(sample_table, 'Experiment_number'))
names(sampleFiles)

fac <- 'Age'
fac1 <- 'Treatment'
fac2 <- 'Experiment_number'

sampleCondition<-pull(sample_table, fac)
sampleCondition
sampleCondition1 <- pull(sample_table,fac1)
sampleCondition1
sampleCondition2 <- pull(sample_table, fac2)
sampleCondition2
length(sampleFiles)
sampleCondition3 <- pull(sample_table, 'Gender')
sampleCondition3
sampleTable <- data.frame(sampleName = names(sampleFiles), fileName = sampleFiles, 
                          gender = as.factor(sampleCondition3),
                          age = as.factor(sampleCondition), condition = as.factor(sampleCondition1),
                          pairID = as.factor(sampleCondition2))
head(sampleTable)
str(sampleTable)

group = as.factor(sampleTable$condition)
group


tx2gene <- read.table("D:/17.SARSC-V2_experiment/tx2genesymbol.txt", header=FALSE, sep="\t")

tx2gene <- dplyr::select(tx2gene, -V3)

head(tx2gene)

count_data <- tximport(files = sampleFiles,
                       type = "salmon",
                       tx2gene = tx2gene,
                       ignoreTxVersion = FALSE)

####################################### ~ Import data into DESeq2 and process#################################
head(sampleTable)
(sampleTable$age)

####Design >>>>>
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                #design = ~ condition + age)
                                design = ~ gender + age + condition + age:condition)
                                
design(dds) #Confirm the design used in the analysis
baselevel <- "Adult" #Setting up reference levels
dds$age <- relevel(dds$age, ref = baselevel)
dds$condition <- relevel(dds$condition, ref = "Mock")

keep <- rowSums(counts(dds)) > 10 #FIltering out low abundance reads
ddsKEEP <- dds[keep,]

ddsKEEP <- DESeq(ddsKEEP) #DESeq analysis
colData(ddsKEEP)
resultsNames(ddsKEEP)
####################################### ~ First round of results and other stats#################################################
ageclass <- myresdir
nonnormalized_counts <- counts(dds)
resultsdir
write.table(nonnormalized_counts, file=paste(resultsdir, "/", ageclass, 
                                             "sarscov2_non_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

normalized_counts <- counts(ddsKEEP, normalized=TRUE)
head(normalized_counts)
t <- as.data.frame(normalized_counts)
head(t)
colnames(t)
ageclass
resultsdir
class(normalized_counts)
write.table(t, file=paste(resultsdir, "/", ageclass, "_hsapiens_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

####################################### ~ Different contrast of the results from DEseq2 analysis####################################
#What is the difference between Kid and Age without treatment?
resultsNames(ddsKEEP)
resAge <- results(ddsKEEP, tidy = F, name = "age_Kid_vs_Adult", alpha = 0.05)
res <- resAge
DGE.df <- as.data.frame(res)
DGE.df
getwd()
write.csv(DGE.df, file = paste0(resultsdir,"/differential_kid_vs_adult_without_treatment.csv"))
res$log2FoldChange
sum(res$padj<0.05, na.rm = T)
res[which(res$padj < 0.05),]
table(is.na(res$padj))
resultsNames(ddsKEEP)

resSig= res[which(res$padj<0.05),] #Significant genes
testgroup <- ageclass
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_kid_vs_adult_without_virus_DEGS_hsapiens.csv", sep = "") )
resultsdir
(resSig_up= resSig[which(resSig$log2FoldChange > 0.58),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_kid_vs_adult_without_virus_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -0.58),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_kid_vs_adult_without_virus_DOWN_DEGS_hsapiens.csv", sep = ""))
rm(res)

#With treatment what is the difference between the kid and the adult.

resCon <- results(ddsKEEP, tidy = F, list(c("age_Kid_vs_Adult", "ageKid.conditionVirus")), alpha = 0.05)
head(resCon)
res <- resCon
res$log2FoldChange
sum(res$padj<0.05, na.rm = T)
table(is.na(res$padj))
resultsNames(ddsKEEP)

resSig= res[which(res$padj<0.05),] #Significant genes
testgroup <- ageclass
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_kid_vs_adult_with_viruseffect_DEGS_hsapiens.csv", sep = "") )
resultsdir
(resSig_up= resSig[which(resSig$log2FoldChange > 0.58),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_kid_vs_adult_with_viruseffect_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -0.58),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_kid_vs_adult_with_viruseffect_DOWN_DEGS_hsapiens.csv", sep = ""))
rm(res)

#The effect of treatment (i.e. Virus) in Kids only
resKid <- results(ddsKEEP, tidy = F, list(c("condition_Virus_vs_Mock", "ageKid.conditionVirus")), alpha = 0.05)
head(resKid)
res <- resKid
#res$log2FoldChange
sum(res$padj<0.05, na.rm = T)
table(is.na(res$padj))
resultsNames(ddsKEEP)

resSig= res[which(res$padj<0.05),] #Significant genes
testgroup <- ageclass
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_kid_effect_of_virus_DEGS_hsapiens.csv", sep = "") )
resultsdir
(resSig_up= resSig[which(resSig$log2FoldChange > 0.58),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_kid_effect_of_virus_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -0.58),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_kid_effect_of_virus_DOWN_DEGS_hsapiens.csv", sep = ""))
rm(res)

####Effect of the virus in adults only
design(dds) #Confirm the design used in the analysis
dds$age <- relevel(dds$age, ref = 'Adult')
dds$condition <- relevel(dds$condition, ref = "Mock")

keep <- rowSums(counts(dds)) > 10 #FIltering out low abundance reads
ddsKEEP <- dds[keep,]

ddsKEEP <- DESeq(ddsKEEP) #DESeq analysis
colData(ddsKEEP)
resultsNames(ddsKEEP)

resAdult = results(ddsKEEP, contrast=c("condition","Virus","Mock"))
dds$age
dds$condition

res <- resAdult
#res$log2FoldChange
sum(res$padj<0.05, na.rm = T)
table(is.na(res$padj))
resultsNames(ddsKEEP)

resSig= res[which(res$padj<0.05),] #Significant genes
testgroup <- ageclass
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_Adult_effect_of_virus_DEGS_hsapiens.csv", sep = "") )
resultsdir
(resSig_up= resSig[which(resSig$log2FoldChange > 0.58),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_Adult_effect_of_virus_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -0.58),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_Adult_effect_of_virus_DOWN_DEGS_hsapiens.csv", sep = ""))
rm(res)

design(dds)
resultsNames(ddsKEEP)
resInt = results(ddsKEEP, name = "ageKid.conditionVirus")

res <- resInt
#res$log2FoldChange
sum(res$padj<0.05, na.rm = T)
table(is.na(res$padj))

resSig= res[which(res$padj<0.05),] #Significant genes
testgroup <- ageclass
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_Ineraction_DEGS_hsapiens.csv", sep = "") )
resultsdir
(resSig_up= resSig[which(resSig$log2FoldChange > 0.58),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_Interaction_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -0.58),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_Interaction_DOWN_DEGS_hsapiens.csv", sep = ""))
rm(res)

####################################### ~ MA plots ~ ####
plotMA(resAge,ylim=c(-2,2))

plotDispEsts( ddsKEEP, ylim = c(1e-6, 1e1) )

hist( resAge$pvalue, breaks=20, col="green" )

resultsNames(ddsKEEP)

resLFC <- lfcShrink(ddsKEEP, coef = 3, type="apeglm")
resLFC
plotMA(resLFC,ylim=c(-5,5))
abline(h=c(-0.5,0.5), col="red")

# create bins using the quantile function
qs <- c( 0, quantile( resAge$baseMean[resAge$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( resAge$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of Â£pÂ£ values less than .01 for each bin
ratios <- tapply( resAge$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

####################################### ~ PCA plots ~ ####
source('mypcaplots.R')
mypcaplot(ddsKEEP)
####################################### ~ Differential gene expression of given genes ~ ######################################
source('mytopgenes.R')
mytopgenes(ddsKEEP)

####################################### ~ Another method of PCA plot ~ ####
graphics.off()
source('myotherpca.R')
myotherpca(ddsKEEP, testgroup)

####################################### ~ GO Enrichment analysis for NCBI based annotations ~ #################################
head(tx2gene)
count_data <- tximport(files = sampleFiles,
                       type = "salmon",
                       tx2gene = tx2gene,
                       ignoreTxVersion = FALSE)
head(count_data)
#import data into DESeq2
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                design = ~ age + condition + condition:age)
dds
baselevel
dds$age <- relevel(dds$age, ref = baselevel)
dds$age
baselevel2 <- "Mock"
dds$condition <- relevel(dds$condition, ref = baselevel2)
dds$condition
rowsum.threshold <- 5
fdr.threshold <- 0.05
rs <- rowSums(counts(ddsKEEP))
rs
ddsG <- ddsKEEP[rs > rowsum.threshold,]
ddsG<-DESeq(ddsG)
resultsNames(ddsG)

source('mygoscript.R')
source('mypathview.R')
ensembl<- useMart("ensembl")
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
sp <- useDataset("hsapiens_gene_ensembl",mart = ensembl)
GTF <- "../../gencode.v36.annotation.gtf"
txdb = makeTxDbFromGFF(GTF, format = "gtf", )
txdb.allcolumns <- transcripts(txdb)
txdb.allcolumns
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
graphics.off()

ddsG <- DESeq(dds)

resG <- results(ddsG, tidy = F, name = "age_Kid_vs_Adult", alpha = 0.05)
resG$log2FoldChange
testgroup <- "age_Kid_vs_Adult"
testgroup
my_go_script(resG, testgroup, 0.05)
mykegg(resG, testgroup)

resG <- results(ddsG, tidy = F, name = "condition_Virus_vs_Mock", alpha = 0.05)
testgroup <- "adults_with_virus_"
my_go_script(resG, testgroup, 0.05)
mykegg(resG, testgroup)

resG <- results(ddsG, tidy = F, list(c("age_Kid_vs_Adult", "ageKid.conditionVirus"))
                , alpha = 0.05)
testgroup <- "kids_vs_adults_with_virus_"
my_go_script(resG, testgroup, 0.05)
mykegg(resG, testgroup)

resG <- results(ddsG, tidy = F, 
                  list(c("condition_Virus_vs_Mock", "ageKid.conditionVirus")), alpha = 0.05)
resG <- results(ddsG, tidy = F, 
                list(c("condition_Virus_vs_Mock", "ageKid.conditionVirus")), alpha = 0.05)
testgroup <- "kids_with_virus_"
my_go_script(resG, testgroup, 0.05)
mykegg(resG, testgroup)

####################################### ~ Volcano Plots with enhancedVolcano####
#This is for the With treatment what is the difference between the kid and the adult. ~ resCon
#resLFC <- lfcShrink(ddsKEEP, coef = 2 , type = 'apeglm')
#This is for the second co-efficient of the design which is "age_Kid_vs_Adult" without the effect of treatment ~ resCon
source('myvocanoplot.R')
ddsCon <- ddsG
resultsNames(ddsCon)
resAge <- results(ddsCon, tidy = F, name = "age_Kid_vs_Adult", alpha = 0.05)
myvolcano(resAge, '_age_kid_vs_adult_without_rx_')
resAge <- results(ddsCon, tidy = F, name = "condition_Virus_vs_Mock", alpha = 0.05)
myvolcano(resAge, '_adult_without_rx_')
resAge <- results(ddsCon, tidy = F, list(c("age_Kid_vs_Adult", "ageKid.conditionVirus")), alpha = 0.05)
myvolcano(resAge, '_adult_vs_kid_rx_')

####################################### ~ Gene Set Enrichment Analysis ~ #############################
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                design = ~ age + condition + age:condition)
dds
baselevel
dds$age <- relevel(dds$age, ref = baselevel)
dds$age
baselevel2 <- "Mock"
dds$condition <- relevel(dds$condition, ref = baselevel2)
dds$condition
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]
ddsKEEP <- DESeq(ddsKEEP)
###Difference between kids and adults without treatment effect
source('myfgseaplot.R')
resAge <- results(ddsKEEP, tidy = T, name = "age_Kid_vs_Adult", alpha = 0.05)
myfgsea(resAge, 'age_kid_vs_adult_without_rx_')
###With treatment what is the different between Kids and adults
resAge <- results(ddsKEEP, tidy = T, list(c("age_Kid_vs_Adult", "ageKid.conditionVirus")), alpha = 0.05)
myfgsea(resAge, '_age_kid_vs_adult_with_rx_')

####################################### ~ Whole transcriptome PCA analysis ~ ####################################
dir.create(paste0(resultsdir))
#Adult_vs_Kids_Mock_samples
MyreadCountMatrix <- read.table(paste0(resultsdir, '/', "analysis1_hsapiens_normalized_counts.txt"), 
                                sep = '\t', header = T)

head(MyreadCountMatrix)
MyreadCountMatrix

df <- (dplyr::select(MyreadCountMatrix, -X))

head(df)
dim(df)
df1 <- df[,1:20]
dim(df1)
dfmock1 <- df1[,!c(TRUE, FALSE)]
dim(dfmock1)
df2 <- df[,21:25]
df2

dfmock <- cbind(dfmock1, df2)
dfmock
metadata <- data.frame(row.names = colnames(df))
metadata$Group <- rep(NA, ncol(df))
metadata
metadata$Group[seq(1,20,2)] <- 'Virus'
metadata$Group[seq(2,20,2)] <- 'Mock'
metadata[c(21:25),] <- 'Mock'
metadata
metadata$CRP <- sample.int(100, size=ncol(df), replace=TRUE)
metadata$ESR <- sample.int(100, size=ncol(df), replace=TRUE)

metadata
metamock <- metadata[metadata$Group == "Mock",]

metadata <- metamock
metadata

metadata$age <- "Adult" 
metadata[6:10,]$age <- "Kid"
metadata[11:14,]$age <- 'Kid'
metadata[15,]$age <- 'Adult'
metadata
pca(dfmock, metadata = metadata)
p <- pca(dfmock, metadata = metadata)

p$rotated
df <- p$rotated
dim(df)
dim(metadata)
metadata$age
df$Age <- metadata$age
df$Age
df
df$Age
resultsdir
pca1_df <- df
write.table(pca1_df, file = paste0(resultsdir, '/pca1_data_frame.txt'), sep = '\t')
library(ggalt)

pca <- ggplot(df, aes(x = PC1, y = PC2, fill = Age))+
  geom_point(aes(colour = Age)) + geom_text(aes(label=rownames(df), hjust=0.5, vjust=0.5)) +
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title = element_text(colour = "black", size = 12, face = "bold"))+
  geom_encircle(alpha = 0.2, show.legend = FALSE)
resultsdir
print(pca)
pdf(file = paste0(resultsdir, "/pca_kids_and_adults_mock.pdf"), width = 6, height = 5)
print(pca)
dev.off()
df <- df %>%
  mutate(color = ifelse(Age == "Kid", "orangered2",
                        ifelse(Age == "Adult", "lightslateblue", "none")))

df <- df[order(df$Age, decreasing = T), ]
df$color <- as.factor(df$color)
col <- as.character(df$color)
pca1 <- ggplot(df, aes(x = PC1, y = PC4, color = color))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = c(0.9, 0.8),
        text = element_text(size = 12),      
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),#element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.background = element_blank(),#element_rect(fill = "grey", size = 0.5),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks = NULL) +  geom_text(aes(label=rownames(df))) +
  geom_encircle(alpha = 0.2, show.legend = FALSE) 
  
pca1 <- pca1 + geom_point(size = 4) + scale_color_manual(values = c("lightslateblue", "orangered2"), 
                                         labels = c("Adult", "Paediatric")) 
print(pca1)
#Adult_vs_Kid_virus_samples
df <- (dplyr::select(MyreadCountMatrix, -X))
df1 <- df[,1:20]
dfmock1 <- df1[,!c(FALSE, TRUE)]
head(dfmock1)
metadata

metadata <- data.frame(row.names = colnames(df1))
metadata
metadata$Group <- rep(NA, ncol(df1))
metadata
metadata$Group[seq(1,20,2)] <- 'Virus'
metadata$Group[seq(2,20,2)] <- 'Mock'
metadata$CRP <- sample.int(100, size=ncol(df1), replace=TRUE)
metadata$ESR <- sample.int(100, size=ncol(df1), replace=TRUE)

metadata
metavirus <- metadata[metadata$Group == "Virus",]

metadata <- metavirus
metadata

metadata$age <- "Adult" 
metadata[6:10,]$age <- "Kid"
metadata

pca(dfmock1, metadata = metadata)
p <- pca(dfmock1, metadata = metadata)

p$rotated
df <- p$rotated
df
df$Age <- rep("Adult", 10)
df[6:10,]$Age <- "Paediatric" 
df$Age
pca2df <- df
write.table(pca2df, file = paste0(resultsdir, '/pca2_data_frame.txt'), sep = '\t')
pca <- ggplot(df, aes(x = PC1, y = PC4, fill = Age)) + geom_text(aes(label=rownames(df))) +
  geom_point(aes(colour = Age)) +
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title = element_text(colour = "black", size = 12, face = "bold"))+
  geom_encircle(alpha = 0.2, show.legend = FALSE)


print(pca)
pdf(file = paste0(resultsdir, "/pca_kids_and_adults_virus.pdf"), width = 6, height = 5)
print(pca)
dev.off()

df <- df %>%
  mutate(color = ifelse(Age == "Paediatric", "orangered2",
                        ifelse(Age == "Adult", "lightslateblue", "none")))

df <- df[order(df$Age, decreasing = T), ]
df$color <- as.factor(df$color)
col <- as.character(df$color)
pca2 <- ggplot(df, aes(x = PC1, y = PC4, color = color))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = c(0.9, 0.8),
        text = element_text(size = 12),      
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),#element_text(color = "black", size = 12, face = "bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.background = element_blank(),#element_rect(fill = "grey", size = 0.5),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks = NULL) +
  geom_encircle(alpha = 0.2, show.legend = FALSE) 

pca2 <- pca2 + geom_point(size = 4) + scale_color_manual(values = c("lightslateblue", "orangered2"), 
                                                 labels = c("Adult", "Paediatric")) 
print(pca2)

####################################### ~ Combining plots ~ ########################################
#source("combine_lots_analysis1.R")
#source('GO_enrichment_plots.R')

####################################### ~ MY Interaction plots ##################################
source("myineract_plot.R") 
source("pilot_config.R") #Interaction analysis 
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                design = ~ age + condition + age:condition)


design(dds) #Verify desing 
ddsmult <- estimateSizeFactors(dds)
ddsres <- DESeq(ddsmult)

results(ddsres, tidy=TRUE) %>%
  arrange(padj) %>%
  head(20) ->
  tophits

tophits

topgenes <- tophits$row[c(1,2,20)]
topgenes

#Visualize interaction
vstexp <- vst(ddsmult, blind=TRUE) #Normalize counts with vst method
p1 <- myinteractplot(vstexp, topgenes[1], "age")
p2 <- myinteractplot(vstexp, topgenes[2], "age")
p3 <- myinteractplot(vstexp, topgenes[3], "age")
grid.arrange(p1, p2, p3, ncol = 3)

#Estimate treatment effect size within each age group
dds <- ddsmult
dds$group <- as.factor(paste0(dds$age,dds$condition))
design(dds) <- ~ group
ddsgrp <- DESeq(dds)

colData(ddsgrp)

ddsmult@design
dds@design

### Results for interaction analysis
r0<-results(ddsres, tidy=TRUE)
### Contrast treatments within "Adult"
r1<-results(ddsgrp, contrast = c("group", "AdultVirus", "AdultMock"), tidy=TRUE)
### Contrast treatments within "Kid"
r2<-results(ddsgrp, contrast = c("group", "KidVirus", "KidMock"), tidy=TRUE)
colData(ddsgrp)
r3 <- results(ddsgrp, contrast = c("group", "KidMock", "AdultMock"), tidy = T)
r3 %>%
  arrange(padj) %>%
  head(20)

r4 <- results(ddsgrp, contrast = c("group", "KidVirus", "AdultVirus"), tidy = T)
r4 %>%
  arrange(padj) %>%
  head(20)

head(r1)
### Look at the results for top two genes
r0 %>% filter(row %in% topgenes)
r1 %>% filter(row %in% topgenes)
r2 %>% filter(row %in% topgenes)
####################################### ~ end of analysis ####
graphics.off()
rm(list = ls())

