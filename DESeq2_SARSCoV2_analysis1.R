####Difference between age groups with or without treatment effect####
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

sampleFiles <- paste0("../../SCOV2_Study_19122020/salmon_res/",pull(sample_table, `ID`),"/", "quant.sf")
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

sampleTable <- data.frame(sampleName = names(sampleFiles), fileName = sampleFiles, 
                          age = as.character(sampleCondition), condition = as.character(sampleCondition1),
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
                                design = ~ age + condition + age:condition)
                                
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
resAge <- results(ddsKEEP, tidy = F, name = "age_Kid_vs_Adult", alpha = 0.05)
res <- resAge
DGE.df <- as.data.frame(res)
DGE.df
getwd()
write.csv(DGE.df, file = "differential_kid_vs_adult_without_treatment.csv")
res$log2FoldChange
sum(res$padj<0.05, na.rm = T)
dge res[which(res$padj < 0.05),]
table(is.na(res$padj))
resultsNames(ddsKEEP)

resSig= res[which(res$padj<0.05),] #Significant genes
testgroup <- ageclass
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_kid_vs_adult_without_virus_DEGS_hsapiens.csv", sep = "") )
resultsdir
(resSig_up= resSig[which(resSig$log2FoldChange > 2),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_kid_vs_adult_without_virus_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -2),]
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
(resSig_up= resSig[which(resSig$log2FoldChange > 2),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_kid_vs_adult_with_viruseffect_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -2),]
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
(resSig_up= resSig[which(resSig$log2FoldChange > 2),])
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_kid_effect_of_virus_UP_DEGS_hsapiens.csv", sep = "" ))
resSig_down= resSig[which(resSig$log2FoldChange < -2),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_kid_effect_of_virus_DOWN_DEGS_hsapiens.csv", sep = ""))
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
rld <- rlog(ddsKEEP)
class(rld)

pdf(paste0(resultsdir, "/", testgroup, "_rlog_PCA_plot.pdf"))
plotPCA(rld)
pcaplot(rld, pcX = 2, pcY = 1, ntop = 1000, intgroup = c("age", "condition", "pairID"))
pcaplot3d(rld)

pcaplot(rld, pcX = 4, pcY = 3, ntop = 5)
pcaplot(rld, pcX = 4, pcY = 3, ntop = 500)
pcaplot(rld, pcX = 3, pcY = 1)
pcaplot(rld, pcX = 3, pcY = 2)

dev.off()

rld
graphics.off()
plotPCA(rld, returnData = F)

pcaplot(rld, pcX = 3, pcY = 4)
pcaplot(rld, pcX = 4, pcY = 3)
pcaplot(rld, pcX = 1, pcY = 3)
pcaplot(rld, pcX = 3, pcY = 1)
pcaplot(rld, pcX = 2, pcY = 3)
pcaplot(rld, pcX = 3, pcY = 2)
pcaplot(rld, pcX = 1, pcY = 2)

####################################### ~ Differential gene expression of given genes ~ ######################################
colData(ddsKEEP)
ddsKEEP$condition
(topGene <- rownames(resAge)[which.min(resAge$padj)])

plotCounts(ddsKEEP, gene=topGene, intgroup=c("condition", "age"), returnData = T)

mygene <- topGene
d <- plotCounts(ddsKEEP,gene = mygene, intgroup = "age", main = paste( myspecies, "gene -", mygene), returnData=TRUE)
p <- ggplot(d, aes(x=age, y=count)) + 
  geom_boxplot(colour = "red", fill = "orange", alpha = 0.2, 
               outlier.colour="red", outlier.shape=8, outlier.size=2, notch=F, notchwidth = 0.5) + 
  geom_point(position=position_jitter(w=0.1,h=0), colour = 'purple', size = 1) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme(
    #panel background elements
    panel.background = element_rect(
      fill = "grey90",
      colour = "black",
      size = 1,
    ),
    legend.position= "bottom",
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=12, face="bold"),
    axis.title.y = element_text(color="#993333", size=12, face="bold")
  ) + 
  ggtitle(paste(myspecies, testgroup, ":condition gene -", mygene)) + xlab(testgroup) + ylab("Noramlized gene count") +labs(fill = sampleCondition) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) + scale_color_grey() +
  scale_fill_manual(values=c("#999999", "#E69F00")) 
print(p)
ggsave(p, file=paste0(resultsdir, "/", testgroup,"_", mygene,".png", sep = ""), width = 14, height = 10, units = "cm")

testgroup
sampleCondition

for (mygene in mygenes) {
  print(mygene)
  d <- plotCounts(ddsKEEP,gene = mygene, intgroup = "age", main = paste( myspecies, "gene -", mygene), returnData=TRUE)
  d
  p <- ggplot(d, aes(x=age, y=count)) + 
    #geom_violin()  +
    geom_boxplot(colour = "red", fill = "orange", alpha = 0.2, 
                 outlier.colour="black", outlier.shape=8, outlier.size=2, notch=F, notchwidth = 0.5) + 
    #geom_dotplot(binwidth = 50, binaxis='y', stackdir='center', dotsize=1)
    geom_point(position=position_jitter(w=0.1,h=0), colour = 'purple', size = 2) + 
    scale_y_log10(breaks=c(25,100,400)) + 
    theme(
      #panel background elements
      panel.background = element_rect(
        fill = "grey90",
        colour = "black",
        size = 1,
      ),
      legend.position= "bottom",
      plot.title = element_text(color="black", size=14, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=12, face="bold"),
      axis.title.y = element_text(color="#993333", size=12, face="bold")
    ) + 
    ggtitle(paste(myspecies, "gene -", mygene)) + xlab(testgroup) + ylab("Noramlized gene count") +labs(fill = sampleCondition) +
    stat_summary(fun=mean, geom="point", shape=23, size=4) + scale_color_grey() +
    scale_fill_manual(values=c("#999999", "#E69F00"))
  ggsave(p, file=paste0(resultsdir, "/", testgroup, mygene,".png", sep = ""), width = 14, height = 10, units = "cm")
  print(p)
}

####################################### ~ Getting the top genes ~ ####
resultsNames(ddsKEEP)
resD <- results(ddsKEEP, name = "age_Kid_vs_Adult", alpha = 0.1)
resD
summary(resD)

resDSort <- resD[order(resD$padj),]
resDSort
topDESeq2 <- resDSort[1:395,]
topDESeq2
write.csv(topDESeq2, file=paste(resultsdir, "/", testgroup, "_hsapiens_topDESeq2.csv", sep = ""))

(topgenes <- head(rownames(resDSort),40))
mat <- assay(rld)[topgenes,]
(mat <- mat -rowMeans(mat))

col.pan <- colorpanel(100, "blue","white","red")
#Non scaled heatmap for topgenes
heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none", trace="none", labRow= "", labCol = sampleCondition)

#Scaled heatmap for topgenes
scaled.mat<-t(scale(t(mat)))
t <- scaled.mat[,c(TRUE, FALSE)]
head(scaled.mat)
sampleCondition
dev.off()
heatmap.2(scaled.mat, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", labRow= "",margins = c(10,8),cexRow=0.5, cexCol=1, keysize=1,labCol = sampleCondition)

pdf(paste0(resultsdir, "/", testgroup,"_VSD_scaled_topgenes_heatmap.pdf"))
heatmap.2(scaled.mat, col=col.pan, Rowv=TRUE, scale="none",
          trace="none",margins = c(10,8),cexRow=0.5, cexCol=1, keysize=1,labCol = sampleCondition)
dev.off()
####################################### ~ Another method of PCA plot ~ ####
pdf(paste0(resultsdir, "/", testgroup,"_VSD_PCA_plot.pdf"))
vsdata <- vst(ddsKEEP, blind=FALSE)
plotPCA(vsdata, intgroup="age")
pcaplot3d(vsdata)
dev.off()

plotPCA(vsdata, intgroup="condition")
# also possible to perform custom transformation:
ddsEST <- estimateSizeFactors(ddsKEEP)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(ddsEST, normalized=TRUE) + 1),
                           colData=colData(ddsKEEP))

# trigger our plotPCA method.
plotPCA( DESeqTransform( se ) )
#Custom transformation with summerized experiment
pdf(paste0(resultsdir, "/", testgroup,"_SE_PCA_plot.pdf"))
plotPCA( DESeqTransform( se ) )
dev.off()

head( order( rowVars( assay(rld) ), decreasing=F ), 25)
topVarGenes <- head(order( rowVars( assay(rld) ), decreasing = T ), 25)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(100),margin=c(10,8), cexRow=0.5, cexCol=1, keysize=1.5, labCol = sampleCondition)

pdf(paste0(resultsdir, "/", testgroup,"_topVargenes_heatmap.pdf"))
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(100),margin=c(10,8), cexRow=0.5, cexCol=1, keysize=1.5, labCol = sampleCondition)
dev.off()

#Heatmap of the count matrix

select <- order(rowMeans(counts(ddsKEEP,normalized=TRUE)),
                decreasing=TRUE)[1:20]

select

#Heatmap of the sample-to-sample distances
sampleDists <- dist( t( assay(vsdata) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
sampleDistMatrix
rownames(sampleDistMatrix) <- paste( vsdata$condition,
                                     vsdata$sizeFactor, sep="-" )
colnames(sampleDistMatrix) <- NULL

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( 
  sampleDistMatrix, trace="none", col=colours,margin=c(10,8), cexRow=0.5, cexCol=1, keysize=1.5, labRow = sampleCondition)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colours)

pdf(paste0(resultsdir, "/", testgroup,"_sample_to_sample_distance.pdf"))

heatmap.2( 
  sampleDistMatrix, trace="none", col=colours,margin=c(10,8), cexRow=0.5, cexCol=1, keysize=1.5, labRow = sampleCondition)
dev.off()

ramp <- 1:3/3

cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print( plotPCA( rld ))

pdf(paste0(resultsdir, "/", testgroup,"_sizefac_condition_PCA_plot.pdf"))
plotPCA( rld, intgroup = "condition" )
dev.off()

colData(ddsKEEP)


(topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 20 ))
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margins = c(10,10),labCol = sampleCondition)

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
keep <- rowSums(counts(dds)) >= 10
ddsKEEP <- dds[keep,]

ddsKEEP$sampleName

rowsum.threshold <- 10
fdr.threshold <- 0.1
rs <- rowSums(counts(ddsKEEP))
rs
ddsG <- ddsKEEP[rs > rowsum.threshold,]
ddsG

ddsG<-DESeq(ddsG)
resultsNames(ddsG)

#What is the difference between Kid and Age without treatment?
resG <- results(ddsG, name = "age_Kid_vs_Adult", independentFiltering=FALSE) # use count threshold instead of IF
assayed.genes <- rownames(resG)
assayed.genes

de.genes <- rownames(resG)[ which(resG$padj < fdr.threshold) ]

(de.genes)
options(width = 84)

class(assayed.genes)
gene.vector=as.integer(assayed.genes%in%de.genes)

gene.vector

(names(gene.vector)=assayed.genes)
(gene.vector)
names(gene.vector)

####################################### ~ Biomart GO::terms matching*******************************########################
# define biomart object
ensembl<- useMart("ensembl")

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

sp <- useDataset("hsapiens_gene_ensembl",mart = ensembl)
sp

EG2KEGG<- getBM(mart = sp, values = gene.vector, attributes = c("hgnc_symbol", "go_id"))
head(EG2KEGG)

geneID2GO <- by(EG2KEGG$hgnc_symbol,
                EG2KEGG$go_id,
                function(x) as.data.frame(x))
head(geneID2GO,15)
class(geneID2GO)
write.csv(as.data.frame(EG2KEGG), file = paste0(resultsdir, "/", testgroup, "_geneID.csv"))

GTF <- "../../gencode.v36.annotation.gtf"

txdb = makeTxDbFromGFF(GTF, format = "gtf", )

(txdb.allcolumns <- transcripts(txdb))

gene.vector

lengthdata <- getlength(genes = names(gene.vector), "hg19",  "geneSymbol")
(names(lengthdata) <- names(gene.vector))

pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
plotPWF(pwf,binsize = 1000)

goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)

class(goResults)
head(goResults)
nrow(goResults)

write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", testgroup, "_gores.csv" ))
enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", testgroup, "enriched_go_", myspecies, ".csv"))

capture.output(for(go in enriched.GO[1:88]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file=paste0(resultsdir, "/", testgroup,"_SigGo.txt"))

goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 

####################################### ~ visualize the top 30 hits*****************************##########################
pdf(paste0(resultsdir, "/", testgroup,"_goseq_enrichment.pdf"))
goResults %>% 
  top_n(40, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

goBPRes %>% 
  top_n(40, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term: BP", colour="p value", size="Count")
dev.off()

goResults %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

goBPRes %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term:BP", colour="p value", size="Count")


write.csv(as.data.frame(goBPRes), file = paste0(resultsdir,"/", testgroup, myspecies, "_go_bp.csv"))

####################################### ~ Volcano Plots with enhancedVolcano####
ddsKEEP
resultsNames(ddsKEEP)
#This is for the With treatment what is the difference between the kid and the adult. ~ resCon
#resLFC <- lfcShrink(ddsKEEP, coef = 2 , type = 'apeglm')
#This is for the second co-efficient of the design which is "age_Kid_vs_Adult" without the effect of treatment ~ resCon
resultsNames(ddsCon)
resAge <- results(ddsCon, tidy = F, name = "age_Kid_vs_Adult", alpha = 0.05)
res <- resAge
res$log2FoldChange
keyvals <- ifelse(
  res$log2FoldChange < -1.5, 'royalblue',
  ifelse(res$log2FoldChange > 1.5, 'gold',
         'black'))

keyvals
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

res$padj
vol1 <- EnhancedVolcano(res,
                        lab = NA, #rownames(res),
                        x = 'log2FoldChange',
                        y = 'padj', 
                        selectLab = NULL,#rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                        #selectLab = c('OASL','CXCL10','ISG15', 'USP41','IFITM1','MX1'),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        title = NULL, #'Volcano plot - Difference in Kids compared to Adults ',
                        subtitle = NULL, #'~ without treatment effect(i.e. without Virus infection)',
                        pCutoff = 0.05,
                        FCcutoff = 1.5,
                        pointSize = 4.0,
                        labSize = 3.0,
                        labCol = 'black',
                        labFace = 'bold',
                        boxedLabels = F,
                        ylim = c(0,2.5),
                        xlim = c(-8,8),
                        colAlpha = 3/4,
                        gridlines.major = F,
                        gridlines.minor = F,
                        col = c("grey30", "palegreen3", "lightslateblue", "orangered2"),
                        legendPosition = c(0.9, 0.9),
                        legendLabels = c("NS", bquote(~Log[2]~'fold change (±1.5)'), 'p-adj < 0.05', bquote(~Log[2]~ 'fold change & p-adj')),
                        legendLabSize = 12,
                        legendIconSize = 4.0,axisLabSize = 12,
                        drawConnectors = F,
                        widthConnectors = 1.0, caption = NULL,
                        colConnectors = 'black') 

vol1 = vol1 + ggplot2::theme(axis.text.y = element_blank(),text = element_text(family = 'serif')) + scale_y_continuous(breaks=NULL) 
print(vol1)

dds
keep <- rowSums(counts(dds)) > 10
ddsCon <- dds[keep,]
ddsCon <- DESeq(ddsCon)
####################################### ~ EnhancedVolcano: Effect of virus in between kids and adults #####
resCon <- results(ddsCon, tidy = F, list(c("age_Kid_vs_Adult", "ageKid.conditionVirus")), alpha = 0.05)
res <- resCon
res$log2FoldChange
keyvals <- ifelse(
  res$log2FoldChange < -1.5, 'royalblue',
  ifelse(res$log2FoldChange > 1.5, 'gold',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

vol2 <- EnhancedVolcano(res,
                        lab = NA, #rownames(res),
                        x = 'log2FoldChange',
                        y = 'padj', 
                        selectLab = NULL, #rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                        #selectLab = c('OASL','CXCL10','ISG15', 'USP41','IFITM1','MX1'),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        title = NULL, #'Volcano plot - Difference in Kids compared to Adults ',
                        subtitle = NULL, #'~ without treatment effect(i.e. without Virus infection)',
                        pCutoff = 0.05,
                        FCcutoff = 1.5,
                        pointSize = 3,
                        labSize = 3.0,
                        labCol = 'black',
                        labFace = 'bold',
                        boxedLabels = F,
                        ylim = c(0,2.5),
                        xlim = c(-10,10),
                        colAlpha = 3/4,
                        gridlines.major = F,
                        gridlines.minor = F,
                        col = c("grey30", "palegreen3", "lightslateblue", "orangered2"),
                        legendPosition = c(0.85, 0.7),
                        legendLabels = c("NS", bquote(~Log[2]~'fold change (±1.5)'), 'p-adj < 0.05', bquote(~Log[2]~ 'fold change & p-adj')),
                        legendLabSize = 8,
                        legendIconSize = 2.5,
                        drawConnectors = F,
                        widthConnectors = 1.0, caption = NULL,
                        colConnectors = 'black') 

vol2
vol2 = vol2 + ggplot2::theme(axis.text.y = element_blank(), axis.title = element_text(size = 12),
                             text = element_text(family = 'serif', face = 'bold')) + 
  scale_y_continuous(breaks=NULL) 
print(vol2)

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

ddsKEEP$sampleName

ddsKEEP <- DESeq(ddsKEEP)
###Difference between kids and adults without treatment effect
resAge <- results(ddsKEEP, tidy = T, name = "age_Kid_vs_Adult", alpha = 0.05)

res2 <- resAge %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat))

res2
ranks <- tibble::deframe(res2)
head(ranks, 20)

pathways.hallmark <- gmtPathways("../../h.all.v7.2.symbols.gmt")

pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks) #nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

t <- fgseaResTidy[fgseaResTidy$padj <= 0.05,]
t$pathway
t$leadingEdge[[17]]
t$leadingEdge[[18]]
t$leadingEdge[[23]]

p1 <-ggplot(t, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme(
    legend.box      = "horizontal",
    legend.key      = element_blank(),
    legend.title    = element_text(family = 'serif'),
    axis.title.y = element_text(family = 'serif'),
    title = element_text(family = 'serif', size = 12, face = 'bold'),
    axis.text.y = element_text(family = 'serif', size = 12),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank()
  ) +
  theme_bw() +
  scale_fill_gradient(low="royalblue", high="royalblue4") +
  theme(legend.position = c(0.8, 0.2))

print(p1)

t <- head(t, 20)
p1 <-ggplot(t, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme(
    legend.box      = "horizontal",
    legend.key      = element_blank(),
    legend.title    = element_text(family = 'serif'),
    axis.title.y = element_text(family = 'serif'),
    title = element_text(family = 'serif', size = 12, face = 'bold'),
    axis.text.y = element_text(family = 'serif', size = 12),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank()
  ) +
  theme_bw() +
  scale_fill_gradient(low="royalblue", high="royalblue4") +
  theme(legend.position = c(0.8, 0.2))

ggsave(p1, file=paste0(resultsdir, "/", testgroup, "hallmark_fgsea_difference_between_age_without_treatment_effect",".png", sep = ""), width = 20, height = 28, units = "cm")
pdf(paste0(resultsdir, "/", testgroup,"hallmark_fgsea_difference_between_age_without_treatment_effect.pdf"))
print(p1)
dev.off()
tiff(file = paste0(resultsdir, "/", testgroup, "hallmark_fgsea_difference_between_age_without_treatment_effect.tiff"), width = 3200, height = 3200, units = "px", res = 800)
print(p1)
dev.off()

p1
###With treatment what is the different between Kids and adults
resCon <- results(ddsKEEP, tidy = T, list(c("age_Kid_vs_Adult", "ageKid.conditionVirus")), alpha = 0.05)

res2 <- resCon %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat))

res2
ranks <- tibble::deframe(res2)
head(ranks, 20)

pathways.hallmark <- gmtPathways("../../h.all.v7.2.symbols.gmt")

pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks) 

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

t <- fgseaResTidy[fgseaResTidy$padj <= 0.05,]

t$pathway
t$leadingEdge[[21]]

p2 <-ggplot(t, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme(
    legend.box      = "horizontal",
    legend.key      = element_blank(),
    legend.title    = element_text(family = 'serif', face = 'bold'),
    axis.title.y = element_text(family = 'serif',hjust = 0.8),
    title = element_text(family = 'serif', size = 10, face = 'bold'),
    axis.text.y = element_text(family = 'serif', size = 8),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank()
  ) +
  theme_bw() +
  scale_fill_gradient(low="royalblue", high="royalblue4") +
  theme(legend.position = c(0.8, 0.2))

print(p2)


t <- head(t, 20)

p2 <-ggplot(t, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme(
    legend.box      = "horizontal",
    legend.key      = element_blank(),
    legend.title    = element_text(family = 'serif', face = 'bold'),
    axis.title.y = element_text(family = 'serif',hjust = 0.8),
    title = element_text(family = 'serif', size = 10, face = 'bold'),
    axis.text.y = element_text(family = 'serif', size = 8),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank()
  ) +
  theme_bw() +
  scale_fill_gradient(low="royalblue", high="royalblue4") +
  theme(legend.position = c(0.8, 0.2))

print(p2)


ggsave(p2, file=paste0(resultsdir, "/", testgroup, "hallmark_fgsea_difference_across_age_group_with_treatment_effect",".png", sep = ""), width = 20, height = 28, units = "cm")
pdf(paste0(resultsdir, "/", testgroup,"hallmark_fgsea_with_treatment.pdf"))
print(p2)
dev.off()

tiff(file = paste0(resultsdir, "/", testgroup, "hallmark_fgsea_difference_across_age_group_with_treatment_effect.tiff"), width = 3200, height = 3200, units = "px", res = 800)
print(p2)
dev.off()

resultsNames(ddsKEEP)

####################################### ~ Whole transcriptome PCA analysis ~ ####################################
library(PCAtools)
library(dplyr)

ageclass <- "mocks"

directory <- "D:/17.SARSC-V2_experiment/"

setwd("D:/17.SARSC-V2_experiment/R-project/SCVO2_KVSA/")
resultsdir <- paste0(directory, "results")
testgroup <- "within_age_PCA_without_the_effect_of_treatment"

resdir <- "/salmon_res/"

dir.create(paste0(resultsdir))
resultsdir
dir.create(paste0( resultsdir, resdir))

resultsdir <- paste0(resultsdir, resdir)
resultsdir
dir.create(paste0(resultsdir, ageclass))
resultsdir <- paste0(resultsdir, ageclass)
resultsdir


#Adult_vs_Kids_Mock_samples
MyreadCountMatrix <- read.table(paste0(resultsdir, '/', ageclass, "_hsapiens_normalized_counts.txt"), 
                                sep = '\t', header = T)

head(MyreadCountMatrix)
MyreadCountMatrix

df <- (dplyr::select(MyreadCountMatrix, -X))

head(df)
df
dfmock <- df[,!c(TRUE, FALSE)]
head(dfmock)

metadata <- data.frame(row.names = colnames(df))
metadata$Group <- rep(NA, ncol(df))
metadata
metadata$Group[seq(1,20,2)] <- 'Virus'
metadata$Group[seq(2,20,2)] <- 'Mock'
metadata$CRP <- sample.int(100, size=ncol(df), replace=TRUE)
metadata$ESR <- sample.int(100, size=ncol(df), replace=TRUE)

metadata
metamock <- metadata[metadata$Group == "Mock",]

metadata <- metamock
metadata

metadata$age <- "Adult" 
metadata[6:10,]$age <- "Kid"
metadata

pca(dfmock, metadata = metadata)
p <- pca(dfmock, metadata = metadata)

p$rotated
df <- p$rotated
df$Age <- rep("Adult", 10)
df[6:10,]$Age <- "Paediatric" 
df$Age
resultsdir
pca1_df <- df
write.table(pca1_df, file = paste0(resultsdir, '/pca1_data_frame.txt'), sep = '\t')
library(ggalt)
pca <- ggplot(df, aes(x = PC1, y = PC4, fill = Age))+
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

df <- df %>%
  mutate(color = ifelse(Age == "Paediatric", "orangered2",
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
  scale_x_continuous(breaks = NULL) +
  geom_encircle(alpha = 0.2, show.legend = FALSE) 
  #stat_ellipse(aes(fill = color), geom = "polygon", level = 0.8, alpha=0.2) +
  
  
pca1 <- pca1 + geom_point(size = 4) + scale_color_manual(values = c("lightslateblue", "orangered2"), 
                                         labels = c("Adult", "Paediatric")) 

print(pca1)

#Adult_vs_Kid_virus_samples
MyreadCountMatrix <- read.table(paste0(resultsdir, '/', ageclass, "_hsapiens_normalized_counts.txt"), 
                                sep = '\t', header = T)

head(MyreadCountMatrix)
MyreadCountMatrix

df <- (dplyr::select(MyreadCountMatrix, -X))

head(df)
df
dfmock <- df[,!c(FALSE, TRUE)]
head(dfmock)

metadata <- data.frame(row.names = colnames(df))
metadata$Group <- rep(NA, ncol(df))
metadata
metadata$Group[seq(1,20,2)] <- 'Virus'
metadata$Group[seq(2,20,2)] <- 'Mock'
metadata$CRP <- sample.int(100, size=ncol(df), replace=TRUE)
metadata$ESR <- sample.int(100, size=ncol(df), replace=TRUE)

metadata
metavirus <- metadata[metadata$Group == "Virus",]

metadata <- metavirus
metadata

metadata$age <- "Adult" 
metadata[6:10,]$age <- "Kid"
metadata

pca(dfmock, metadata = metadata)
p <- pca(dfmock, metadata = metadata)

p$rotated
df <- p$rotated
df
df$Age <- rep("Adult", 10)
df[6:10,]$Age <- "Paediatric" 
df$Age
pca2df <- df
write.table(pca2df, file = paste0(resultsdir, '/pca2_data_frame.txt'), sep = '\t')
pca <- ggplot(df, aes(x = PC1, y = PC4, fill = Age))+
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
#stat_ellipse(aes(fill = color), geom = "polygon", level = 0.8, alpha=0.2) +


pca2 <- pca2 + geom_point(size = 4) + scale_color_manual(values = c("lightslateblue", "orangered2"), 
                                                 labels = c("Adult", "Paediatric")) 
print(pca2)
####################################### ~ Combining plots ~ ########################################
library(ggpubr)
figure1 <-ggarrange(ggarrange(pca1, vol1, nrow = 2), 
          p1, ncol = 2 )


figure1 
#annotate_figure(figure1,
#                top = text_grob("Visualizing Tooth Growth",  size = 14, family = "serif"),
#                bottom = text_grob("ToothGrowth data set", size = 14, family = "Times")
#)

pdf(paste0(resultsdir, "/figure01.pdf"))
print(figure1)
dev.off()

tiff(file = paste0(resultsdir, "/figure1.tiff"), width = 3200, height = 3200, units = "px", res = 800)
print(figure1)
dev.off()
resultsdir
ggsave(figure1, file=paste0(resultsdir, "/figure1",".png", sep = ""), width = 30, height = 34, units = "cm")

graphics.off()
figure2 <-ggarrange(ggarrange(pca2, vol2, nrow = 2), 
                    p2, ncol = 2, widths = c(0.3, 0.7))


#source('GO_enrichment_plots.R')
#figure2 <- ggarrange(figure2, goplots, nrow = 2, heights = c(0.5, 0.5))
figure2

figure2 <-ggarrange(ggarrange(pca2, vol2, nrow = 2), 
                    p2, ncol = 2 )

figure2
#figure2 + theme(text = element_text(size = 8))
pdf(paste0(resultsdir, "/figure02.pdf"))
print(figure2)
dev.off()

tiff(file = paste0(resultsdir, "/figure2.tiff"), width = 3200, height = 3200, units = "px", res = 800)
print(figure2)
dev.off()
resultsdir
ggsave(figure2, file=paste0(resultsdir, "/figure2",".png", sep = ""), width = 30, height = 34, units = "cm")
print(figure2)
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
#graphics.off()
#rm(list = ls())
