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
####################################### ~ DECLARE YOUR VARIABLES HERE ~ #####
myspecies <- "H.sapiens"
baselevel <-  "Mock" 
baselevel2 <- "Virus"
ageclass <- "analysis2"
testgroup <- "/salmon_res/"

mygenes <- c("MYD88", "TRAIP", "LRRFIP1", "PIK3IP1", "TRAM1", "DDX58", "TLR3") #, "TLR7") 
#"TICAM1", "IRF3", "IRAK4", "IRAK2", "TRAF6", "IFIH1", "DHX58",
#"MAVS")
directory <- "D:/17.SARSC-V2_experiment/"
setwd("D:/17.SARSC-V2_experiment/R-project/SCVO2_KVSA/")
resultsdir <- paste0(directory, "results")
dir.create(paste0(resultsdir))
dir.create(paste0( resultsdir, testgroup))
resultsdir <- paste0(resultsdir,testgroup)
dir.create(paste0(resultsdir, ageclass))
resultsdir <- paste0(resultsdir, ageclass)
resultsdir

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
head(sample_table)
sampleCondition<-pull(sample_table, fac)
sampleCondition
sampleCondition1 <- pull(sample_table,fac1)
sampleCondition1
sampleCondition2 <- pull(sample_table, fac2)
sampleCondition2
length(sampleFiles)

sampleTable <- data.frame(sampleName = names(sampleFiles), fileName = sampleFiles, 
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

#import data into DESeq2
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                design = ~ condition)

dds

baselevel <- "Adult"
dds$age <- relevel(dds$age, ref = baselevel)
dds$age
dds$condition <- relevel(dds$condition, ref = "Mock")
dds$condition
keep <- rowSums(counts(dds)) > 10

ddsKEEP <- dds[keep,]

ddsKEEP

ddsKEEP <- estimateSizeFactors(ddsKEEP)
ddsKEEP <- estimateDispersions(ddsKEEP)
ddsKEEP <- nbinomWaldTest(ddsKEEP)

colData(ddsKEEP)

resultsNames(ddsKEEP)
ageclass <- 'analysis2'
nonnormalized_counts <- counts(dds)
write.table(nonnormalized_counts, file=paste(resultsdir, "/", ageclass, "sarscov2_non_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

normalized_counts <- counts(ddsKEEP, normalized=TRUE)
head(normalized_counts)
t <- as.data.frame(normalized_counts)
head(t)
colnames(t)
ageclass
resultsdir
write.table(t, file=paste(resultsdir, "/", ageclass, "_hsapiens_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

head(results(ddsKEEP, tidy = T, name = "condition_Virus_vs_Mock"))


#What is the difference between Kid and Age without treatment?
resCon <- results(ddsKEEP, tidy = F, name = "condition_Virus_vs_Mock")
resCon$log2FoldChange
(testgroup <- ageclass)

resultsNames(ddsKEEP)
resSig= resCon[which(resCon$padj<0.1),]
resSig

(resSig[order(resSig$log2FoldChange),])
write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", testgroup, "_DEGS_hsapiens.csv", sep = "") )

(resSig_up= resSig[which(resSig$log2FoldChange > 2),])
resSig_up=resSig_up[order(resSig_up$log2FoldChange),]
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", testgroup,"_UP_DEGS_hsapiens.csv", sep = "" ))


resSig_down= resSig[which(resSig$log2FoldChange < -2),]
resSig_down=resSig_down[order(resSig_down$log2FoldChange),]
head(resSig_down)
write.csv( as.data.frame(resSig_down), file=paste(resultsdir, "/", testgroup, "_DOWN_DEGS_hsapiens.csv", sep = ""))

####################################### ~ MA plots ~ ####
plotMA(resCon,ylim=c(-2,2))
plotMA(resCon,ylim=c(-5,5))
#plotMA(resCon,ylim=c(-2,2))

plotDispEsts( ddsKEEP, ylim = c(1e-6, 1e1) )
hist( resCon$pvalue, breaks=20, col="green" )

resLFC <- lfcShrink(ddsKEEP, coef = 2, type="apeglm")
resLFC
plotMA(resLFC, ylim = c(-1,1))
abline(h=c(-0.5,0.5), col="red")

# create bins using the quantile function
qs <- c( 0, quantile( resCon$baseMean[resCon$baseMean > 0], 0:7/7 ) )

# "cut" the genes into the bins
bins <- cut( resCon$baseMean, qs )

# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))

# calculate the ratio of £p£ values less than .01 for each bin
ratios <- tapply( resCon$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )

# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

####################################### ~ PCA plots ~ ####
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
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
(topGene <- rownames(resCon)[which.min(resCon$padj)])

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
  #scale_fill_brewer(palette="BuPu") +
  ggtitle(paste(myspecies, testgroup, ": gene -", mygene)) + xlab(testgroup) + ylab("Noramlized gene count") +labs(fill = sampleCondition) +
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
resD <- results(ddsKEEP, name = "condition_Virus_vs_Mock", alpha = 0.1)
resD
summary(resD)

resDSort <- resD[order(resD$padj),]
resDSort

topDESeq2 <- resDSort[1:395,]
topDESeq2
write.csv(topDESeq2, file=paste(resultsdir, "/", testgroup, "_hsapiens_topDESeq2.csv", sep = ""))

(topgenes <- head(rownames(resDSort),200))
mat <- assay(rld)[topgenes,]
(mat <- mat -rowMeans(mat))

col.pan <- colorpanel(100, "blue","white","red")

#Non scaled heatmap for topgenes
heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none", trace="none", labRow= "", labCol = sampleCondition)

#Scaled heatmap for topgenes
scaled.mat<-t(scale(t(mat)))
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

# the call to DESeqTransform() is needed to trigger our plotPCA method.
plotPCA( DESeqTransform( se ) )

#Custom transformation with summerized experiment
pdf(paste0(resultsdir, "/", testgroup,"_SE_PCA_plot.pdf"))
plotPCA( DESeqTransform( se ) )
dev.off()

head( order( rowVars( assay(rld) ), decreasing=F ), 25)
topVarGenes <- head(order( rowVars( assay(rld) ), decreasing = T ), 25)

#heatmap.2( assay(rld)[ topVarGenes, ],scale="row",trace="none", dendrogram="column",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(256))

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

##################################### ~ GO Enrichment analysis for NCBI based annotations ~ #################################
head(tx2gene)

count_data <- tximport(files = sampleFiles,
                       type = "salmon",
                       tx2gene = tx2gene,
                       ignoreTxVersion = FALSE)

head(count_data)

#import data into DESeq2

dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                design = ~ condition)
dds
baselevel
dds$age <- relevel(dds$age, ref = baselevel)
baselevel2 <- "Mock"
dds$condition <- relevel(dds$condition, ref = baselevel2)
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

#######***************************Biomart GO::terms matching*******************************########################
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

library(GenomicFeatures)
library(goseq)

txdb = makeTxDbFromGFF(GTF, format = "gtf", )

(txdb.allcolumns <- transcripts(txdb))

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

library(GO.db)
capture.output(for(go in enriched.GO[1:88]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file=paste0(resultsdir, "/", testgroup,"_SigGo.txt"))

goResults[goResults$ontology == 'BP',]
goBPRes <- goResults[goResults$ontology == 'BP',] 


##########*************************visualize the top 30 hits*****************************##########################
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

####################################Gene Set Enrichment Analysis#############################
dds <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sampleTable,
                                design = ~ condition)
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
resCon <- results(ddsKEEP, tidy = T, name = "condition_Virus_vs_Mock")

res2 <- resCon %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat))

res2
ranks <- tibble::deframe(res2)
head(ranks, 20)

library(fgsea)
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

p <-ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

print(p)
ggsave(p, file=paste0(resultsdir, "/", testgroup, "hallmark_fgsea_effect_of_treatment_irrespective_of_age",".png", sep = ""), width = 20, height = 28, units = "cm")
pdf(paste0(resultsdir, "/", testgroup,"hallmark_fgsea_without_treatment.pdf"))
print(p)
dev.off()

#### ~ end of analysis ####
graphics.off()
rm(list = ls())
