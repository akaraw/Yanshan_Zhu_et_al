####################################### ~ Loading libraries ~ #####
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
library(edgeR)
####################################### ~ DECLARE YOUR VARIABLES HERE ~ #####
directory <- "D:/17.SARSC-V2_experiment/"
setwd("D:/17.SARSC-V2_experiment/R-project/SCVO2_KVSA/")
DATADIR <- "D:/17.SARSC-V2_experiment/SCOV2_Study_19122020/salmon_quant/"
METAFILE <- "D:/17.SARSC-V2_experiment/YZ_experiment.txt"
RESDIR <- "D:/17.SARSC-V2_experiment/Davids_folder/limma-voom/"
dir.create(RESDIR, recursive = T)

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

fac1 <- 'Age'
fac2 <- 'Treatment'
fac3 <- 'Gender'

sampleAge<-pull(sample_table, fac1)
sampleCondition <- pull(sample_table, fac2)
sampleGender <- pull(sample_table, fac3)
sampleCondition
sampleAge
sampleGender
length(sampleFiles)

sampleTable <- data.frame(sampleName = names(sampleFiles), fileName = sampleFiles, 
                          condition = as.factor(sampleCondition), age = as.factor(sampleAge),
                          gender = as.factor(sampleGender))
head(sampleTable)
str(sampleTable)
group = as.factor(paste0(sampleTable$age, sampleTable$condition))
group
gender = as.factor(sampleTable$gender)
gender
tx2gene <- read.table("D:/17.SARSC-V2_experiment/tx2genesymbol.txt", header=FALSE, sep="\t")
head(tx2gene)
tx2entrezid <- read.table("D:/17.SARSC-V2_experiment/gencode.v36.metadata.EntrezGene", sep = "\t", header = F)
head(tx2entrezid)
tx2gene <- dplyr::select(tx2gene, -V3)
head(tx2gene)
count_data <- tximport(files = sampleFiles,
                       type = "salmon",
                       tx2gene = tx2gene,
                       ignoreTxVersion = FALSE)

dim(count_data$counts)
cts <- count_data$counts
dim(cts)
normMat <- count_data$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
rownames(cts)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, genes = rownames(cts))
y <- scaleOffset(y, normMat)

#Filtering
keep <- filterByExpr(y, group = group)
y <- y[keep,]

y$samples$group <- group
y$samples$gender <- gender
y$samples$group
y$samples$gender
y$genes$genes

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$genes)
y <- y[!d,]
nrow(y)

y$samples$lib.size <- colSums(y$counts)
y$genes$EntrezGene <- NULL
y$genes$genes
y <- calcNormFactors(y)
y$samples

#The first step of an analysis should be to examine the samples for outliers and for other
#relationships. The function plotMDS produces a plot in which distances between samples
#correspond to leading biological coefficient of variation (BCV) between those samples:
plotMDS(y, top = 10000)

# filtering
#group
y$samples
design <- model.matrix(~0+group+gender, data=y$samples)
colnames(design)
levels(gender)
levels(group)
colnames(design) <- c(levels(group), 'Gen')
design
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
rownames(design) 
y$samples
y$tagwise.dispersion
limma::plotMDS(y)
plotBCV(y)
levels(y$samples$group)

#Classic edgeR method
et <- exactTest(y,pair = c(3,4))
topTags(et)

#counts per millions
cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)
cpms

#Testing for DE genes
group
design <- model.matrix(~0+group, data = y$samples)
design

colnames(design) <- levels(group)
design
#design[,2]
#colnames(design)
#dim(design[2])
rownames(design) <- colnames(y)
design

#Testing for DE genes - testing procedures for determining DE using quasi-likelihood (QL) F-test
fit <- glmQLFit(y, design)
fit

colnames(fit)
AVvsAM <- makeContrasts(AMvsAV = AdultVirus-AdultMock, levels = design)
AVvsAM
qlf <- glmQLFTest(fit, contrast = AMvsAV)
topTags(qlf)

KVvsKM <- makeContrasts(KVvsKM = KidVirus-KidMock, levels = design)
KVvsKM
qlf <- glmQLFTest(fit, contrast = KVvsKM)
topTags(qlf)

#Testing for DE genes - testing procedures for determining DE using likelihood ratio test (LRT)
design <- model.matrix(~group, data = y$samples)
design
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
o <- order(qlf$table$PValue)
cpm(y)[o[1:10],]

summary(decideTests(qlf, lfc = 1.5))
plotMD(fit)
abline(h=c(-1,1), col="blue")

fit <- glmFit(y, design)
colnames(fit)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
#colnames(design)
summary(decideTests(lrt, lfc = 1.5))
plotMD(lrt)
abline(h=c(-1,1), col="blue")

####Starting Limma-Voom####
sample_table
fac <- 'Age'
fac1 <- 'Treatment'
fac2 <- 'Experiment_number'
fac3 <- 'Gender'
sampleCondition<-pull(sample_table, fac)
sampleCondition
sampleCondition1 <- pull(sample_table,fac1)
sampleCondition1
sampleCondition2 <- pull(sample_table, fac2)
sampleCondition2
sampleCondition3 <- pull(sample_table, fac3)
sampleCondition3
length(sampleFiles)

sampleTable <- data.frame(sampleName = names(sampleFiles), fileName = sampleFiles, 
                          gender = as.factor(sampleCondition3), 
                          age = as.character(sampleCondition), condition = as.character(sampleCondition1),
                          pairID = as.factor(sampleCondition2))
group = as.factor(paste(sampleTable$age, sampleTable$condition, sep = '.'))
group
sampleTable$age.cond <- as.factor(paste(sampleTable$age, sampleTable$condition, sep = ".")) 
pairID <- as.factor(sampleTable$pairID)
pairID
head(sampleTable)
str(sampleTable)
cts <- count_data$counts
cts
normMat <- count_data$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
# Combining effective library sizes with the length factors, and calculating offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
rownames(cts)
# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, genes = rownames(cts))
y <- scaleOffset(y, normMat)
#Filtering
keep <- filterByExpr(y, group = group)
y <- y[keep,]
gender
group
y$samples$group <- group
y$samples$gender <- gender
y$samples$pairID <- pairID
y$samples$group
y$samples$gender
y$samples$pairID
y$genes$genes
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$genes)
y <- y[!d,]
nrow(y)
y$samples$lib.size <- colSums(y$counts)
y$genes$EntrezGene <- NULL
y$genes$genes
y <- calcNormFactors(y)
y$samples
group
gender
y$samples$group
design = model.matrix(~0 + group + gender + pairID, data = y$samples)#sampleTable)
par(mfrow=c(1,2))
v <- voom(y,design,plot = TRUE) ##change counts.1y to your data set ###
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
design
colnames(design[,c(1:5)]) <- c(levels(group), 'genderM')
colnames(design)
levels(group)
design
cont.fit = makeContrasts(A=(groupAdult.Virus - groupAdult.Mock),  
                            B=(groupKid.Virus - groupKid.Mock),
                            C=(groupKid.Virus - groupKid.Mock) - (groupAdult.Virus - groupAdult.Mock),
                            levels=design)
cont.matrix
v <- voom(y,design,plot = TRUE)
v
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")

pairID <- y$samples$pairID
pairID
v$design
design
corfit <- duplicateCorrelation(v, design, block=pairID)
#model fit
#Here you are fitting your model and accounting for the within subject correlation#
fit <- lmFit(v,design,correlation=corfit$consensus.correlation, block = pairID)
fit
cfit <- contrasts.fit(fit,contrasts=cont.fit)
cfit <- eBayes(cfit, robust=TRUE)
results <- decideTests(cfit, p.value = 0.05)
results
vennDiagram(results, 
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))
print(summary(results))
#  Ordinary t-statistic
ordinary.t <- cfit$coef / cfit$stdev.unscaled / cfit$sigma
ordinary.t
#  Q-Q plots of t statistics
#  Points off the line may be differentially expressed
par(mfrow=c(1,2))
qqt(ordinary.t, df=cfit$df.residual, main="Ordinary t")
abline(0,1)
qqt(cfit$t, df=cfit$df.total,main="Moderated t")
abline(0,1)

#results table 
statsA = topTable(cfit, coef = "A", sort.by = "p", number = Inf) #These are your stats for contrast A
statsA
#i.e Adult.mock - Adult-virus. Do a head on the statsA object so you can see the table
statsA = statsA[statsA$adj.P.Val <= 0.05,] #FDR filter - Note, you may want to apply a FC cutoff in addtiion eg statsA = statsA[statsA$adj.P.Val <= 0.05 & abs(statsA$logFC)>=1.5,]
dim(statsA)
statsB = topTable(cfit, coef = "B", sort.by = "p", number = Inf)
dim(statsB)
statsB = statsB[statsB$adj.P.Val <= 0.05,] #FDR filter
dim(statsB)
statsC = topTable(cfit, coef = "C", sort.by = "p", number = Inf)
statsC = statsC[statsC$adj.P.Val <= 0.05,] #FDR filter
dim(statsC)

dev.off()
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")

par(mfrow=c(2,2))
plotMD(cfit,coef='C',status=results[,"C"]) #plot MD
volcanoplot(cfit,coef="C",highlight=100,names=cfit$genes$genes) #plot volcano you can do a separate plot for each coefficient

genes <- as.vector(rownames(statsC)[1:2]) #top 2 genes to plot

pheno <- v$targets$group

for(i in 1:length(genes)){
  boxplot(v$E[genes[i], ] ~ pheno, 
          las = 2, main = paste(statsC$symbol[i]))
  stripchart(v$E[genes[i], ] ~ pheno, 
             vertical = TRUE, method = "jitter", 
             pch = 21, col = nice.col, bg = "bisque", 
             add = TRUE)
}

RESDIR
write.csv(statsA, file.path(RESDIR, 'adults.csv'))
write.csv(statsB, file.path(RESDIR, 'kids.csv'))
write.csv(statsC, file.path(RESDIR, 'interaction.csv'))
rownames(cfit)

#4.0 Gene sets testing
##Ontology
###Here for the ontology testing you can substitute the contrast of interest so you can do it for contrast 'A' and 'B' and 'C'. You dont need to do all this ontology, maybe only camera competetive gene sets testing.
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)
converted <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                   filters = "hgnc_symbol",
                   values = rownames(cfit),
                   mart = mart)
converted
rownames(cfit)
m <- match(row.names(cfit), converted$hgnc_symbol)
table(is.na(m))
annot <- converted[m[!is.na(m)],]
table(is.na(annot))
annot = annot[match(rownames(cfit), annot$hgnc_symbol),]
annot
cfit$genes$entrez_id <- annot$entrezgene_id
GTF <- "../../gencode.v36.annotation.gtf"
gtf <- rtracklayer::readGFFAsGRanges(GTF)
gtf <- as.data.frame(gtf)
m <- match(rownames(cfit),gtf$gene_name)
table(is.na(m))
annot <- gtf[m[!is.na(m)],] 
dim(annot)
dim(cfit)
cfit
cfit$genes$start <- annot$start
cfit$genes$end <- annot$end
gene_length <- cfit$genes$end - cfit$genes$start
table(is.na(gene_length))

library(topGO)
goC <- goana(cfit, coef="C",species = "Hs", 
             geneid = cfit$genes$entrez_id,covariate=gene_length, FDR = 0.05)
goC
goCT <- topGO(goC, ontology="BP")
goCT
goCTu <- topGO(goC, sort = 'up', ontology = 'BP')
goCTu  
goCTd <- topGO(goC, sort = 'down', ontology = 'BP')
goCTd
goB <- goana(cfit, coef="B",species = "Hs", 
             geneid = cfit$genes$entrez_id,covariate=gene_length, FDR = 0.05)
goBT <- topGO(goB, ontology="BP")
goBT
goA <- goana(cfit, coef="A",species = "Hs", 
             geneid = cfit$genes$entrez_id,covariate=gene_length, FDR = 0.05)
goAT <- topGO(goA, ontology="BP")
goAT
statsC
source('../../R-script_folder/my_fgseaplot.R')
myfgsea(statsA, "A")
myfgsea(statsB, "B")
myfgsea(statsC, "C")

##Camera competetive gene sets test
# Load in the human hallmark gene sets
#load(file.path(DATADIR,'human_H_v5p2.rdata'))
#map entrez ids between gene sets and voom object
#H.ind <- ids2indices(Hs.H, annot$entrez)
#gst.camera <- camera(v,index=H.ind,design=design,contrast = cont.fit[,3],inter.gene.cor=0.01)
#gst.camera[1:10,]
#table(gst.camera$FDR < 0.05)
#write to file
#write.csv(gst.camera,file.path(RESDIR,"camera_baseline_cd3_int.csv"))
##Roast self contained test
# Load in the human hallmark gene sets
#load(file.path(DATADIR,'human_c7_v5p2.rdata'))
#map entrez ids between gene sets and voom object
#c7.ind <- ids2indices(Hs.c7, annot$entrez)
#tcell <- grep("CD4_",names(c7.ind))
#tcell.rst <- roast(v,index=c7.ind[tcell],design=design,contrast=cont.fit[,3],nrot=999)
#table(tcell.rst$FDR < 0.1)
#tcell.rst[1:7,]
#barcodeplot(cfit$t[,3], index=c7.ind[["GSE2770_IL12_ACT_VS_ACT_CD4_TCELL_6H_DN"]], main="T-statistic: IL12_ACT_VS_ACT_CD4_TCELL_6H_DN")
#write to file
#write.csv(tcell.rst[1:7,],file.path(RESDIR,"roast_baseline_cd3_int.csv"))
#R session info

#####topGO package analysis ####
require(topGO)
require(org.Hs.eg.db)
require(ggplot2)

####
genesA <- statsA$adj.P.Val
names(genesA) <- statsA$genes
head(genesA) #Just to double check
tail(genesA)
range(genesA)
genesB <- statsB$adj.P.Val
names(genesB) <- statsB$genes
genesC <- statsC$adj.P.Val
names(genesC) <- statsC$genes

#Adult.Virus - Adult.Mock
genes <- genesA
selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)

results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

require(ggplot2)
ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#Kid.Virus - Kid.Mock
genes <- genesB
selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)

results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#(Kid.Virus - Kid.Mock) - (Adult.Virus - Adult.Mock)
genes <- genesC
selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)

results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#https://www.biostars.org/p/350710/
sessionInfo()




