library(PCAtools)
library(dplyr)

ageclass <- "mocks"

directory <- "D:/17.SARSC-V2_experiment/"

setwd("D:/17.SARSC-V2_experiment/R-project/SCVO2_KVSA/")
resultsdir <- paste0(directory, "results")
testgroup <- "PCA_in_Kids_with_the_effect_of_treatment"

resdir <- "/salmon_res/"

dir.create(paste0(resultsdir))
resultsdir
dir.create(paste0( resultsdir, resdir))

resultsdir <- paste0(resultsdir, resdir)
resultsdir
dir.create(paste0(resultsdir, ageclass))
resultsdir <- paste0(resultsdir, ageclass)
resultsdir

MyreadCountMatrix <- read.table(paste0(resultsdir, '/', testgroup, "_hsapiens_normalized_counts.txt"), 
                              sep = '\t', header = T)

head(MyreadCountMatrix)

df <- (dplyr::select(MyreadCountMatrix, -X))
head(df)

dfmock <- df[,!c(TRUE, FALSE)]
head(dfmock)
dfage <- df[,11:20]
head(dfage)
metadata <- data.frame(row.names = colnames(df))
metadata$Group <- rep(NA, ncol(df))
metadata
metadata$Group[seq(1,20,2)] <- 'Virus'
metadata$Group[seq(2,20,2)] <- 'Mock'
metadata$CRP <- sample.int(100, size=ncol(df), replace=TRUE)
metadata$ESR <- sample.int(100, size=ncol(df), replace=TRUE)

metadata$age <- "Adult" 

metadata[11:20,]$age <- "Kid"
metadata

metaage <- metadata[metadata$age == "Adult",]
metaage

pca(dfage, metadata = metaage)
p <- pca(dfage, metadata = metaage)

resultsdir
pdf(file = paste0(resultsdir, '/', testgroup, '_PCA_plot.pdf'))
screeplot(p)
biplot(p)
biplot(p, colby = 'Group', shape = 'Group')
biplot(p, colby = 'Group', colkey = c(Mock = 'blue', Virus = 'red'),
       legendPosition = 'right', showLoadings = T, showLoadingsNames = T)

biplot(p, x = 'PC1', y = 'PC2', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC1', y = 'PC4', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC4', y = 'PC1', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC10', y = 'PC9', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC10', y = 'PC8', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC8', y = 'PC10', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC10', y = 'PC7', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

biplot(p, x = 'PC10', y = 'PC6', colby = 'Group', colkey = c(Mock='blue', Virus='red'),
       shape = 'Group', shapekey = c(Mock=10, Virus=21), legendPosition = 'right',
       pointSize = 3, encircle = T, ellipse = F, labSize = 3, labhjust = 0.7, labvjust = 0.5)

p$components

pairsplot(p, triangle = F, shape = 'Group', pointSize = 3,
          legendPosition = 'bottom', axisLabSize = 6, vline = c(0,0), hline = c(0,0),
          title = paste('PCA pairsplot for ', testgroup), titleLabSize = 12, colby = 'Group',
          colkey = c(Mock = 'blue', Virus = 'red'))

dfmock
project.pca <- prcomp(t(dfmock))

project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
barplot(project.pca.proportionvariances, cex.names=1, 
        xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), 
        ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

par(cex=1.0, cex.axis=0.8, cex.main=0.8)

pairs(project.pca$x[,1:5], col=c("blue", "red"), 
      main="Principal components analysis bi-plot\nPCs 1-5 : in pairs", pch=c(6, 16))
pairs(project.pca$x[,6:10], col=c("blue", "red"), 
      main="Principal components analysis bi-plot\nPCs 6-10 : in pairs", pch=c(6, 16))

dev.off()

