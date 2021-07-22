#res <- statsCnf
#testgroup <- "C"
library(EnhancedVolcano)
myvolcano <- function(res, testgroup){
  keyvals <- ifelse(
    res$logFC< -1.5, 'royalblue',
    ifelse(res$logFC > 1.5, 'gold',
           'black'))
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'high'
  names(keyvals)[keyvals == 'black'] <- 'mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'low'
  
  vol2 <- EnhancedVolcano(res,
                          lab = NA, #rownames(res),
                          x = 'logFC',
                          y = 'adj.P.Val', 
                          selectLab = NULL, #rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                          #selectLab = c('OASL','CXCL10','ISG15', 'USP41','IFITM1','MX1'),
                          xlab = 'fold change',
                          title = NULL, #'Volcano plot - Difference in Kids compared to Adults ',
                          subtitle = NULL, #'~ without treatment effect(i.e. without Virus infection)',
                          pCutoff = 0.05,
                          FCcutoff = 1.5,
                          pointSize = 5,
                          labSize = 3.0,
                          labCol = 'black',
                          labFace = 'bold',
                          boxedLabels = F,
                          ylim = c(0,3),
                          xlim = c(-7.5,7.5),
                          colAlpha = 3/4,
                          gridlines.major = F,
                          gridlines.minor = F,
                          col = c("grey30", "palegreen3", "lightslateblue", "orangered2"),
                          legendPosition = 'none', #c(0.15, 0.15),
                          legendLabels = NA, #c("NS", 'log fold change (±1.5)', 'p-adj < 0.05', 
                                           #'log fold change (±1.5) & p-adj < 0.05'),
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
  pdf(file = paste0(RESDIR, testgroup, ".pdf"), width = 6, height = 5)
  print(vol2)
  dev.off()

}
