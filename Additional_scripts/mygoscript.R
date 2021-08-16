my_go_script <- function(resG, testgroup, fdr.threshold) {
  graphics.off()
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
  EG2KEGG<- getBM(mart = sp, values = gene.vector, attributes = c("hgnc_symbol", "go_id"))
  head(EG2KEGG)
  
  geneID2GO <- by(EG2KEGG$hgnc_symbol,
                  EG2KEGG$go_id,
                  function(x) as.data.frame(x))
  head(geneID2GO,15)
  class(geneID2GO)
  write.csv(as.data.frame(EG2KEGG), file = paste0(resultsdir, "/", testgroup, "_geneID.csv"))
  
  gene.vector
  
  lengthdata <- getlength(genes = names(gene.vector), "hg19",  "geneSymbol")
  (names(lengthdata) <- names(gene.vector))
  
  pwf<- nullp(gene.vector,id = geneID2GO,bias.data = lengthdata)
  plotPWF(pwf,binsize = 1000)
  
  goResults <- goseq(pwf,gene2cat = geneID2GO, method = "Wallenius", repcnt = 2000)
  
  head(goResults)
  nrow(goResults)
  
  write.csv( as.data.frame(goResults), file=paste0(resultsdir, "/", testgroup, "_gores.csv" ))
  enriched.GO=goResults$category[goResults$over_represented_pvalue<0.05]
  
  class(enriched.GO)
  head(enriched.GO)
  length(enriched.GO)
  write.csv(as.character(enriched.GO), file = paste0(resultsdir, "/", testgroup, "enriched_go_", myspecies, ".csv"))
  
  capture.output(for(go in enriched.GO[1:20]) { print(GOTERM[[go]])
    cat("--------------------------------------\n")
  }
  , file=paste0(resultsdir, "/", testgroup,"_SigGo.txt"))
  
  goResults[goResults$ontology == 'BP',]
  goBPRes <- goResults[goResults$ontology == 'BP',] 
  goBPRes$over_represented_pvalue
  goBPRes_up <- goBPRes[which(goBPRes$over_represented_pvalue < 0.05),]
  goBPRes_down <- goBPRes[which(goBPRes$under_represented_pvalue < 0.05),]
  goBPRes_down
  dim(goBPRes_up)
  ####################################### ~ visualize the top 30 hits*****************************##########################
  graphics.off()
  pdf(file = paste0(resultsdir, "/", testgroup,"_goseq_enrichment.pdf"))
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
  
  goBPRes_up %>% 
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
  paste0(resultsdir, "/", testgroup,"_goBPseq_enrichment.pdf")
  graphics.off()
  testgroup
  pdf(file = paste0(resultsdir, "/", testgroup,"_goBPseq_enrichment.pdf"), width = 7, height = 6)
  fig <- goBPRes_up %>% 
    top_n(75, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, y=term)) + 
    geom_bar(stat = "identity", fill = "mediumpurple3") +
    theme(text = element_text(size = 10, family = "serif")) +
    labs(x="Hits (%)", y=NULL, colour="p value", size="Count") +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(fig)
  dev.off()
  
  #write.csv(as.data.frame(goBPRes), file = paste0(resultsdir,"/", testgroup, "_go_bp.csv"))
  
  ####################################### ~ topGO mapping ~ ############################################
  bg_genes <- row.names(ddsKEEP)
  length(bg_genes)
  de.genes
  candidate_list <- de.genes
  length(candidate_list)
  db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
  go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genes, mart=db)
  head(go_ids)
  
  # build the gene 2 GO annotation list (needed to create topGO object)
  gene_2_GO=unstack(go_ids[,c(1,2)])
  
  # remove any candidate genes without GO annotation
  keep = candidate_list %in% go_ids[,2]
  keep =which(keep==TRUE)
  candidate_list=candidate_list[keep]
  length(candidate_list)
  # make named factor showing which genes are of interest
  library(topGO)
  geneList=factor(as.integer(bg_genes %in% candidate_list))
  names(geneList)= bg_genes
  GOdata=new('topGOdata', ontology='BP', 
             allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
  
  # define test using the classic algorithm with fisher (refer to [1] if you want to understand how the different algorithms work)
  #classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
  
  # define test using the weight01 algorithm (default) with fisher
  weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher') 
  
  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO=usedGO(GOdata)
  options(scipen = 999)
  all_res=GenTable(GOdata, weightFisher=weight_fisher_result, 
                   orderBy='weightFisher', topNodes=length(allGO))
  
  #performing BH correction on our p values
  p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
  p.adj
  
  # create the file with all the statistics from GO analysis
  all_res_final=cbind(all_res,p.adj)
  all_res_final=all_res_final[order(all_res_final$p.adj),]
  
  all_res_final
  #get list of significant GO before multiple testing correction
  keep <- which(all_res_final$weightFisher<=0.001)
  keep
  all_res_final[keep,]
  
  results.table.p <- as.data.frame(all_res_final[keep,])
  dim(results.table.p)
  results.table.p
  #get list of significant GO after multiple testing correction
  head(all_res_final$p.adj)
  results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
  results.table.bh
  
  #save first top 50 ontolgies sorted by adjusted pvalues
  write.table(all_res_final[1:50,],paste0(resultsdir, "/", testgroup, "summary_topGO_analysis.csv"),
              sep=",",quote=FALSE,row.names=FALSE)
  showSigOfNodes(GOdata, score(weight_fisher_result), 
                 useInfo = "all", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
  
  myterms = results.table.p$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
  mygenes = genesInTerm(GOdata, myterms)
  mygenes
  var=c()
  for (i in 1:length(myterms))
  {
    myterm=myterms[i]
    mygenesforterm= mygenes[myterm][[1]]
    mygenesforterm=paste(mygenesforterm, collapse=',')
    var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
  }
  
  write.table(var, paste0(resultsdir,"/", testgroup, "_genetoGOmapping.txt"),sep="\t",quote=F)
}
