myfgsea <- function(stat, contrast) {
  stat <- stats %>% 
    dplyr::select(genes, logFC) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(genes) %>% 
    summarize(stat=mean(logFC))
  
  stat
  ranks <- tibble::deframe(stat)
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
  
  t <- fgseaResTidy[fgseaResTidy$padj <= 0.05,]
  t$pathway
  t$leadingEdge[[17]]
  t$leadingEdge[[18]]
  t$leadingEdge[[23]]
  t
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
  pdf(file = paste0(RESDIR, contrast), width = 6, height = 5)
  print(p1)
  dev.off()
  
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
  
  #print(p1)
  
}
