{
  library(VennDiagram)
  vd.param <- list(filename = NULL, height = 100, width = 100, scale = T, 
                   cat.fontface = 2, cat.pos = c(-10, 10), cat.dist = c(0.02, 0.05), margin = .08,
                   main.pos = c(0.5, 0.95), main.cex = 2, cat.fontfamily = "sans", cat.cex = 1.5, cex = 2,
                   main.fontface = 2, main.fontfamily = 'sans', fontfamily = "sans")
  
  vg.ttest.ccl <- lapply(1:8, function(i) {
    hotblock_ttest <- get(paste('hotblock_ttest', i, sep=''))
    res.ccl <- res.purity.ccl.list[[i]]
    #res.cpd <- get(paste('res.purity.cpd', i, sep=''))
    
    hps <- hotblock_ttest %>% filter(qval < 0.05)
    
    vd.ccl.data <- list(x = list(hotblock = hps$col.node %>% unique() %>% as.character(),
                                 feat_enriched = res.ccl$node %>% unique() %>% na.omit() %>% as.character()),
                        category.names = c('Nodes in hot block', 'Nodes with \nenriched feature'),
                        main = data.nms[i], fill = c('red', 'green'))
    # vd.cpd.data <- list(x = list(hotblock = hps$row.node %>% unique() %>% as.character(),
    #                              feat_enriched = res.cpd$node %>% unique() %>% as.character()),
    #                     category.names = c('Nodes in hot block', 'Nodes with enriched feature'),
    #                     main = data.nms[i])
    # 
    do.call(venn.diagram, c(vd.param, vd.ccl.data))
    #vd.cpd <- do.call(venn.diagram, c(vd.param, vd.cpd.data))
    
  })
  vg.ttest.cpd <- lapply(1:8, function(i) {
    hotblock_ttest <- get(paste('hotblock_ttest', i, sep=''))
    res.cpd <- res.purity.cpd.list[[i]]
    
    hps <- hotblock_ttest %>% filter(qval < 0.05)
    vd.cpd.data <- list(x = list(hotblock = hps$row.node %>% unique() %>% as.character(),
                                 feat_enriched = res.cpd$node %>% unique() %>% as.character()),
                        category.names = c('Nodes in \nhot block', 'Nodes with \nenriched feature'),
                        main = data.nms[i], fill = c('red', 'blue'))
    
    v <- do.call(venn.diagram, c(vd.param, vd.cpd.data))
    
  })
  grid.newpage()
  grid.draw(v)
  # pdf('venn_nodes_in_hps_or_w_feat_scaled.pdf', height = 8, width = 12)
  # multi_gl(vg.ttest.ccl, 3, 3)
  # multi_gl(vg.ttest.cpd, 3, 3)
  # dev.off()
}