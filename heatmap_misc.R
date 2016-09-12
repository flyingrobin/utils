my_palette <- colorRampPalette(c("red", "white", "blue"))(256)
image(1-cors,  col = my_palette, zlim=c(0, 2))
filled.contour(z = 1-cors,  col= my_palette, nlevels = 256)

heatmap.2(mx, dendrogram = 'none', Rowv = hc_cpd$order, Colv = hc_ccl$order, trace = 'none',
          col = my_palette, cexRow = .5, cexCol = .5, margins = c(6, 8),
          symm=F, symkey=F, symbreaks=F)


# multi-panel heatmap
multi_heatmap <- function(graph.list, mfrow){
  # Args:
  #   graph.list: a list of graphs
  #   mfrow: panel layout as mfrow parameter in par()
  grid.newpage()
  lay <- grid.layout(nrow = mfrow[1], ncol = mfrow[2])
  pushViewport(viewport(layout = lay))
  for(i in 1:length(graph.list))
  grid.draw(editGrob(g, vp=viewport(layout.pos.row = 1,layout.pos.col = 1, clip=TRUE)))
  grid.draw(editGrob(g, vp=viewport(layout.pos.row = 1, layout.pos.col = 2, clip=TRUE)))
  upViewport(1)
}


heatmap(data[[1]])
g1 <- grid.grab()
heatmap.2(data.rand, dendrogram = 'none', Rowv = F, Colv = F, trace = 'none')
g2 <- grid.grab()


# Other plots
# Draw heatmap without clustering
# 
# heatmap.2(d_ccl, dendrogram = 'none', Rowv = F, Colv = F, 
#           trace = 'none', scale = 'row',
#           col = my_palette, labCol = "", labRow = "",
#           margins = c(3, 3))
# 
# heatmap.2(mx, 
#           symm=F, symkey=F, symbreaks=F, 
#           margins = c(5,5), trace = 'n', 
#           col = my_palette, labCol = "", labRow = "")
# # with clustering and dendrogram
# heatmap.2(d_ccl, trace = 'none',
#           col = my_palette, cexRow = .5, cexCol = .5, margins = c(3, 8))

library(pheatmap)
