library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)
library(gridExtra)
multi_gl <- function(gl.list, nrow, ncol, title = NULL){
  grid.newpage()
  # setup layout
  nrow <- ifelse(is.null(title), nrow, nrow + 1)
  gl <- grid.layout(nrow = nrow, ncol = ncol)
  nrow_start = 1

  # setup viewports
  for(i in 1:nrow){
    for(j in 1:ncol){
      assign(paste('vp', i, j, sep = ''), viewport(layout.pos.col = j, layout.pos.row = i))
    }
  }
  
  # plotting 
  pushViewport(viewport(layout=gl))
  if(!is.null(title)){
    grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol),
              gp=gpar(fontsize=20, fontfamily = "sans"))
    nrow_start <- 2
  }
  
  k <- 1
  for(i in nrow_start:nrow){
    for(j in 1:ncol){
      if(k > length(gl.list)) break # alreadt drew all plots
      
      vp <- get(paste('vp', i, j, sep = ''))
      pushViewport(vp) # start new plot
      
      # if the plot is the initial plot
      if(i == 1 & j == 1) par(new=TRUE, fig=gridFIG()) 
      
      grid.draw(gl.list[[k]]) # draw plot
      popViewport() # done with this viewport
      k <- k+1
    }
  }
}
