library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)
library(gridExtra)
multi_gl <- function(gl.list, nrow, ncol){
  grid.newpage()
  # setup layout
  gl <- grid.layout(nrow = nrow, ncol = ncol)
  
  # setup viewports
  for(i in 1:nrow){
    for(j in 1:ncol){
      assign(paste('vp', i, j, sep = ''), viewport(layout.pos.col = j, layout.pos.row = i))
    }
  }
  
  # plotting 
  pushViewport(viewport(layout=gl))
  k <- 1
  for(i in 1:nrow){
    for(j in 1:ncol){
      if(k > length(gl.list)) break # alread drew all plots
      
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
