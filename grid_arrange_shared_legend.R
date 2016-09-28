
library(ggplot2)
library(gridExtra)
library(grid)

# extract legend from a ggplot
# http://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot

gg_legend_extract <- function(a.gplot, legend.panel = FALSE){
  legend.position <- ifelse(legend.panel, 'left', 'bottom')
  g <- ggplotGrob(a.gplot + theme(legend.position = legend.position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  legend
}


# grid_arrange_shared_legend
# Share a legend between multiple plots using grid.arrange
# from http://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., grobs = NULL, nrow = 1 , ncol = 1, title = NULL) {
  if(is.null(grobs)){
    plots <- list(...)
  }else{
    plots <- grobs
  }
  # extract legend from the first ggplot
  legend <- gg_legend_extract(plots[[1]])
  
  lheight <- sum(legend$height)
  grid.arrange(do.call(arrangeGrob, c(grobs = lapply(plots, function(x) x + theme(legend.position="none")),
                       nrow = nrow, ncol = ncol)),
               legend,
               nrow = 2, ncol = 1,
               heights = unit.c(unit(0.9, "npc") - lheight, lheight),
               top = textGrob(title, gp=gpar(fontsize = 20)),
               padding = unit(2, "line"))
}

# use user supplied legend
grid_arrange_custom_legend <- function(..., grobs = NULL, legend, nrow = 1 , ncol = 1, title = NULL) {
  if(is.null(grobs)){
    plots <- list(...)
  }else{
    plots <- grobs
  }
 
  # use user supplied legend
  lheight <- sum(legend$height)
  grid.arrange(do.call(arrangeGrob, c(grobs = lapply(plots, function(x) x + theme(legend.position="none")),
                                      nrow = nrow, ncol = ncol)),
               legend,
               nrow = 2, ncol = 1,
               heights = unit.c(unit(0.9, "npc") - lheight, lheight),
               top = textGrob(title, gp=gpar(fontsize = 20)),
               padding = unit(2, "line"))
}