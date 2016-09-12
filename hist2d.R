library(MASS)
library(RColorBrewer)

hist2d <- function(X, Y, n = 25){
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)

  df = data.frame(x = X, y = Y, n = n)
  df = na.omit(df)

  h1 <- hist(df$x, breaks=n, plot=F)
  h2 <- hist(df$y, breaks=n, plot=F)
  top <- max(h1$counts, h2$counts)
  k <- kde2d(df$x, df$y, n=n)

  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
  image(k, col=r) #plot the image
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)

}

# An even nicer plot
gghist2d <- function(X, Y){
  library("ggExtra")
  df <- data.frame(X, Y)
  sp2 <- ggplot(df,aes(X, Y)) + geom_point()
  # Marginal histogram plot
  ggMarginal(sp2 + theme_gray(), type = "histogram",
             fill = "steelblue", col = "darkblue")
}


# hist2d with factor groups,
# go to http://rforpublichealth.blogspot.hk/2014/02/ggplot2-cheatsheet-for-visualizing.html
