library(corrplot)
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram

plot_corr_mx <- function(corr_mx){
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  
  corrplot(corr_mx, method="color", col=col(200),  
           type="lower",order="hclust", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation
           mar = c(1,1,1,1),
           tl.cex = .6, number.cex = 0.6,
           # Combine with significance
           # p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           # hide correlation coefficient on the principal diagonal
           diag=T 
  )
}
