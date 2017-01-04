library(corrplot)
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
plot_corr_mx <- function(corr_mx){
  corrplot(corr_mx, method="color", col=col(200),  
           type="lower",# order="hclust", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation
           mar = c(4,4,4,4),
           # Combine with significance
           # p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           # hide correlation coefficient on the principal diagonal
           diag=FALSE 
  )
}
