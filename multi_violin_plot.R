
multi_violin_plot <- function(df){
  # Args:
  #   df: a data frame with each row as one group and each column as one feature
  melted <- melt(df)
  ggplot(melted, aes(x = variable, y = log10(value))) + 
    geom_violin(aes(fill = variable)) + coord_flip() +
    stat_summary(fun.y=mean, geom="point", size=2, color="red") + 
    theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1))
}

