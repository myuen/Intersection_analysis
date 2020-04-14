library(ggplot2)

# Function to make Volcano plot
makeVolcanoPlot <- function(df, pCutoff, lfcCutoff, pTitle) {
  
  minX <- round(min(df$logFC))
  maxX <- round(max(df$logFC))
  
  g <- ggplot(
    df, 
    aes(x = logFC, 
        y = -log10(adj.P.Val),
        color = focus_term,
        shape = focus_term)) + 
    geom_point(size = 0.9) + 
    
    scale_colour_manual(values = c("black", "grey")) +
    
    scale_x_continuous(breaks = c(minX, -2, 0, 2, maxX)) +
    
    scale_y_continuous(breaks = c(0, 1.3), labels = c(0, 0.05)) +
    
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    
    geom_abline(aes(intercept = -log10(pCutoff), slope = 0), colour = "blue") + 
    geom_vline(xintercept = lfcCutoff, colour = "blue") + 
    geom_vline(xintercept = -(lfcCutoff), colour = "blue") +
    
    labs(title = pTitle,
         x = expression(paste('Log'[2], " Fold Change")), 
         y = expression(paste('-Log'[10], " Adj. p-Value"))) +
    
    guides(colour = guide_legend("", keywidth = .5),
           shape = guide_legend("", keywidth = .5))
  
  return (g)
}
