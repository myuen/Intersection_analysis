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
    geom_point(size = 2.5) + 
    
    scale_colour_manual(values = c("black", "gray50")) +
    
    # Add break points on x-axis to show min, max, log2FC threshold
    scale_x_continuous(breaks = c(minX, -2, 0, 2, maxX)) +
    
    # Add break point to show adj p-value threshold
    # -log(0.05) = 1.3
    scale_y_continuous(breaks = c(0, 1.3), labels = c(0, 0.05)) +
    
    theme_classic() +
    
    # Define text size and remove grid lines
    theme(text = element_text(size = 16),
          legend.text = element_text(size = 16),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    
    # Add line to show log2FC and adj p-value threshold
    geom_abline(aes(intercept = -log10(pCutoff), slope = 0), colour = "blue") + 
    geom_vline(xintercept = lfcCutoff, colour = "blue") + 
    geom_vline(xintercept = -(lfcCutoff), colour = "blue") +
    
    # Add title for plot, x and y-axis
    labs(title = pTitle,
         x = expression(paste('Log'[2], 'Fold Change')), 
         y = expression(paste('-Log'[10], 'Adj. p-Value'))) +
    
    # Define legend
    guides(colour = guide_legend("", keywidth = 1.5),
           shape = guide_legend("", keywidth = 2.5))
  
  return (g)
}
