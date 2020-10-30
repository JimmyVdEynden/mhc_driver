# Format input
# A tibble: 6 x 8
# cancer_type    OR        p `2.5 %` `97.5 %` ci_low ci_high             q
# <chr>       <dbl>    <dbl>   <dbl>    <dbl>  <dbl>   <dbl>         <dbl>
#   1 BLCA         1.15 3.32e- 4  0.0639    0.218  1.07     1.24 0.000921     
# 2 BRCA         1.23 6.93e-10  0.142     0.274  1.15     1.32 0.00000000578
# 3 CESC         1.09 1.10e- 1 -0.0196    0.192  0.981    1.21 0.125        
# 4 COAD         1.20 3.13e- 8  0.116     0.243  1.12     1.27 0.000000156  
# 5 ESCA         1.18 2.87e- 2  0.0176    0.321  1.02     1.38 0.0399       
# 6 GBM          1.15 3.99e- 3  0.0456    0.240  1.05     1.27 0.00636   

plot_OR<- function(results,title="title") {
  plot_df = results %>%
    arrange(desc(OR))
  
  ggplot(plot_df, aes(x = OR, y = cancer_type)) + # the x-position of the dot is given by the column OR
    theme_bw() + # theme without background
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # rotate the x-labels so they are better visible
    
    geom_point(size = 1.5) + # size of the dot that indicates the odds ratio
    geom_errorbarh(aes(xmax = ci_high, xmin = ci_low), size = .5, height = .2) + # add an error bar for the confidence interval
    
    geom_vline(xintercept = 1, linetype = "dashed") +
    
    scale_y_discrete(limits = plot_df$cancer_type) + # where to put breaks on the y-axis
    scale_x_continuous(breaks = scales::breaks_width(0.2)) +
    
    xlab("Odds ratio") + # label for the x-axis
    ylab("") + # label for the y-axis
    ggtitle(label = title)
}

theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

