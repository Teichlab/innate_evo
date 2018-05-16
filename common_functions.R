#####################################
# Formatting of significance values
formatSignificance <- function(p){
  if(p<0.001){
    "***"
  } else if(p<0.01){
    "**"
  } else if(p<0.05) {
    "*"
  } else if(p>0.9) {
    ""
  } else {
    sprintf("%.4s",p)
  }
}


#####################################
# multi-plots
# accessory function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#####################################
# Calculate empirical p-values
calc_empirical_p <- function(larger_group_values, smaller_group_values) {
  nboot=10000
  small_group_size = length(smaller_group_values) 
  large_group_size = length(larger_group_values) 
  if(small_group_size>large_group_size){
    calc_empirical_p(smaller_group_values, larger_group_values)
  } else {
    #subsample_larger_group_values = as.matrix(sample(larger_group_values, small_group_size, replace = FALSE))
    total_p_value=0
    for (i in 1:nboot) {
      subsample_larger_group_values = sample(larger_group_values, small_group_size, replace = FALSE)
      if (median(subsample_larger_group_values) > median(smaller_group_values)) {
        total_p_value=total_p_value+1
      }
    }
    total_p_value=total_p_value/nboot
    min(total_p_value, 1-total_p_value)
  }
}



`%!in%` <- function(x,y) !(x %in% y)

