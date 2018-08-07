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

formatSignificance <- function(p){
  #sprintf(formatC(signif(p,digits=2), digits=2,format="fg", flag="#"))
  prettyNum(p,format="g", digits=2)
}

#formatSignificance(0.00312312)


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



normapprox <- function(va,vb){
  thep <- pnorm(0,mean(va)-mean(vb),sqrt(var(va)+var(vb)))
  return(min(thep,1-thep))
}

#####################################
# Calculate empirical p-values - either sb2 (standard permutation) or mw (mann-whitney) or empirical (original in the paper)
calc_empirical_p <- function(larger_group_values, smaller_group_values, usetest="mw") {   #"mw"
  nboot=10000
  if(length(smaller_group_values) > length(larger_group_values)){
    return(calc_empirical_p(smaller_group_values, larger_group_values))
  }
    
  #usetest <- "sb"
#  usetest <- "mw"

  if (usetest=="mw"){  #used in paper, mann-whitney
    
    # Alternative for calculating p-value using mann-whitney U test
    v <- wilcox.test(larger_group_values, smaller_group_values,alternative="less",paired=FALSE)
    total_p_value <- v$p.value
    return(min(total_p_value, 1-total_p_value))
    
  } else if (usetest=="sb"){  #straight bootstrap (boot+gaussian)


    if(all(smaller_group_values %in% larger_group_values)) {   
      ######## Here one group is a proper subgroup. 
      #Compare to other random subgroups of the same size, do not include these in the large group

      small_group_size = length(smaller_group_values) 
      large_group_size = length(larger_group_values) 

      va <- 1:nboot
      vb <- 1:nboot
      for (i in 1:nboot) {
    	  takesmall <- sample(1:large_group_size, small_group_size)
    	  va[i] <- median(larger_group_values[takesmall]) 
    	  vb[i] <- median(smaller_group_values)
      }
      total_p_value <- normapprox(va,vb)
      
      print(222)
      print(total_p_value)
      return(min(total_p_value, 1-total_p_value))
      
    } else {
      ######## One group is not a subset

      small_group_size = length(smaller_group_values) 
      large_group_size = length(larger_group_values) 
      
      va <- 1:nboot
      vb <- 1:nboot
      for (i in 1:nboot) {
        takesmall <- sample(1:small_group_size, small_group_size, replace = TRUE)
        takelarge <- sample(1:large_group_size, large_group_size, replace = TRUE)
        va[i] <- median(larger_group_values[takelarge])
        vb[i] <- median(smaller_group_values[takesmall])
      }
      total_p_value <- normapprox(va,vb)
      
      print(111)
      print(total_p_value)
      return(min(total_p_value, 1-total_p_value))

    }
  } else if (usetest=="sb2"){  ### Straight bootstrap, no normal distribution approximation
    
    
    if(all(smaller_group_values %in% larger_group_values)) {
      ######## Here one group is a proper subgroup. 
      #Compare to other random subgroups of the same size, do not include these in the large group
      
      small_group_size = length(smaller_group_values) 
      large_group_size = length(larger_group_values) 

      
      print("66666666666666")
      
      ad <- median(larger_group_values[larger_group_values %!in% smaller_group_values]) - median(smaller_group_values)
      print(ad)
      print()
      
      total_p_value=0
      for (i in 1:nboot) {
        takesmall <- sample(1:large_group_size, small_group_size)
        
        bootd <- median(larger_group_values[-takesmall]) - median(larger_group_values[takesmall])
#        print(bootd)
        
        if(ad<bootd){
        #if (median(larger_group_values[takesmall]) > median(smaller_group_values)) {
#        if (median(larger_group_values[takesmall]) > median(smaller_group_values)) {
          total_p_value=total_p_value+1
        }
      }
      total_p_value=total_p_value/nboot
      
      return(min(total_p_value, 1-total_p_value))
      
    } else {
      ######## One group is not a subset
      
      small_group_size = length(smaller_group_values) 
      large_group_size = length(larger_group_values) 
      
      total_p_value=0
      for (i in 1:nboot) {
        takesmall <- sample(1:small_group_size, small_group_size, replace = TRUE)
        takelarge <- sample(1:large_group_size, large_group_size, replace = TRUE)
        if (median(larger_group_values[takelarge]) > median(smaller_group_values[takesmall])) {
          total_p_value=total_p_value+1
        }
      }
      total_p_value=total_p_value/nboot
      
      return(min(total_p_value, 1-total_p_value))

    }
  }
}






#####################################
### "not in"
`%!in%` <- function(x,y) {
  !(x %in% y)
}



