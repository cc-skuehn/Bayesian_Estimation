### Helper function for plotting

# dev.off() resets par and other device-specific settings!

addTrans <- function(color,trans)
  {
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

mult_plot <- function(xdat,nofvar=2,transp=100,xdatclust="1",marg=2)
  {
  
  # Multiplot configuration
  par(mfrow=c(nofvar,nofvar))
  # Small margins on all sides
  par(mar=rep(marg,4))
  #par(new=TRUE)
  x1=xdatclust==1
  x2=xdatclust==2
  x3=xdatclust==3
  xdatclust2=xdatclust
  xdatclust2[x1]="blue"
  xdatclust2[x2]="black"
  xdatclust2[x3]="red"
#  xdatclust2[x1]="black"
#  xdatclust2[x2]="red"
#  xdatclust2[x3]="blue"
  for (i in 1:nofvar){
    for (j in 1:nofvar){
      # Plot core points
      par(pch=20, cex=.5)
      #par(pch=".")

      if (i==j){
        #plot(xdat[,c(i,j)],col=addTrans(xdatclust2,transp),xlim=c(0,1),ylim=c(0,1))
        #hist(xdat[,i],col=addTrans(xdatclust2,transp),xlim=c(0,1),100)
        hist(xdat[,i],col=addTrans(xdatclust2,transp),xlim=c(min(xdat[,i]),max(xdat[,i])),100)
        mtext(toString(colnames(xdat)[i]),side=3,line=-3.0, 
              at=par("usr")[1]+0.45*diff(par("usr")[1:2]),
              cex=0.75)
      } else {
        plot(xdat[,c(i,j)],col=addTrans(xdatclust2,transp))
        #plot(xdat[,c(i,j)],col=addTrans(xdatclust2,transp),xlim=c(0,1),ylim=c(0,1))
        } # end if
    } # end j
  } # end i
  #par(new=TRUE)
} # end function

# scaling functions
linMap <- function(x, from=0, to=1)
  {
  if (min(x)==max(x)) {
    print('Min=Max, set everything to 1, no scaling')
    return (x<-rep(1,length(x)))
  } else {
    (x - min(x)) / max(x - min(x)) * (to - from) + from    
  }
} # end linMap

scale_array <- function(xarr,from=0,to=1)
  {
  xarr_un<-as.matrix(xarr)
  noff <- length(xarr_un[1,])
  nofp <- length(xarr_un[,1])
  xarr_sc <- matrix(data=NA,nrow=nofp,ncol=noff)
  for (i in 1:noff){
    xarr_sc[,i] = linMap(xarr_un[,i])
  }
  return(xarr_sc)
} # end scale_array

# Dummy multiplot for variables via stringnames
myplot = function(df, x_string, y_string) 
  {
  ggplot(df, aes_string(x = x_string, y = y_string)) + geom_point()
} # end myplot

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
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
} # end ggplot multplot

# plot regression infos
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

# for nicer ticks
number_ticks <- function(n) {function(limits) pretty(limits, n)}
