#library(devtools)
#install_github("thijsjanzen/isoSIM")
#library(isoSIM)


print.individual = function(x, ...) {
  print("Individual with two Chromosomes")
  v1 <- paste("Chromosome 1:", length(x$chromosome1)-2,"junctions")
  v2 <- paste("Chromosome 2:", length(x$chromosome2)-2,"junctions")
  print(v1)
  print(v2)
  
  invisible(x)
}

summary.individual = function(x, ...) {
  print("Individual with two Chromosomes")
  v1 <- paste("Chromosome 1:", length(x$chromosome1)-2,"junctions")
  v2 <- paste("Chromosome 2:", length(x$chromosome2)-2,"junctions")
  print(v1)
  print(v2)
  
  invisible(x)
}

plot.individual = function(x, ...) {
  
  par(mfrow=c(2,1)); par(mar=c(2,2,2,2))
  plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",bty="n");
  for(i in 1:length(x$chromosome1[,1])) {
    xleft <- x$chromosome1[i,1]
    xrght <- 1;
    if(i < length(x$chromosome1[,1])) xrght <- x$chromosome1[i+1,1]
    colourToPlot <- "red"
    if(x$chromosome1[i,2] == 1) colourToPlot <- "blue"
    
    rect(xleft = xleft, xright = xrght, ybottom = 0, ytop =1, col = colourToPlot, border = NULL)
  }
  
  plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",bty="n");
  for(i in 1:length(x$chromosome2[,1])) {
    xleft <- x$chromosome2[i,1]
    xrght <- 1;
    if(i < length(x$chromosome2[,1])) xrght <- x$chromosome2[i+1,1]
    colourToPlot <- "red"
    if(x$chromosome2[i,2] == 1) colourToPlot <- "blue"
    
    rect(xleft = xleft, xright = xrght, ybottom = 0, ytop =1, col = colourToPlot, border = NULL)
  }
}


#vx <- create_population(100, 2, 10, 1, 42, FALSE)

#pop <- vx$population

create_pop_class <- function(pop) {
  
  whole_pop <- list()
  cntr <- 1
  
  chrom1 <- c()
  chrom2 <- c();
  
  indic_chrom1 <- 1;
  indic_chrom2 <- 0;
  addIndiv = FALSE
  
  for(i in seq(from = 1,to=length(pop), by = 2)) {
    focal <- pop[c(i,i+1)]
    if(focal[2] == -1) {
      if(indic_chrom2 == 0) {
        indic_chrom2 <- 1
        indic_chrom1 <- 0
      } else {
        addIndiv <- TRUE;
      }
    } else {
      if(indic_chrom1 == 1) {
        chrom1 <- rbind(chrom1, focal)
      }
      if(indic_chrom2 == 1) {
        chrom2 <- rbind(chrom2, focal)
      }
    }
      
    if(addIndiv == TRUE) {
      indiv <- list(chromosome1 = chrom1, 
                    chromosome2 = chrom2)
      
      class(indiv) <- "individual"
      #print(indiv);
      
      whole_pop[[cntr]] <- indiv
      cntr <- cntr + 1
      
      indic_chrom2 <- 0
      indic_chrom1 <- 1
      addIndiv <- FALSE
      chrom1 <- c()
      chrom2 <- c()
    }
    
  }
  class(whole_pop) <- "population"
}