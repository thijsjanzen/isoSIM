print.individual <- function(x, ...) {
  print("Individual with two Chromosomes")
  v1 <- paste("Chromosome 1:",
              length(x$chromosome1) - 2,
              "junctions")
  v2 <- paste("Chromosome 2:", 
              length(x$chromosome2) - 2,
              "junctions")
  print(v1)
  print(v2)
}

print.population <- function(x, ...) {
  v1 <- paste("Population with",
              length(x),
              "individuals")
  print(v1)
}

plot.individual <- function(x, ...) {
  numColors <- 1
  alleles_chrom1 <- unique(x$chromosome1[, 2])
  alleles_chrom2 <- unique(x$chromosome2[, 2])
  num_colors <- length( unique( c(alleles_chrom1, alleles_chrom2)))
  color_palette <- grDevices::rainbow(num_colors)

  par(mfrow=c(2, 1))
  par(mar=c(2, 2, 2, 2))
  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in 1:length(x$chromosome1[, 1])) {
    xleft <- x$chromosome1[i, 1]
    xrght <- 1;
    if (i < length(x$chromosome1[, 1])) {
      xrght <- x$chromosome1[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome1[i, 2]
    colour_to_plot <- color_palette[colour_index]

    rect(xleft = xleft, 
         xright = xrght, 
         ybottom = 0, 
         ytop = 1, 
         col = colour_to_plot, 
         border = NULL)
  }

  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in 1:length(x$chromosome2[, 1])) {
    xleft <- x$chromosome2[i, 1]
    xrght <- 1;
    if (i < length(x$chromosome2[, 1])) {
      xrght <- x$chromosome2[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome2[i,2]
    colour_to_plot <- color_palette[colour_index]

    rect(xleft = xleft, 
         xright = xrght, 
         ybottom = 0, 
         ytop = 1, 
         col = colour_to_plot, 
         border = NULL)
  }
}

create_pop_class <- function(pop) {

  whole_pop <- list()
  cntr <- 1

  chrom1 <- c()
  chrom2 <- c();

  indic_chrom1 <- 1;
  indic_chrom2 <- 0;
  addIndiv = FALSE

  for (i in seq(from = 1,to = length(pop), by = 2)) {
    focal <- pop[c(i, i + 1)]
    if (focal[2] == -1) {
      if(indic_chrom2 == 0) {
        indic_chrom2 <- 1
        indic_chrom1 <- 0
      } else {
        addIndiv <- TRUE;
      }
    } 

    if (indic_chrom1 == 1) {
      chrom1 <- rbind(chrom1, focal)
    }
    if (indic_chrom2 == 1) {
      chrom2 <- rbind(chrom2, focal)
    }

    if (addIndiv == TRUE) {
      indiv <- list(chromosome1 = chrom1, 
                    chromosome2 = chrom2)

      class(indiv) <- "individual"

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
  return(whole_pop)
}

findtype <- function(chrom, pos) {
  
  if (length(chrom) == 2) {
    # if the chromosome has no junctions
    # the entire chromosome is the same type
    return(chrom[1, 2])
  }
  
  chromtype <- -1

  for (i in 2:length(chrom[, 1])) {
    
    if (chrom[i, 1] == pos) {
      chromtype <- chrom[i, 2]
      break;
    }

    if (chrom[i, 1] > pos) {
      chromtype <- chrom[i - 1, 2]
      break;
    }
  }

  if (chromtype < 0) {
    if (chrom[length(chrom[,1]), 1] < pos) {
      chromtype <- chrom[length(chrom[, 1]), 2]
    }
  }

  return(chromtype)
}