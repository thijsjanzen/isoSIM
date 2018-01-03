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
  alleles_chrom1 <- unique(x$chromosome1[, 2])
  alleles_chrom2 <- unique(x$chromosome2[, 2])
  num_colors <- 1 + max(alleles_chrom1, alleles_chrom2)    #length(unique(c(alleles_chrom1, alleles_chrom2)))
  color_palette <- grDevices::rainbow(num_colors)

  par(mfrow = c(2, 1))
  par(mar = c(2, 2, 2, 2))
  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in seq_along(x$chromosome1[, 1])) {
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
         border = NA)
  }

  plot(NA,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")

  for (i in seq_along(x$chromosome2[, 1])) {
    xleft <- x$chromosome2[i, 1]
    xrght <- 1
    if (i < length(x$chromosome2[, 1])) {
      xrght <- x$chromosome2[i + 1, 1]
    }
    colour_index <- 1 + x$chromosome2[i, 2]
    colour_to_plot <- color_palette[colour_index]

    rect(xleft = xleft,
         xright = xrght,
         ybottom = 0,
         ytop = 1,
         col = colour_to_plot,
         border = NA)
  }
}

plot_chromosome <- function(chrom, xmin, xmax) {
  alleles <- unique(chrom[, 2])
  num_colors <- length(unique(alleles))
  color_palette <- grDevices::rainbow(num_colors)
  
  par(mfrow = c(1, 1))
  par(mar = c(2, 2, 2, 2))
  plot(NA,
       xlim = c(xmin, xmax),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       bty  = "n")
  
  for (i in seq_along(chrom[, 1])) {
    xleft <- chrom[i, 1]
    xrght <- 1;
    if (i < length(chrom[, 1])) {
      xrght <- chrom[i + 1, 1]
    }
    colour_index <- 1 + chrom[i, 2]
    colour_to_plot <- color_palette[colour_index]
    
    rect(xleft = xleft,
         xright = xrght,
         ybottom = 0,
         ytop = 1,
         col = colour_to_plot,
         border = NA)
  }
}



create_pop_class <- function(pop) {
  whole_pop <- list()
  cntr <- 1
  chrom1 <- c()
  chrom2 <- c()
  indic_chrom <- 1
  add_indiv <- FALSE

  for (i in seq(from = 1, to = length(pop), by = 2)) {
    focal <- pop[c(i, i + 1)]
    if (indic_chrom == 1) {
      chrom1 <- rbind(chrom1, focal)
    }

    if (indic_chrom == 2) {
      chrom2 <- rbind(chrom2, focal)
    }

    if (focal[2] == -1) {
      if (indic_chrom == 1) {
        indic_chrom <- 2
      } else {
        add_indiv <- TRUE;
      }
    }

    if (add_indiv == TRUE) {
      indiv <- list(chromosome1 = chrom1,
                    chromosome2 = chrom2)

      class(indiv) <- "individual"

      whole_pop[[cntr]] <- indiv
      cntr <- cntr + 1

      indic_chrom <- 1
      add_indiv <- FALSE
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

  a <- which(chrom[, 1] == pos)
  if (length(a) > 0) {
    chromtype <- chrom[a[1], 2]
  } else {
    b <- which(chrom[, 1] > pos)
    chromtype <- chrom[b[1] - 1, 2]
  }

  if (chromtype[[1]] < 0) {
    if (chrom[length(chrom[, 1]), 1] < pos) {
      chromtype <- chrom[length(chrom[, 1]), 2]
    }
  }

  return(chromtype[[1]])
}