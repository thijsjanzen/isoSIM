calc_heterozygosity_indiv <- function(indiv) {

  if (length(indiv$chromosome1) > 100) {
    #now we have to get the genetic type at each stretch
    heterozygosity <- 0
    chrom1 <- 2
    chrom2 <- 2
    left <- 0
    right <- 1
    while (chrom1 < length(indiv$chromosome1[, 1]) &&
           chrom2 < length(indiv$chromosome2[, 1]) ) {
      is_hetero <- 0
      if (indiv$chromosome1[chrom1 - 1, 2] !=
         indiv$chromosome2[chrom2 - 1, 2]) {
         is_hetero <- 1
      }

      if (indiv$chromosome1[chrom1, 1] < indiv$chromosome2[chrom2, 1]) {
         right  <- indiv$chromosome1[chrom1, 1]
         chrom1 <- chrom1 + 1
      } else {
         right  <- indiv$chromosome2[chrom2, 1]
         chrom2 <- chrom2 + 1
      }

      heterozygosity <- heterozygosity + (right - left) * is_hetero
      left <- right
    }

    return(heterozygosity)
  } else {
    pos <- unique(c(0,
                    indiv$chromosome1[, 1],
                    indiv$chromosome2[, 1], 1))

    pos <- sort(pos)
    left <- 0
    right <- pos[1]
    heterozygosity <- 0
    for (i in 2:length(pos)) {
      left <- right
      right <- pos[i]
      type1 <- isoSIM::findtype(indiv$chromosome1, left)
      type2 <- isoSIM::findtype(indiv$chromosome2, left)
      if (type1 != type2) {
        heterozygosity <- heterozygosity + (right - left)
      }
    }
    return(heterozygosity)
  }
}

calculate_dist_junctions <- function(pop) {
  get_num_junctions <- function(indiv) {
    v1 <- length(indiv$chromosome1[, 1]) - 1
    v2 <- length(indiv$chromosome2[, 1]) - 1 #subract one for start
    return(c(v1, v2))
  }

  vx <- as.numeric(sapply(pop, get_num_junctions))

  return(vx)
}

plot_dist_junctions <- function(pop) {
  junct <- calculate_dist_junctions(pop)
  vx <- table(junct)
  barplot(vx)
}

calc_allele_frequencies <- function(indiv, alleles) {

  for (i in seq_along(indiv$chromosome1[, 1])) {
    left <- indiv$chromosome1[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome1[, 1])) {
      right <- indiv$chromosome1[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome1[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  for (i in seq_along(indiv$chromosome2[,1])) {
    left <- indiv$chromosome2[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome2[,1])) {
      right <- indiv$chromosome2[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome2[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  alleles <- alleles / sum(alleles)
  return(alleles)
}