create_random_markers <- function() {
  markers <- c();
  while (length(markers) < number_of_markers) {
    markers <- runif(number_of_markers, 0, 1)
    #if (length(which(markers == 0.0))) {
    if (sum(markers == 0.0)) {
      markers <- markers[- (markers == 0.0)] #remove borders
    }
    if (sum(markers == 1.0)) {
      markers <- markers[- (markers == 1.0)]
    }
    #remove duplicates
    if (length(which(duplicated(markers)))) {
      markers <- markers[-which(duplicated(markers))]
    }
  }
  markers <- sort(markers)
  return(markers)
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
  junct <- isoSIM::calculate_dist_junctions(pop)
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

  for (i in seq_along(indiv$chromosome2[, 1])) {
    left <- indiv$chromosome2[i, 1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome2[, 1])) {
      right <- indiv$chromosome2[i + 1, 1]
    }

    allele <- 1 + indiv$chromosome2[i, 2]
    alleles[allele] <- alleles[allele] + (right - left)
  }

  alleles <- alleles / sum(alleles)
  return(alleles)
}