create_random_markers <- function(number_of_markers) {
  markers <- c();
  while (length(markers) < number_of_markers) {
    markers <- runif(number_of_markers, 0, 1)
    #remove duplicates
    which_dupl <- which(duplicated(markers))
    if (length(which_dupl)) {
      markers <- markers[-which_dupl]
    }
  }
  markers <- sort(markers)
  return(markers)
}

calculate_dist_junctions <- function(pop) {
  get_num_junctions <- function(indiv) {
    v1 <- length(indiv$chromosome1[, 1]) - 2
    v2 <- length(indiv$chromosome2[, 1]) - 2 #subract one for start
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

calculate_marker_frequency <- function(pop, location) {

  fun_chrom <- function(indiv) {
    return(c(findtype(indiv$chromosome1, location),
             findtype(indiv$chromosome2, location)))
  }

  types <- unlist(lapply(pop, fun_chrom))

  vv <- tibble::as.tibble(table(types))
  colnames(vv) <- c("ancestor", "frequency")
  vv$frequency <- vv$frequency / sum(vv$frequency)
  return(vv)
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