calculate_heterozygosity_and_freq_table <- function(pop,
                                                 number_of_founders) {
  pop_for_cpp <- c()
  for (i in seq_along(pop)) {
    x <- pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }

  results <- calculate_summaryStats(pop_for_cpp,
                                    number_of_founders)
  return(results)
}

calculate_heterozygosity <- function(pop) {

  #first we have to unwind all the individuals into one large vector
  pop_for_cpp <- c()
  for (i in seq_along(pop)) {
     x <- pop[[i]]$chromosome1
     chrom1 <- as.vector(t(x))
     x <- pop[[i]]$chromosome2
     chrom2 <- as.vector(t(x))
     pop_for_cpp <- c(pop_for_cpp,
                      chrom1,
                      chrom2)
  }

  heterozygosity <- isoSIM::calc_heterozygosity_cpp(pop_for_cpp)
  return(heterozygosity)
}

create_loci_matrix <- function(pop1,
                               pop2,
                               number_of_founders,
                               number_of_markers,
                               random_markers) {

  all_loci <- matrix(nrow = length(pop1) + length(pop2),
                     ncol = 1 + number_of_markers, 0)
  all_loci[, 1] <- c(rep(1, length(pop1)), rep(2, length(pop1)))
  colnames(all_loci) <- c("population", 1:number_of_markers)

  markers <- seq(1e-9, 1 - (1e-9), length.out = number_of_markers)
  if (random_markers) {
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
  }

  for (x in seq_along(markers)) {
    focal_marker <- markers[x]
    for (i in seq_along(pop1)) {
      allele_1 <- 10 + isoSIM::findtype(pop1[[i]]$chromosome1, focal_marker)
      allele_2 <- 10 + isoSIM::findtype(pop1[[i]]$chromosome2, focal_marker)
      final_allele <- paste0(allele_1, allele_2)
      all_loci[i, x + 1] <- as.numeric(final_allele)
    }

    for (i in seq_along(pop2)) {
      allele_1 <- 10 + isoSIM::findtype(pop2[[i]]$chromosome1, focal_marker)
      allele_2 <- 10 + isoSIM::findtype(pop2[[i]]$chromosome2, focal_marker)
      final_allele <- paste0(allele_1, allele_2)
      all_loci[length(pop1) + i, x + 1] <- as.numeric(final_allele)
    }
  }

  return(all_loci)
}

hierfstat_basic_stats <- function(pop1,
                              pop2,
                              number_of_founders,
                              number_of_markers = 100,
                              random_markers = FALSE) {

  number_of_markers <- round(number_of_markers)

  all_loci <- create_loci_matrix(pop1, pop2,
                                 number_of_founders, number_of_markers,
                                 random_markers)

  hierf_sum_overall <- hierfstat::basic.stats(as.data.frame(all_loci))$overall

  return(hierf_sum_overall)
}


hierfstat_fst_wc <- function(pop1,
                             pop2,
                             number_of_founders,
                             sampled_individuals,
                             number_of_markers = 100,
                             random_markers = FALSE) {

  number_of_markers <- round(number_of_markers)

  all_loci <- create_loci_matrix(pop1[sample(1:length(pop1), sampled_individuals)], 
                                 pop2[sample(1:length(pop2), sampled_individuals)],
                                 number_of_founders, number_of_markers,
                                 random_markers)

  hierf_wc <- hierfstat::wc(as.data.frame(all_loci))

  return(hierf_wc$FST)
}