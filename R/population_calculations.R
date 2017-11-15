
calc_heterozygosity <- function(indiv) {
  
  #first, get a list of all junction positions
  pos <- unique(c(0,
                  indiv$chromosome1[, 1], 
                  indiv$chromosome2[, 1],
                  1))
  
  pos <- sort(pos)
  
  #now we have to get the genetic type at each stretch
  
  left <- 0
  right <- pos[1]
  heterozygosity <- 0
  for (i in 2:length(pos)) {
    left <- right;
    right <- pos[i]
    type1 <- findtype(indiv$chromosome1, left)
    type2 <- findtype(indiv$chromosome2, left)
    
    if (type1 != type2) {
      heterozygosity <- heterozygosity + (right - left)
    }
  }
  
  return(heterozygosity)
}

calculate_pop_heterozygosity <- function(pop) {
  a <- sapply(pop, calc_heterozygosity)
  return( mean(a))
}

calculate_dist_junctions <- function(pop) {
  get_num_junctions <- function(indiv) {
    v1 <- length(indiv$chromosome1[, 1]) - 1
    v2 <- length(indiv$chromosome2[, 1]) - 1 #subract one for start
    return( c(v1, v2))
  }
  
  vx <- as.numeric(sapply(pop,get_num_junctions))
  
  return(vx)
}

plot_dist_junctions <- function(pop) {
  junct <- calculate_dist_junctions(pop)
  vx <- table(junct)
  barplot(vx)
}

calc_allele_frequencies <- function(indiv) {
  alleles <- rep(0,1+max(indiv$chromosome1[,2],indiv$chromosome2[,2]))
  
  for (i in 1:length(indiv$chromosome1[,1])) {
    left <- indiv$chromosome1[i,1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome1[,1])) right <- indiv$chromosome1[i+1,1]
    
    allele <- 1 + indiv$chromosome1[i,2]
    alleles[allele] <- alleles[allele] + (right - left)
  }
  
  for (i in 1:length(indiv$chromosome2[,1])) {
    left <- indiv$chromosome2[i,1]
    right <- 1
    if (i + 1 <= length(indiv$chromosome2[,1])) right <- indiv$chromosome2[i+1,1]
    
    allele <- 1 + indiv$chromosome2[i,2]
    alleles[allele] <- alleles[allele] + (right - left)
  }
  
  alleles <- alleles / sum(alleles)
  
  
  return(alleles)
}





calculate_fst <- function(pop1, 
                          pop2) {
  
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 100
  morgan <- 1
  overlap <- 0.0
  write_to_file <- FALSE
  
  vx <- create_two_full_populations(pop_size, number_of_founders, 
                                    run_time, morgan, 42, 
                                    overlap, write_to_file)
  
  pop1 <- vx$Population_1
  pop2 <- vx$Population_2
  
  h_a <- calculate_pop_heterozygosity(pop1)
  h_b <- calculate_pop_heterozygosity(pop2)
  
  combined_pop <- c(pop1, pop2)
  class(combined_pop) <- "population"
  h_c <- calculate_pop_heterozygosity(combined_pop)
  
  fst <- (h_c - h_a)/(h_c)
  
  
  
  
  
  require(BEDASSLE)
  
  
  
  
  
  
  
}
