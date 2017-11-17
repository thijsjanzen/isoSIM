calculate_basic_stats <- function(pop1,
                                  pop2,
                                  number_of_founders) {
  
  digits <- 4
  
  pop_size <- length(pop1)
  
  # we calculate some statistics in C++ for speed 
  # (in R this takes ages for highly fragmented chromosomes)
  sum_stats_pop1 <- calculate_heterozygosity_and_freq_table(pop1, number_of_founders)
  sum_stats_pop2 <- calculate_heterozygosity_and_freq_table(pop2, number_of_founders)
  
  freq_pop1 <- sum_stats_pop1$freq_pop
  freq_pop2 <- sum_stats_pop2$freq_pop
  
  final_table <- cbind(1:(2*number_of_founders), freq_pop1, freq_pop2)
  colnames(final_table) <- c("x","1","2")
  
  p <- list(X1 = final_table, X2 = final_table)
  n <- matrix(ncol = 2, nrow = 2, pop_size)
  
  #observed heterozygosity
  Ho_1 <- sum_stats_pop1$Hst
  Ho_2 <- sum_stats_pop2$Hst
  sHo <- rbind( c(Ho_1,Ho_2), c(Ho_1,Ho_2))
  
  mHo <- apply(sHo, 1, mean, na.rm = TRUE)
  
  sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
  sp2 <- matrix(unlist(sp2), nrow = 2, byrow = TRUE)
  sp2 <- sp2[,-1]
  Hs <- (1 - sp2 - sHo/2/n)
  Hs <- n/(n - 1) * Hs
  Fis = 1 - sHo/Hs
  np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
  mn <- apply(n, 1, fun <- function(x) {
    np <- sum(!is.na(x))
    np/sum(1/x[!is.na(x)])
  })
  msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
  mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
  mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
  mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
  Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
  mFis = 1 - mHo/mHs
  Dst <- Ht - mHs
  Dstp <- np/(np - 1) * Dst
  Htp = mHs + Dstp
  Fst = Dst/Ht
  Fstp = Dstp/Htp
  Dest <- Dstp/(1 - mHs)
  res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst, 
                          Fstp, mFis, Dest))
  names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp", "Fst", 
                  "Fstp", "Fis", "Dest")

  is.na(res) <- do.call(cbind, lapply(res, is.infinite))
  overall <- apply(res, 2, mean, na.rm = TRUE)
  overall[7] <- overall[4]/overall[3]
  overall[8] <- overall[6]/overall[5]
  overall[9] <- 1 - overall[1]/overall[2]
  overall[10] <- overall[6]/(1 - overall[2])
  names(overall) <- names(res)
  all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round, 
                                                    digits), Ho = round(sHo, digits), Hs = round(Hs, digits), 
                  Fis = round(Fis, digits), perloc = round(res, digits), 
                  overall = round(overall, digits))
  all.res
}

calculate_heterozygosity_and_freq_table <- function(pop,
                                                 number_of_founders) {
  
  pop_for_cpp <- c();
  for(i in 1:length(pop)) {
    x <- pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    
    
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }
  
  results <- calculate_summaryStats(pop_for_cpp, number_of_founders);
  return(results);
}



calculate_heterozygosity <- function(pop) {

  #first we have to unwind all the individuals into one large vector
  pop_for_cpp <- c();
  for(i in 1:length(pop)) {
     x <- pop[[i]]$chromosome1
     chrom1 <- as.vector(t(x))
     x <- pop[[i]]$chromosome2
     chrom2 <- as.vector(t(x))
     
     
     pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }
  
  heterozygosity <- calc_heterozygosity_cpp(pop_for_cpp)
  return(heterozygosity)
}

hierf_basic_stats <- function(pop1, 
                              pop2, 
                              number_of_founders,
                              number_of_markers = 100,
                              markers_random = FALSE) {
  
  pop_size <- length(pop1)
  
  all_loci <- matrix(nrow= 2 * pop_size, ncol= 1 + number_of_markers, 0);
  all_loci[,1] <- c(rep(1, pop_size), rep(2, pop_size))
  colnames(all_loci) <- c("population", 1:number_of_markers)
  
  markers <- seq(1e-9, 1-(1e-9), number_of_markers)
  if(markers_random) {
    markers <- c();
    while(length(markers) < number_of_markers) {
      markers <- runif(number_of_markers, 0, 1)
      if(length(which(markers == 0.0))) {
          markers <- markers[- which(markers == 0.0)] #remove borders
      }
      if(length(which(markers == 1.0))) {
        markers <- markers[- which(markers == 1.0)]
      }
      #remove duplicates
      if(length(which(duplicated(markers)))) {
        markers <- markers[-which(duplicated(markers))]
      }
    }
    markers <- sort(markers)
  }
  
  for(x in 1:length(markers)) {
    focal_marker <- markers[x]
    for(i in 1:length(pop1)) {
      
      allele_1 <- 10 + findtype(pop1[[i]]$chromosome1, focal_marker)
      allele_2 <- 10 + findtype(pop1[[i]]$chromosome2, focal_marker)
      final_allele <- paste0(allele_1,allele_2)
      all_loci[i, x + 1] <- as.numeric(final_allele)
    }
    
    for(i in 1:length(pop2)) {
      
      allele_1 <- 10 + findtype(pop2[[i]]$chromosome1, focal_marker)
      allele_2 <- 10 + findtype(pop2[[i]]$chromosome2, focal_marker)
      
      final_allele <- paste0(allele_1,allele_2)
      all_loci[pop_size + i, x + 1] <- as.numeric(final_allele)
    }
  }
  
  require(hierfstat)
  
  hierf_sum <- hierfstat::basic.stats(all_loci)
  
  return(hierf_sum)
}