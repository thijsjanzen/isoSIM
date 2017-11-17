calculate_basic_stats <- function(pop1,
                                  pop2,
                                  number_of_founders) {
  
  digits <- 4
  cat("starting estimating p table\n")
  
  
  pop_size <- length(pop1)
  
  sum_stats_pop1 <- calculate_heterozygosity_and_freq_table(pop1, number_of_founders)
  sum_stats_pop2 <- calculate_heterozygosity_and_freq_table(pop2, number_of_founders)
  
  freq_pop1 <- sum_stats_pop1$freq_pop
  freq_pop2 <- sum_stats_pop2$freq_pop
  
  cat("allele table pop 1\n")
  
  final_table <- cbind(1:(2*number_of_founders), freq_pop1, freq_pop2)
  colnames(final_table) <- c("x","1","2")
  
  p <- list(X1 = final_table, X2 = final_table)
  n <- matrix(ncol = 2, nrow = 2, pop_size)
  
  cat("generated p table\n")
  
  #observed heterozygosity
  cat("calculating hst pop1\n")
  Ho_1 <- sum_stats_pop1$Hst
  cat("calculating hst pop2\n")
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






