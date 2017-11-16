calculate_basic_stats <- function(pop1,
                                  pop2,
                                  number_of_founders) {
  
  
  cat("starting estimating p table\n")
  pop_size <- length(pop1)
  cat("allele table pop 1\n")
  freq_pop1 <- matrix(nrow=length(pop1), ncol = number_of_founders * 2, 0)
  for(i in 1:length(pop1)) {
    vx <- calc_allele_frequencies(pop1[[i]], rep(0, number_of_founders * 2))
    freq_pop1[i,] <- vx;
  }
  freq_pop1 <- colMeans(freq_pop1)
  
  
  cat("allele table pop 2\n")
  freq_pop2 <- matrix(nrow=length(pop2), ncol = number_of_founders * 2, 0)
  
  for(i in 1:length(pop2)) {
    vx <- calc_allele_frequencies(pop2[[i]], rep(0, number_of_founders * 2))
    freq_pop2[i,] <- vx;
  }
  freq_pop2 <- colMeans(freq_pop2)
  final_table <- cbind(1:(2*number_of_founders), freq_pop1, freq_pop2)
  colnames(final_table) <- c("x","1","2")
  
  p <- list(X1 = final_table, X2 = final_table)
  n <- matrix(ncol = 2, nrow = 2, pop_size)
  
  cat("generated p table\n")
  
  #observed heterozygosity
  cat("calculating hst pop1\n")
  Ho_1 <- calculate_pop_heterozygosity(pop1)
  cat("calculating hst pop2\n")
  Ho_2 <- calculate_pop_heterozygosity(pop2)
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
  if (diploid) {
    rownames(sHo) <- loc.names
    rownames(Fis) <- loc.names
  }
  
  is.na(res) <- do.call(cbind, lapply(res, is.infinite))
  overall <- apply(res, 2, mean, na.rm = TRUE)
  overall[7] <- overall[4]/overall[3]
  overall[8] <- overall[6]/overall[5]
  overall[9] <- 1 - overall[1]/overall[2]
  overall[10] <- overall[6]/(1 - overall[2])
  names(overall) <- names(res)
  if (!diploid) {
    overall[-2] <- NA
  }
  all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round, 
                                                    digits), Ho = round(sHo, digits), Hs = round(Hs, digits), 
                  Fis = round(Fis, digits), perloc = round(res, digits), 
                  overall = round(overall, digits))
  class(all.res) <- "basic.stats"
  all.res
}