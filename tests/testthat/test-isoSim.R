context("sim_inf_chrom")


test_that("sim_inf_chrom use", {
  N <- 100
  H0 <- 0.5
  maxT <- 100
  M <- 1

  sim_inf_chrom(popSize = N, 
                Hzero = H0, 
                maxTime = maxT, 
                size_in_Morgan = M, 
                markers = -1, 
                seed = 42)
  
  
  #test with a number of markers
  sim_inf_chrom(popSize = N, 
                Hzero = H0, 
                maxTime = maxT, 
                size_in_Morgan = M, 
                markers = 1000, 
                seed = 42)
  
  #test with extremely high recombination rate
  sim_inf_chrom(popSize = N, 
                Hzero = H0, 
                maxTime = maxT, 
                size_in_Morgan = 20, 
                markers = -1, 
                seed = 42)
  
  #test with extremely low recombination rate
  sim_inf_chrom(popSize = N, 
                Hzero = H0, 
                maxTime = maxT, 
                size_in_Morgan = 1e-7, 
                markers = -1, 
                seed = 42)
})
