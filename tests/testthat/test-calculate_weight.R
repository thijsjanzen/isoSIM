context("sim_inf_chrom")


test_that("sim_inf_chrom use", {
  N <- 100
  Hzero <- 0.5
  maxtime <- 100
  size_in_Morgan <- 1
  markers <- -1
  seed <- 42
  sim_inf_chrom(popSize, Hzero, maxTime, size_in_Morgan, markers, seed)
})
