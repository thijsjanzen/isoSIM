context("sim_inf_chrom")

test_that("sim_inf_chrom use", {
  pop_size <- 100
  init_heterozygosity <- 0.5
  max_time <- 100
  morgan <- 1

  sim_inf_chrom(pop_size = pop_size,
                initial_heterozygosity = init_heterozygosity,
                total_runtime = max_time,
                morgan = morgan,
                markers = -1,
                seed = 42)

  #test with a number of markers
  sim_inf_chrom(pop_size = pop_size,
                initial_heterozygosity = init_heterozygosity,
                total_runtime = max_time,
                morgan = morgan,
                markers = 1000,
                seed = 42)

  #test with extremely high recombination rate
  sim_inf_chrom(pop_size = pop_size,
                initial_heterozygosity = init_heterozygosity,
                total_runtime = max_time,
                morgan = 20,
                markers = -1,
                seed = 42)

  #test with extremely low recombination rate
  sim_inf_chrom(pop_size = pop_size,
                initial_heterozygosity = init_heterozygosity,
                total_runtime = max_time,
                morgan = 1e-7,
                markers = -1,
                seed = 42)
})