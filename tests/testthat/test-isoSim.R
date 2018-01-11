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

test_that("sim_inf_chrom accuracy", {
  pop_size <- 50
  init_heterozygosity <- 0.5
  max_time <- 1000
  morgan <- 1

  found <- c();
  for (r in 1:100) {
    vx <- sim_inf_chrom(pop_size = pop_size,
                initial_heterozygosity = init_heterozygosity,
                total_runtime = max_time,
                morgan = 1,
                markers = -1,
                seed = r)
    found <- rbind(found, as.numeric(vx$avgJunctions))
  }

  found <- colMeans(found)

  K <- 2 * init_heterozygosity * pop_size * morgan
  pred <- K - K * (1 - init_heterozygosity * morgan / K) ^ (0:(max_time - 1))

  rel_error <- abs(found[2:length(found)] / pred[2:length(pred)] - 1)

  for (i in 1:length(rel_error)) {
    expect_equal(rel_error[i], expected =  0, tolerance = 0.1)
  }

  pop_size <- 50
  init_heterozygosity <- 0.1
  max_time <- 500
  morgan <- 1

  found <- c();
  for (r in 1:100) {
    vx <- sim_inf_chrom(pop_size = pop_size,
                        initial_heterozygosity = init_heterozygosity,
                        total_runtime = max_time,
                        morgan = 1,
                        markers = -1,
                        seed = r)
    found <- rbind(found, as.numeric(vx$avgJunctions))
  }

  found <- colMeans(found)

  K <- 2 * init_heterozygosity * pop_size * morgan
  pred <- K - K * (1 - init_heterozygosity * morgan / K) ^ (0:(max_time - 1))

  rel_error <- abs(found[2:length(found)] / pred[2:length(pred)] - 1)

  for (i in 1:length(rel_error)) {
    expect_equal(rel_error[i], expected =  0, tolerance = 0.1)
  }
})