context("create_populations")

test_that("create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  write_to_file <- TRUE

  vx <- create_population(pop_size, number_of_founders,
                    run_time, morgan, 42, write_to_file)
})

test_that("create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  write_to_file <- FALSE

  vx <- create_population(pop_size, number_of_founders,
                    run_time, morgan, 42, write_to_file)

  print(vx)
  print(vx[[1]])
  expect_equal(length(vx), pop_size)
  plot(vx[[1]])
})

test_that("calculate_heterozygosity", {
  pop_size <- 1000
  number_of_founders <- 2
  run_time <- 1
  morgan <- 1
  write_to_file <- FALSE

  vx <- create_population(pop_size, number_of_founders,
                               run_time, morgan, 42, write_to_file)

  avg_hetero <- 0
  for (i in 1:pop_size) {
    avg_hetero <- avg_hetero + calc_heterozygosity_indiv(vx[[i]])
  }
  avg_hetero <- avg_hetero / pop_size

  expect_equal(avg_hetero, 0.5, tolerance = 0.01)

  avg_hetero2 <- calculate_heterozygosity(vx)
  expect_equal(avg_hetero, avg_hetero2)

  pop_size <- 1000
  number_of_founders <- 2
  run_time <- 500
  morgan <- 1
  write_to_file <- FALSE

  vx <- create_population(pop_size, number_of_founders,
                               run_time, morgan, 42, write_to_file)

  avg_hetero <- 0
  for (i in 1:pop_size) {
    avg_hetero <- avg_hetero + calc_heterozygosity_indiv(vx[[i]])[[1]]
  }
  avg_hetero <- avg_hetero / pop_size

  expect_hetero <- 2 * 0.5 * 0.5 * (1 - 1 / (2 * pop_size)) ^ run_time

  expect_equal(avg_hetero, expect_hetero, tolerance = 0.05)

  avg_hetero2 <- calculate_heterozygosity(vx)
  expect_equal(avg_hetero, avg_hetero2, tolerance = 0.01)

  avg_hetero3 <- calculate_heterozygosity_and_freq_table(vx,
                                                         number_of_founders)

  testthat::expect_equal(avg_hetero2, avg_hetero3$Hst, tolerance = 0.01)
  testthat::expect_equal(avg_hetero3$freq_pop[1], 0.5, tolerance = 0.05)
  testthat::expect_equal(avg_hetero3$freq_pop[2], 0.5, tolerance = 0.05)

  testthat::expect_equal(avg_hetero3$freq_pop[1] + avg_hetero3$freq_pop[2],
                         1.0, tolerance = 0.05)
})

test_that("calculate_dist_junctions", {
  found <- c();
  pop_size <- 100
  run_time <- 100
  morgan <- 1
  for (r in 1:100) {
    vx <- create_population(pop_size, 2,
                                 run_time, morgan, r, FALSE)
    junct <- calculate_dist_junctions(vx)
    mean(junct)
    found <- c(found, mean(junct))
    cat(r, mean(junct), "\n")
  }
  K <- 2 * pop_size * 0.5
  expected <- K - K * (1 - 0.5 / K) ^ run_time

  expect_equal(mean(found), expected, tolerance = 1)

  vx <- create_population(pop_size, 2,
                               run_time, morgan, r, FALSE)
  plot_dist_junctions(vx)
})

test_that("calculate_allele_frequencies", {
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 5
  morgan <- 1
  write_to_file <- FALSE

  found <- c();
  for (r in 1:100) {
    vx <- create_population(pop_size, 2,
                               run_time, morgan, r, FALSE)
    for (i in 1:pop_size) {
      found <- rbind(found, calc_allele_frequencies(vx[[i]],
                            alleles = rep(0, number_of_founders * 2))
                    )
    }
  }

  v <- colMeans(found)
  expect_equal(v[[1]], 0.5, tolerance = 0.01)

  expect_equal(v[[2]], 0.5, tolerance = 0.01)

  found <- c();
  for (r in 1:100) {
    vx <- create_population(pop_size, 4,
                                 run_time, morgan, r, FALSE)
    for (i in 1:pop_size) {
      found <- rbind(found,
        calc_allele_frequencies(vx[[i]],
          alleles = rep(0, number_of_founders * 2))
      )
    }
  }

  v <- mean(colMeans(found))
  expect_equal(v, 0.25, tolerance = 0.01)
})

test_that("create_two_populations", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5
  write_to_file <- FALSE

  vx <- create_two_populations(pop_size, number_of_founders,
                         run_time, morgan, 42,
                         overlap, write_to_file)
})

test_that("create_two_full_populations", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5
  write_to_file <- FALSE

  vx <- create_two_populations(pop_size, number_of_founders,
                         run_time, morgan, 42,
                         overlap, write_to_file)

  expect_equal(length(vx$Population_1), pop_size)
  expect_equal(length(vx$Population_2), pop_size)

  print(vx$Population_1)
  print(vx$Population_2)
})

test_that("basic stats", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.0
  write_to_file <- FALSE

  vx <- create_two_populations(pop_size, number_of_founders,
                                    run_time, morgan, 42,
                                    overlap, write_to_file)

  pop1 <- vx$Population_1
  pop2 <- vx$Population_2

  a <- hierfstat_basic_stats(pop1, pop2,
                             number_of_markers = 100,
                             random_markers = TRUE)
})

test_that("stats", {
  pop_size <- 100
  number_of_founders <- 20
  run_time <- 1
  morgan <- 1
  overlap <- 0.1
  write_to_file <- FALSE

  vx <- create_two_populations(pop_size, number_of_founders,
                                    run_time, morgan, 42,
                                    overlap, write_to_file)

  pop1 <- vx$Population_1
  pop2 <- vx$Population_2

  number_of_markers <- 100
  sampled_individuals <- 10
  v1 <- isoSIM::calculate_fst(pop1, pop2,
                                 sampled_individuals,
                         number_of_markers, random_markers = TRUE)

  v2 <- isoSIM::calculate_fst(pop1, pop2,
                                 sampled_individuals,
                         number_of_markers, random_markers = FALSE)

  testthat::expect_equal(v1, v2, tolerance = 0.05)

  pop_size <- 100
  number_of_founders <- 10
  run_time <- 10000
  morgan <- 1
  overlap <- 0.0
  write_to_file <- FALSE

  vx <- create_two_populations(pop_size, number_of_founders,
                                    run_time, morgan, 42,
                                    overlap, write_to_file)

  pop1 <- vx$Population_1
  pop2 <- vx$Population_2

  number_of_markers <- 100
  v1 <- calculate_fst(pop1, pop2,
                         number_of_markers, random_markers = TRUE)

  testthat::expect_equal(1.0, v1, tolerance = 0.01)
})

test_that("continue_from_file", {
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 10
  morgan <- 1
  write_to_file <- FALSE

  vx1 <- create_population(pop_size, number_of_founders,
                    run_time, morgan, 42, write_to_file)

  total_runtime <- 100
  vx2 <- simulate_from_population("population_1.pop",
                           total_runtime,
                           morgan,
                           -1, 42)
})

test_that("migration",{
  pops_migration <- create_two_populations_migration(pop_size = 100,
                                                             number_of_founders = 10,
                                                             total_runtime = 1000,
                                                             morgan = 1,
                                                             seed = 1234,
                                                             migration = 0.0,
                                                             write_to_file = FALSE)

  testthat::expect_equal(length(pops_migration$Population_1), 100)
  testthat::expect_equal(length(pops_migration$Population_2), 100)
  testthat::expect_equal(length(pops_migration$Population_1),
                         length(pops_migration$Population_2))

  FST <- calculate_fst(pops_migration$Population_1,
                          pops_migration$Population_2,
                          sampled_individuals = 10,
                          number_of_markers = 100,
                          random_markers = TRUE)

  testthat::expect_equal(FST, 1.0, tolerance = 0.05)
})
