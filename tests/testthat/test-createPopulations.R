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

test_that("expected_number_junctions", {
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
  require(junctions)
  expected <- junctions::number_of_junctions(N = pop_size,
                                             H_0 = 0.5,
                                             C = 1,
                                             t = run_time)

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


test_that("fst", {
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

test_that("create_population_from_individuals", {
  two_populations <- create_two_populations(pop_size = 100,
                                            number_of_founders = 20,
                                            total_runtime = 5,
                                            morgan = 1,
                                            seed = 42,
                                            overlap = 0.25,
                                            write_to_file = FALSE)

  isofemale_1 <- create_iso_female(source_pop = two_populations$Population_1,
                                   n = 1,
                                   inbreeding_pop_size = 100,
                                   run_time = 1000000,
                                   morgan = 1)

  isofemale_2 <- create_iso_female(source_pop = two_populations$Population_2,
                                   n = 1,
                                   inbreeding_pop_size = 100,
                                   run_time = 1000000,
                                   morgan = 1)

  mixed_population <- create_population_from_individuals(list(isofemale_1[[1]],
                                                              isofemale_2[[1]]),
                                                         pop_size = 100,
                                                         total_runtime = 100,
                                                         morgan = 1,
                                                         seed = 42,
                                                         write_to_file = FALSE)


  a1 <- list(isofemale_1[[1]],
             isofemale_2[[1]])

  a2 <- c(isofemale_1[[1]],
               isofemale_2[[1]])


  isofemales <- create_iso_female(source_pop = two_populations$Population_1,
                                   n = 5,
                                   inbreeding_pop_size = 100,
                                   run_time = 1000000,
                                   morgan = 1)

  mixed_population_2 <- create_population_from_individuals(isofemales,
                                                         pop_size = 100,
                                                         total_runtime = 100,
                                                         morgan = 1,
                                                         seed = 42,
                                                         write_to_file = FALSE)

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

  fst <- calculate_fst(pops_migration$Population_1,
                          pops_migration$Population_2,
                          sampled_individuals = 10,
                          number_of_markers = 100,
                          random_markers = TRUE)

  testthat::expect_equal(fst, 1.0, tolerance = 0.05)
})
