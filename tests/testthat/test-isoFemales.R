context("isoFemale creation")

test_that("create_isofemale", {

  pop_size <- 100
  number_of_founders <- 2
  run_time <- 100
  morgan <- 1
  write_to_file <- FALSE

  pop <- create_population(pop_size, number_of_founders,
                           run_time, morgan, 42, write_to_file)

  females <- create_iso_female(pop, n = 1)

  testthat::expect_equal(length(females), 1)
})

test_that("create_population_from_isofemales", {

  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5
  write_to_file <- FALSE

  vx <- isoSIM::create_two_populations(pop_size, number_of_founders,
                               run_time, morgan, 42,
                               overlap, write_to_file)

  female_1 <- isoSIM::create_iso_female(vx$Population_1, n = 1)
  female_2 <- isoSIM::create_iso_female(vx$Population_2, n = 1)

  vy <- isoSIM::create_population_from_individuals(female_1[[1]], female_2[[1]],
                                     pop_size, run_time,
                                     morgan,
                                     seed = 666,
                                     write_to_file)
  testthat::expect_equal(length(vy), pop_size)

  isoSIM::plot_chromosome(female_1[[1]]$chromosome1, 0, 1)
})

test_that("cpp classes", {
  test_fish_functions()
})