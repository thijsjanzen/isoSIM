context("create_populations")

test_that("create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 10
  morgan <- 1

  vx <- create_population(pop_size, number_of_founders,
                    run_time, morgan, 42)

  testthat::expect_true(verify_population(vx))
})

test_that("create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1

  vx <- create_population(pop_size, number_of_founders,
                    run_time, morgan, 42)

  testthat::expect_true(verify_population(vx))

  print(vx)
  print(vx[[1]])
  expect_equal(length(vx), pop_size)
  plot(vx[[1]])
})


test_that("migration",{
  pops_migration <- create_two_populations_migration(pop_size = 100,
                                                     number_of_founders = 10,
                                                     total_runtime = 1000,
                                                     morgan = 1,
                                                     seed = 1234,
                                                     migration = 0.0)

  testthat::expect_true(verify_population(pops_migration$Population_1))
  testthat::expect_true(verify_population(pops_migration$Population_2))


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

  pops_migration <- create_two_populations_migration(pop_size = 100,
                                                     number_of_founders = 10,
                                                     total_runtime = 1000,
                                                     morgan = 1,
                                                     seed = 1234,
                                                     migration = 0.1)

  testthat::expect_true(verify_population(pops_migration$Population_1))
  testthat::expect_true(verify_population(pops_migration$Population_2))
})