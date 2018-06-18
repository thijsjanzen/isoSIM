context("simulate")

test_that("simulate create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 10
  morgan <- 1

  vx <- create_population(pop_size, number_of_founders,
                          run_time, morgan, 42)

  testthat::expect_true(verify_population(vx))

  vy <- simulate_admixture(pop_size = pop_size, number_of_founders = number_of_founders,
                 total_runtime = run_time, morgan = morgan, seed = 42)

  testthat::expect_true(verify_population(vy$population))

  testthat::expect_true(all.equal(vx, vy$population))
})


test_that("simulate create_population_selection", {
  select_matrix <- matrix(ncol = 5, nrow = 1)
  s <- 0.1
  select_matrix[1, ] <- c(0.05, 1.0, 1+0.5*s, 1+s, 0)

  selected_pop <- create_population_selection(pop_size = 100,
                                              number_of_founders = 10,
                                              total_runtime = 100,
                                              morgan = 1,
                                              select_matrix,
                                              seed = 1234)

  testthat::expect_true(verify_population(selected_pop$population))

  selected_pop2 <- simulate_admixture(pop_size = 100,
                            number_of_founders = 10,
                            total_runtime = 100,
                            morgan = 1,
                            select_matrix = select_matrix,
                            seed = 1234)

  testthat::expect_true(verify_population(selected_pop2$population))

  testthat::expect_true(all.equal(selected_pop$population,
                                  selected_pop2$population))
})

test_that("simulate continue selection", {
  sourcepop <- create_population(pop_size = 100,
                                 number_of_founders = 10,
                                 total_runtime = 1000,
                                 morgan = 1,
                                 seed = 123)

  select_matrix <- matrix(ncol = 5, nrow = 1)
  s <- 0.1
  select_matrix[1, ] <- c(0.05, 1.0, 1+0.5*s, 1+s, 0)

  selected_pop <- select_population(sourcepop, select_matrix,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1233)

  selected_pop2 <- simulate_admixture(input_population = sourcepop,
                            select_matrix = select_matrix,
                            pop_size = 100,
                            total_runtime = 100,
                            morgan = 1,
                            seed = 1233)

  testthat::expect_true(verify_population(selected_pop2$population))

  testthat::expect_true(all.equal(selected_pop$population,
                                  selected_pop2$population))
})