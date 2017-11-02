context("create_populations")

test_that("create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1

  create_population(pop_size, number_of_founders, 
                    run_time, morgan, 42)
})

test_that("create_two_populations", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5
  create_two_populations(pop_size, number_of_founders, 
                         run_time, morgan, 42, overlap)
})

test_that("continue_from_file", {
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 10
  morgan <- 1
  
  create_population(pop_size, number_of_founders, 
                    run_time, morgan, 42)
  
  total_runtime = 100
  simulate_from_population("population_1.pop",
                           total_runtime,
                           morgan,
                           -1, 42)
})
  
  