context("create_populations")

test_that("create_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  write_to_file <- TRUE

  create_population(pop_size, number_of_founders, 
                    run_time, morgan, 42, write_to_file)
})

test_that("create_full_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  write_to_file <- FALSE
  
  vx <- create_full_population(pop_size, number_of_founders, 
                    run_time, morgan, 42, write_to_file)
  
  print(vx$Population)
  print(vx$Population[[1]])
  expect_equal(length(vx$Population), pop_size)
  plot(vx$Population[[1]])
})

test_that("calculate_heterozygosity", {
  pop_size <- 1000
  number_of_founders <- 2
  run_time <- 1
  morgan <- 1
  write_to_file <- FALSE
  
  vx <- create_full_population(pop_size, number_of_founders, 
                               run_time, morgan, 42, write_to_file)
  
  avg_hetero <- 0
  for(i in 1:pop_size) {
    avg_hetero <- avg_hetero + calc_heterozygosity(vx$Population[[i]])
  }
  avg_hetero <- avg_hetero / pop_size
  
  expect_equal(avg_hetero, 0.5, tolerance=0.01)
  
  avg_hetero2 <- calculate_pop_heterozygosity(vx)
  expect_equal(avg_hetero, avg_hetero2)
})




test_that("create_two_populations", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5
  write_to_file <- FALSE
  
  create_two_populations(pop_size, number_of_founders, 
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
  
  vx <- create_two_full_populations(pop_size, number_of_founders, 
                         run_time, morgan, 42, 
                         overlap, write_to_file)
  
  expect_equal(length(vx$Population_1), pop_size)
  expect_equal(length(vx$Population_2), pop_size)
  
  print(vx$Population_1)
  print(vx$Population_2)
})


test_that("continue_from_file", {
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 10
  morgan <- 1
  write_to_file <- FALSE
  
  create_population(pop_size, number_of_founders, 
                    run_time, morgan, 42, write_to_file)
  
  total_runtime = 100
  simulate_from_population("population_1.pop",
                           total_runtime,
                           morgan,
                           -1, 42)
})
  
  