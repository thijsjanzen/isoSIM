context("create_populations")

test_that("create_population", {
  pop_size = 100
  number_of_founders = 10
  run_time = 100
  morgan = 1

  create_population(pop_size, number_of_founders, 
                    run_time, morgan, 42)
})

test_that("create_two_populations", {
  pop_size = 100
  number_of_founders = 10
  run_time = 100
  morgan = 1
  overlap = 0.5
  create_two_populations(pop_size, number_of_founders, 
                         run_time, morgan, 42, overlap)
})