context("isoFemale creation")

test_that("create_isofemale", {
  
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 100
  morgan <- 1
  write_to_file <- FALSE
  
  pop <- create_full_population(pop_size, number_of_founders, 
                           run_time, morgan, 42, write_to_file)
  
  females <- create_isoFemale(pop, n = 1)
  
  expect_equal(length(females), 1)
})