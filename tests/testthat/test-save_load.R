context("create_populations")


test_that("save_population", {
  pop_size <- 100
  number_of_founders <- 10
  run_time <- 10
  morgan <- 1

  vx <- create_population(pop_size, number_of_founders,
                          run_time, morgan, 42)

  testthat::expect_true(verify_population(vx))

  save_population(vx, file_name = "test.pop")

  vy <- load_population(file_name = "test.pop")

  testthat::expect_true(verify_population(vy))

  testthat::expect_equal(length(vx), length(vy))

  for(i in 1:length(vx)) {
    testthat::expect_true(all.equal(vx[[i]], vy[[i]]))
  }
})