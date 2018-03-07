context("allele_frequencies")


test_that("calculate_allele_frequencies", {
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 5
  morgan <- 1

  found <- c();
  for (r in 1:100) {
    vx <- create_population(pop_size, 2,
                            run_time, morgan, r)

    testthat::expect_true(verify_population(vx))

    for (i in 1:pop_size) {
      found <- rbind(found, calc_allele_frequencies(vx[[i]],
                                                    alleles = rep(0, number_of_founders * 2))
      )
    }
  }

  v <- colMeans(found)
  expect_equal(v[[1]], 0.5, tolerance = 0.05)

  expect_equal(v[[2]], 0.5, tolerance = 0.05)

  found <- c();
  for (r in 1:100) {
    vx <- create_population(pop_size, 4,
                            run_time, morgan, r)

    testthat::expect_true(verify_population(vx))

    for (i in 1:pop_size) {
      found <- rbind(found,
                     calc_allele_frequencies(vx[[i]],
                                             alleles = rep(0, number_of_founders * 2))
      )
    }
  }

  v <- mean(colMeans(found))
  expect_equal(v, 0.25, tolerance = 0.05)
})