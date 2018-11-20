context("create_admixed_individuals")

testthat::test_that("create_admixed_individuals", {
  admixed_pop <- create_admixed_individuals(num_individuals = 10,
                                            population_size = 100,
                                            number_of_founders = 2,
                                            size_in_morgan = 1)

  testthat::expect_true( is(admixed_pop$population, "population") )

  expected_num_junctions <- junctions::number_of_junctions(N = 100, t = 1e6)
  found <- c()
  for(i in 1:length(admixed_pop$population)) {
    found <- c(found, length(admixed_pop$population[[1]]$chromosome1) - 2)
    found <- c(found, length(admixed_pop$population[[1]]$chromosome2) - 2)
  }
  avg_junctions <- mean(found)
  testthat::expect_equal(expected_num_junctions, avg_junctions, tolerance = 5)
})
