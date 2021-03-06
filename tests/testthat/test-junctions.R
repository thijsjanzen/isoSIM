
test_that("expected_number_junctions", {

  test_expected_junction_number <- function(pop_size,
                                            run_time,
                                            morgan,
                                            replicates)
  {
    cat(pop_size, run_time, morgan, "\n")
    found <- c()
    for (r in 1:replicates) {
      vx <- create_population(pop_size, 2,
                              run_time, morgan, r)

      testthat::expect_true(verify_population(vx))

      junct <- calculate_dist_junctions(vx)
      found <- c(found, mean(junct))
    }

    require(junctions)
    expected <- junctions::number_of_junctions(N = pop_size,
                                               H_0 = 0.5,
                                               C = morgan,
                                               t = run_time)

    testthat::expect_equal(mean(found), expected, tolerance = 1)
    cat(pop_size, run_time, morgan, mean(found), expected,"\n")
  }

  test_expected_junction_number(pop_size = 100, run_time = 100,
                                morgan = 1, replicates = 100)

  test_expected_junction_number(pop_size = 100, run_time = 100,
                                morgan = 0.5, replicates = 100)

  test_expected_junction_number(pop_size = 100, run_time = 100,
                                morgan = 3, replicates = 100)

  test_expected_junction_number(pop_size = 1000, run_time = 20,
                                morgan = 1, replicates = 100)

  vx <- create_population(pop_size = 1000, number_of_founders = 2,
                          total_runtime = 5, morgan = 1, seed = 666)

  plot_dist_junctions(vx)
})