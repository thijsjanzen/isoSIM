context("selection two alleles")

test_that("create population", {


  s <- 0.1

  selection_matrix = matrix(ncol=3, nrow = 1)
  selection_matrix[1,] = c(0.5, 0, s)

  selected_pop_1 <- create_population_selection(pop_size = 1000,
                                                number_of_founders = 10,
                                                total_runtime = 100,
                                                morgan = 1,
                                                select_matrix = selection_matrix,
                                                seed = 12345,
                                                track_frequency = TRUE)

  selection_matrix = matrix(ncol=5, nrow = 1)
  selection_matrix[1,] = c(0.5, 1.0, 1 + 0.5 * s, 1 + s, 0)

  selected_pop_2 <- create_population_selection_twoalleles(pop_size = 1000,
                                                           number_of_founders = 10,
                                                           total_runtime = 100,
                                                           morgan = 1,
                                                           select_matrix = selection_matrix,
                                                           seed = 12345,
                                                           track_frequency = TRUE)

  testthat::expect_true(all.equal(selected_pop_1, selected_pop_2));
})