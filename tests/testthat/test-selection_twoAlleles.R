context("selection two alleles")

test_that("create population", {
  s <- 0.1

  selection_matrix <- matrix(ncol=3, nrow = 1)
  selection_matrix[1,] <- c(0.5, 0, s)

  selected_pop_1 <- create_population_selection(pop_size = 1000,
                                                number_of_founders = 10,
                                                total_runtime = 100,
                                                morgan = 1,
                                                select_matrix =
                                                      selection_matrix,
                                                seed = 12345,
                                                track_frequency = TRUE)

  selection_matrix <- matrix(ncol=5, nrow = 1)
  selection_matrix[1,] <- c(0.5, 1.0, 1 + 0.5 * s, 1 + s, 0)

  selected_pop_2 <-
    create_population_selection_twoalleles(pop_size = 1000,
                                           number_of_founders = 10,
                                           total_runtime = 100,
                                           morgan = 1,
                                           select_matrix = selection_matrix,
                                           seed = 12345,
                                           track_frequency = TRUE)

  testthat::expect_true(all.equal(selected_pop_1, selected_pop_2))
})

test_that("select population two_alleles", {
  select_matrix <- matrix(ncol = 5, nrow = 1)
  s <- 0.1
  select_matrix[1, ] <- c(0.05, 1.0, 1+0.5*s, 1+s, 0)

  selected_pop <- isoSIM::create_population_selection_twoalleles(pop_size = 100,
                                                      number_of_founders = 10,
                                                      total_runtime = 100,
                                                      morgan = 1,
                                                      select_matrix,
                                                      seed = 1234)

  testthat::expect_equal(length(selected_pop$population), 100)
  testthat::expect_true(verify_population(selected_pop$population))
})

test_that("select on population", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                         number_of_founders = 10,
                                         total_runtime = 1000,
                                         morgan = 1,
                                         seed = 123)

  testthat::expect_true(verify_population(sourcepop))

  select_matrix <- matrix(ncol = 5, nrow = 1)
  s <- 0.1
  select_matrix[1, ] <- c(0.05, 1.0, 1+0.5*s, 1+s, 0)

  selected_pop <- select_population_twoalleles(sourcepop, select_matrix,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1233)

  testthat::expect_equal(length(selected_pop$population), 100)
  testthat::expect_true(verify_population(selected_pop$population))

  selected_pop <- select_population_twoalleles(sourcepop, select_matrix,
                                               pop_size = 100,
                                               total_runtime = 100,
                                               morgan = 1,
                                               seed = 1233,
                                               track_frequency = TRUE)

  testthat::expect_equal(length(selected_pop$population), 100)
  testthat::expect_true(verify_population(selected_pop$population))
})


test_that("select population two_alleles multiple markers", {
  select_matrix <- matrix(ncol = 5, nrow = 2)
  s <- 0.1
  select_matrix[1, ] <- c(0.25, 1.0, 1+0.5*s, 1+s, 0)
  select_matrix[1, ] <- c(0.75, 1.0, 1, 1+s,  1)

  selected_pop <- isoSIM::create_population_selection_twoalleles(pop_size = 100,
                                                                 number_of_founders = 10,
                                                                 total_runtime = 100,
                                                                 morgan = 1,
                                                                 select_matrix,
                                                                 seed = 1234)

  testthat::expect_equal(length(selected_pop$population), 100)
  testthat::expect_true(verify_population(selected_pop$population))
})



test_that("selection abuse", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                         number_of_founders = 2,
                                         total_runtime = 100,
                                         morgan = 1,
                                         seed = 123)

  testthat::expect_true(verify_population(sourcepop))

  select_matrix <- matrix(ncol = 5, nrow = 3)
  s <- 0.1
  select_matrix[1, ] <- c(0.05, 1.0, 1+0.5*s, 1+s, 0)
  select_matrix[2, ] <- c(0.15, 1.0, 1+0.5*s, 1+s, 0)

  testthat::expect_error(
    select_population_twoalleles(sourcepop, select_matrix,
                      pop_size = 1000,
                      total_runtime = 1000,
                      morgan = 1,
                      seed = 1234),
    "Can't start, there are NA values in the selection matrix!"
  )

  testthat::expect_error(
    create_population_selection_twoalleles(pop_size = 100,
                                number_of_founders = 10,
                                total_runtime = 10,
                                morgan = 1,
                                select_matrix,
                                seed = 1234),
    "Can't start, there are NA values in the selection matrix!"
  )

  select_matrix <- matrix(ncol = 3, nrow = 3)
  select_matrix[1, ] <- c(0.0, NA, 0)
  select_matrix[2, ] <- c(NA, 1.0, 1)


  testthat::expect_error(
    select_population_twoalleles(sourcepop, select_matrix,
                      pop_size = 1000,
                      total_runtime = 1000,
                      morgan = 1,
                      seed = 1234),
    "Can't start, there are NA values in the selection matrix!"
  )

  testthat::expect_error(
    create_population_selection_twoalleles(pop_size = 100,
                                number_of_founders = 10,
                                total_runtime = 10,
                                morgan = 1,
                                select_matrix,
                                seed = 1234),
    "Can't start, there are NA values in the selection matrix!"
  )


  select_matrix <- matrix(ncol = 5, nrow = 2)
  s <- 0.1
  select_matrix[1, ] <- c(0.05, 1.0, 1+0.5*s, 1+s, 0)
  select_matrix[2, ] <- c(0.15, 1.0, 1+0.5*s, 1+s, 0)

  testthat::expect_error(
    create_population_selection_twoalleles(pop_size = 100,
                                number_of_founders = 10,
                                total_runtime = 10,
                                morgan = 1,
                                select_matrix,
                                seed = 1234,
                                track_frequency = TRUE),
    "Can not track the frequency of more than one marker"
  )

  testthat::expect_error(
    select_population_twoalleles(sourcepop, select_matrix,
                      pop_size = 1000,
                      total_runtime = 1000,
                      morgan = 1,
                      seed = 1234,
                      track_frequency = TRUE),
    "Can not track the frequency of more than one marker"
  )
})