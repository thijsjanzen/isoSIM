context("selection")

test_that("select population", {
  select_matrix <- matrix(ncol = 3, nrow = 2)
  select_matrix[1, ] <- c(0.05, 0.1, 0)
  select_matrix[2, ] <- c(0.15, 0.5, 1)

  selected_pop <- isoSIM::create_population_selection(pop_size = 100,
                                                      number_of_founders = 10,
                                                      total_runtime = 100,
                                                      morgan = 1,
                                                      select_matrix,
                                                      selection = 1,
                                                      seed = 1234)

  testthat::expect_equal(length(selected_pop), 100)
  testthat::expect_true(verify_population(selected_pop))
})


test_that("select on population", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                    number_of_founders = 10,
                                    total_runtime = 1000,
                                    morgan = 1,
                                    seed = 123)

  testthat::expect_true(verify_population(sourcepop))

  select_matrix <- matrix(ncol = 3, nrow = 2)
  select_matrix[1, ] <- c(0.05, 0.1, 0)
  select_matrix[2, ] <- c(0.15, 0.5, 1)

  selected_pop <- select_population(sourcepop, select_matrix,
                                    selection = 1,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1233)

  testthat::expect_equal(length(selected_pop), 100)
  testthat::expect_true(verify_population(selected_pop))
})

test_that("allele frequencies", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                              number_of_founders = 2,
                                              total_runtime = 100,
                                              morgan = 1,
                                              seed = 123)

  testthat::expect_true(verify_population(sourcepop))

  select_matrix <- matrix(ncol = 3, nrow = 2)
  select_matrix[1, ] <- c(0.0, 0.5, 0)
  select_matrix[2, ] <- c(0.5, 1.0, 1)

  selected_pop <- isoSIM::select_population(sourcepop, select_matrix,
                                            selection = 1,
                                            pop_size = 1000,
                                            total_runtime = 1000,
                                            morgan = 1,
                                            seed = 12345)

  testthat::expect_true(verify_population(selected_pop))

  freq_output <- calculate_allele_frequencies(selected_pop,
                                                   step_size = 0.01)

  #ggplot(freq_output, aes(x = location, y = frequency, col = as.factor(ancestor))) + geom_line()

  testthat::expect_equal(length(unique(freq_output$ancestor)), 2)

  a <- subset(freq_output, freq_output$location < 0.5)
  require(magrittr)
  b <- a  %>%
        dplyr::group_by(as.factor(ancestor)) %>%
        dplyr::summarise("mean_freq" = mean(frequency))

  #testthat::expect_equal(b$mean_freq[[1]], 1.0, tolerance = 0.05)
  testthat::expect_gt(b$mean_freq[[1]], b$mean_freq[[2]])


  a <- subset(freq_output, freq_output$location > 0.5)
  b <- a %>%
    dplyr::group_by(as.factor(ancestor)) %>%
    dplyr::summarise("mean_freq" = mean(frequency))

  testthat::expect_gt(b$mean_freq[[2]], b$mean_freq[[1]])
  testthat::expect_equal(sum(b$mean_freq), 1)


  select_matrix <- matrix(ncol = 3, nrow = 1)

  under_selection <- 1

  select_matrix[1, ] <- c(0.2, 0.4, under_selection)

  selected_pop <- isoSIM::select_population(sourcepop, select_matrix,
                                            selection = 2,
                                            pop_size = 1000,
                                            total_runtime = 1000,
                                            morgan = 1,
                                            seed = 12345)

  testthat::expect_true(verify_population(selected_pop))

  freq_output <- calculate_allele_frequencies(selected_pop,
                                  step_size = 0.01)

  a <- subset(freq_output, location > 0.2 & location < 0.4)
  b <- a %>%
      dplyr::group_by(as.factor(ancestor)) %>%
      dplyr::summarise("mean_freq" = mean(frequency))
  v <- which.max(b$mean_freq)
  testthat::expect_equal(v, under_selection + 1) #returns ancestor + 1
  testthat::expect_equal(sum(b$mean_freq), 1, tolerance = 0.01)
})

test_that("selection abuse", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                              number_of_founders = 2,
                                              total_runtime = 100,
                                              morgan = 1,
                                              seed = 123)

  testthat::expect_true(verify_population(sourcepop))

  select_matrix <- matrix(ncol = 3, nrow = 3)
  select_matrix[1, ] <- c(0.0, 0.5, 0)
  select_matrix[2, ] <- c(0.5, 1.0, 1)


  testthat::expect_error(
        select_population(sourcepop, select_matrix,
                                            selection = 5,
                                            pop_size = 1000,
                                            total_runtime = 1000,
                                            morgan = 1,
                                            seed = 1234)

  )

  select_matrix <- matrix(ncol = 3, nrow = 3)
  select_matrix[1, ] <- c(0.0, NA, 0)
  select_matrix[2, ] <- c(NA, 1.0, 1)


  testthat::expect_error(
    select_population(sourcepop, select_matrix,
                      selection = 5,
                      pop_size = 1000,
                      total_runtime = 1000,
                      morgan = 1,
                      seed = 1234)

  )
})

test_that("selection vector", {
  under_selection <- 0
  select_matrix <- matrix(ncol = 3, nrow = 1)
  select_matrix[1, ] <- c(0.5, under_selection, 1)

  selected_pop <- create_population_selection_markers(select_matrix,
                                                      pop_size = 100,
                                                      number_of_founders = 10,
                                                      total_runtime = 100,
                                                      morgan = 1,
                                                      seed = 1234)

  freq_output <- calculate_allele_frequencies(selected_pop,
                                              step_size = 0.001)

 # ggplot(freq_output, aes(x = location, y = frequency, col = as.factor(ancestor))) +
#    geom_line()


  testthat::expect_equal(length(selected_pop), 100)
  testthat::expect_true(verify_population(selected_pop))

  a <- subset(freq_output, location > 0.45 & location < 0.55)
  b <- a %>%
    dplyr::group_by(as.factor(ancestor)) %>%
    dplyr::summarise("mean_freq" = mean(frequency))
  v <- which.max(b$mean_freq)
  testthat::expect_equal(v, under_selection + 1) #returns ancestor + 1
  testthat::expect_equal(sum(b$mean_freq), 1, tolerance = 0.01)
})

