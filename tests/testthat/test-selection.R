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
                                                      seed = 1234,
                                                      write_to_file = FALSE)

  testthat::expect_equal(length(selected_pop), 100)
})


test_that("select on population", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                    number_of_founders = 10,
                                    total_runtime = 1000,
                                    morgan = 1,
                                    seed = 123,
                                    write_to_file = FALSE)

  select_matrix <- matrix(ncol = 3, nrow = 2)
  select_matrix[1, ] <- c(0.05, 0.1, 0)
  select_matrix[2, ] <- c(0.15, 0.5, 1)

  selected_pop <- isoSIM::select_population(sourcepop, select_matrix,
                                    selection = 0.1,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1234,
                                    write_to_file = FALSE)

  testthat::expect_equal(length(selected_pop), 100)
})

test_that("allele frequencies", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                              number_of_founders = 2,
                                              total_runtime = 100,
                                              morgan = 1,
                                              seed = 123)

  select_matrix <- matrix(ncol = 3, nrow = 2)
  select_matrix[1, ] <- c(0.0, 0.5, 0)
  select_matrix[2, ] <- c(0.5, 1.0, 1)

  selected_pop <- isoSIM::select_population(sourcepop, select_matrix,
                                            selection = 5,
                                            pop_size = 1000,
                                            total_runtime = 1000,
                                            morgan = 1,
                                            seed = 1234)

  freq_output <- calculate_allele_frequencies(selected_pop,
                                                   number_of_founders = 2,
                                                   step_size = 0.01)

  a <- subset(freq_output, freq_output$location < 0.5)
  require(magrittr)
  b <- a  %>%
        dplyr::group_by(as.factor(ancestor)) %>%
        dplyr::summarise("mean_freq" = mean(frequency))

  testthat::expect_equal(b$mean_freq[[1]], 1.0, tolerance = 0.05)

  a <- subset(freq_output, freq_output$location > 0.5)
  b <- a %>%
    dplyr::group_by(as.factor(ancestor)) %>%
    dplyr::summarise("mean_freq" = mean(frequency))
  testthat::expect_equal(b$mean_freq[[2]], 1.0, tolerance = 0.05)


  number_founders <- 20
  sourcepop <- isoSIM::create_population(
                          pop_size = 10000,
                          number_of_founders = number_founders,
                          total_runtime = 1,
                          morgan = 1,
                          seed = 123)

  freq_output <- calculate_allele_frequencies(sourcepop,
                                  number_of_founders = number_founders,
                                  step_size = 0.01)

  b <- freq_output %>%
    dplyr::group_by(as.factor(ancestor)) %>%
    dplyr::summarise("mean_freq" = mean(frequency))

  testthat::expect_equal(mean(b$mean_freq), 1 / number_founders, tolerance = 0.01)

  number_founders <- 5
  sourcepop <- isoSIM::create_population(pop_size = 1000,
                                  number_of_founders = number_founders,
                                  total_runtime = 1000,
                                  morgan = 1,
                                  seed = 123)

  freq_output <- calculate_allele_frequencies(sourcepop,
                                  number_of_founders = number_founders,
                                  step_size = 0.01)
  b <- freq_output %>%
    dplyr::group_by(as.factor(ancestor)) %>%
    dplyr::summarise("mean_freq" = mean(frequency))

  testthat::expect_equal(mean(b$mean_freq), 1 / number_founders, tolerance = 0.01)

  number_founders <- 20
  sourcepop <- isoSIM::create_population(pop_size = 1000,
                                number_of_founders = number_founders,
                                total_runtime = 1,
                                morgan = 1,
                                seed = 123)

  select_matrix <- matrix(ncol = 3, nrow = 1)

  under_selection <- 1

  select_matrix[1, ] <- c(0.2, 0.4, under_selection)

  selected_pop <- isoSIM::select_population(sourcepop, select_matrix,
                                            selection = 2,
                                            pop_size = 100,
                                            total_runtime = 1000,
                                            morgan = 1,
                                            seed = 12345)

  freq_output <- calculate_allele_frequencies(selected_pop,
                                  number_of_founders = number_founders,
                                  step_size = 0.001)

  a <- subset(freq_output, location > 0.2 & location < 0.4)
  b <- a %>%
      dplyr::group_by(as.factor(ancestor)) %>%
      dplyr::summarise("mean_freq" = mean(frequency))
  v <- which.max(b$mean_freq)
  testthat::expect_equal(v, under_selection + 1) #returns ancestor + 1
})

test_that("selection abuse", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                              number_of_founders = 2,
                                              total_runtime = 100,
                                              morgan = 1,
                                              seed = 123)

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
})
