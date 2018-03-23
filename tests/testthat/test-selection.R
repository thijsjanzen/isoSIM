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

  select_matrix <- matrix(ncol = 3, nrow = 2)
  select_matrix[1, ] <- c(0.05, 0.1, 0)
  select_matrix[2, ] <- c(0.15, 0.5, 1)

  selected_pop <- select_population(sourcepop, select_matrix,
                                    selection = 1,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1233)

  testthat::expect_equal(length(selected_pop$population), 100)
  testthat::expect_true(verify_population(selected_pop$population))
})

test_that("selection abuse", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                              number_of_founders = 2,
                                              total_runtime = 100,
                                              morgan = 1,
                                              seed = 123)

  testthat::expect_true(verify_population(selected_pop$population))

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

  v <- calculate_marker_frequency(selected_pop$population, 0.5)

  a <- subset(v, v$ancestor == 0)
  testthat::expect_equal(a$frequency, 1)


  under_selection <- 0
  select_matrix <- matrix(ncol = 3, nrow = 1)
  select_matrix[1, ] <- c(0.5, under_selection, 0)

  selected_pop <- create_population_selection_markers(select_matrix,
                                                      pop_size = 1000,
                                                      number_of_founders = 20,
                                                      total_runtime = 100,
                                                      morgan = 1,
                                                      seed = 1234)

  v <- calculate_marker_frequency(selected_pop$population, 0.5)
  a <- which.max(v$frequency)
  testthat::expect_false(a == under_selection)
})


test_that("selection vector on existing population", {

sourcepop <- isoSIM::create_population(pop_size = 100,
                                       number_of_founders = 2,
                                       total_runtime = 100,
                                       morgan = 1,
                                       seed = 123)

testthat::expect_true(verify_population(sourcepop))

under_selection <- 0
select_matrix <- matrix(ncol = 3, nrow = 1)
select_matrix[1, ] <- c(0.5, under_selection, 1)

selected_pop <- isoSIM::select_population_markers(sourcepop,
                                                  select_matrix,
                                          pop_size = 1000,
                                          total_runtime = 1000,
                                          morgan = 1,
                                          seed = 12345)

testthat::expect_true(verify_population(selected_pop$population))

freq_output <- calculate_allele_frequencies(selected_pop$population,
                                            step_size = 0.001)

 #ggplot(freq_output, aes(x = location, y = frequency, col = as.factor(ancestor))) +
#    geom_line()


testthat::expect_equal(length(selected_pop), 1000)
testthat::expect_true(verify_population(selected_pop))

a <- subset(freq_output, location > 0.45 & location < 0.55)
b <- a %>%
  dplyr::group_by(as.factor(ancestor)) %>%
  dplyr::summarise("mean_freq" = mean(frequency))
v <- which.max(b$mean_freq)
testthat::expect_equal(v, under_selection + 1) #returns ancestor + 1
testthat::expect_equal(sum(b$mean_freq), 1, tolerance = 0.01)


under_selection <- 0
select_matrix <- matrix(ncol = 3, nrow = 2)
select_matrix[1, ] <- c(0.25, under_selection, 1)
select_matrix[2, ] <- c(0.75, under_selection + 1, 1)


selected_pop <- isoSIM::select_population_markers(sourcepop,
                                                  select_matrix,
                                                  pop_size = 1000,
                                                  total_runtime = 1000,
                                                  morgan = 1,
                                                  seed = 12345)

testthat::expect_true(verify_population(selected_pop$population))

freq_output <- calculate_allele_frequencies(selected_pop$population,
                                            step_size = 0.001)

testthat::expect_equal(length(selected_pop$population), 1000)
testthat::expect_true(verify_population(selected_pop$population))

a <- subset(freq_output, location > 0.20 & location < 0.30)
b <- a %>%
  dplyr::group_by(as.factor(ancestor)) %>%
  dplyr::summarise("mean_freq" = mean(frequency))
v <- which.max(b$mean_freq)
testthat::expect_equal(v, under_selection + 1) #returns ancestor + 1
testthat::expect_equal(sum(b$mean_freq), 1, tolerance = 0.01)

a <- subset(freq_output, location > 0.70 & location < 0.80)
b <- a %>%
  dplyr::group_by(as.factor(ancestor)) %>%
  dplyr::summarise("mean_freq" = mean(frequency))
v <- which.max(b$mean_freq)
testthat::expect_equal(v, under_selection + 2) #returns ancestor + 1
testthat::expect_equal(sum(b$mean_freq), 1, tolerance = 0.01)

})

test_that("track frequencies", {

  sourcepop <- isoSIM::create_population(pop_size = 100,
                                         number_of_founders = 2,
                                         total_runtime = 100,
                                         morgan = 1,
                                         seed = 123)

  testthat::expect_true(verify_population(sourcepop))

  under_selection <- 0
  select_matrix <- matrix(ncol = 3, nrow = 1)
  select_matrix[1, ] <- c(0.5, under_selection, 1)

  selected_pop <- isoSIM::select_population_markers(sourcepop,
                                                    select_matrix,
                                                    pop_size = 1000,
                                                    total_runtime = 1000,
                                                    morgan = 1,
                                                    seed = 12345,
                                                    track_frequency = TRUE)

  testthat::expect_true(verify_population(selected_pop$population))

  v <- subset(selected_pop$frequencies$time > 500 &
              selected_pop$frequencies$ancestor == under_selection)

  mean_freq <- mean(v$frequency)
  testthat::expect_equal(mean_freq, 1.0)
})





