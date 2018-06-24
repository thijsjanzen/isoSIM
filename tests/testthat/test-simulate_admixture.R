devtools::install_github("thijsjanzen/isoSIM")
library(isoSIM)

context("simulate_Admixture")

test_that("simulate admixture abuse", {




  vx <- isoSIM::simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = FALSE,
                           track_frequency = FALSE,
                           multiplicative_selection = TRUE)

  select_matrix <- matrix(NA, nrow=1,ncol=5)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            seed = 42,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            track_frequency = FALSE,
                                            multiplicative_selection = TRUE))

  select_matrix <- matrix(NA, nrow=1,ncol=3)
  testthat::expect_error(simulate_admixture(pop_size = 100,
                                            number_of_founders = 2,
                                            total_runtime = 100,
                                            morgan = 1,
                                            seed = 42,
                                            select_matrix = select_matrix,
                                            progress_bar = TRUE,
                                            track_junctions = FALSE,
                                            track_frequency = FALSE,
                                            multiplicative_selection = TRUE))

  track_freq <- c(0.4, 0.6, 100)
  vx <- isoSIM::simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = FALSE,
                           track_frequency = track_freq,
                           multiplicative_selection = TRUE)

  vx <- isoSIM::simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = TRUE,
                           track_frequency = FALSE,
                           multiplicative_selection = TRUE)


  vx <- isoSIM::simulate_admixture(pop_size = 100,
                           number_of_founders = 2,
                           total_runtime = 100,
                           morgan = 1,
                           seed = 42,
                           select_matrix = NA,
                           progress_bar = TRUE,
                           track_junctions = TRUE,
                           track_frequency = track_freq,
                           multiplicative_selection = TRUE)

})
