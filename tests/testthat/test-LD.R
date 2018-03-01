context("LD stats")

test_that("calculate_average_LD", {
  pop_size <- 100
  number_of_founders <- 2
  run_time <- 1
  morgan <- 1
  write_to_file <- FALSE

  pop1 <- create_population(pop_size, number_of_founders,
                                 run_time, morgan, 42, write_to_file)

  number_of_markers <- 2
  markers <- c(0.25, 0.26)

  all_loci <- matrix(nrow = length(pop1), ncol = 2 * number_of_markers, 0);

  for (x in 1:length(markers)) {
    focal_marker <- markers[x]
    for (i in 1:length(pop1)) {
      allele_1 <- 1 + isoSIM::findtype(pop1[[i]]$chromosome1, focal_marker)
      allele_2 <- 1 + isoSIM::findtype(pop1[[i]]$chromosome2, focal_marker)

      index <- (x - 1) * 2 + 1

      all_loci[i, index]     <- as.numeric(allele_1)
      all_loci[i, index + 1] <- as.numeric(allele_2)
    }
  }

  g1 <- all_loci[, 1:2]
  g2 <- all_loci[, 3:4]

  vv <- calculate_average_LD(g1, g2)
  expect_equal(vv$LD, 1)
})

test_that("calculate_LD_matrix", {
  pop_size <- 100
  number_of_founders <- 2
  sampled_individuals <- pop_size
  run_time <- 1
  morgan <- 1
  write_to_file <- FALSE

  pop1 <- create_population(pop_size, number_of_founders,
                                 run_time, morgan, 42, write_to_file)

  vv <- calculate_LD(pop1, sampled_individuals,
                     number_of_markers = 10, random_markers = TRUE)

  vv1 <- as.vector(vv$LD_matrix[!is.na(vv$LD_matrix)])
  vv2 <- as.vector(vv$dist_matrix[!is.na(vv$dist_matrix)])

  linear_model <- lm(vv1 ~ vv2)
  testthat::expect_equal(linear_model$coefficients[[1]], 1, tolerance = 0.05)

  #it should at least be negative
  testthat::expect_equal(linear_model$coefficients[[2]], -0.5, tolerance = 0.49)
})