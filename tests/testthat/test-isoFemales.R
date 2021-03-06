context("isoFemale creation")

test_that("create_isofemale", {

  pop_size <- 100
  number_of_founders <- 2
  run_time <- 100
  morgan <- 1
  write_to_file <- FALSE

  pop <- create_population(pop_size, number_of_founders,
                           run_time, morgan, 42)

  testthat::expect_true(verify_population(pop))

  females <- create_iso_female(pop, n = 5, run_time = 3000)

  females <- create_iso_female(pop, n = 1, run_time = 3000)


  testthat::expect_equal(length(females), 1)
})

test_that("create_population_from_isofemales", {

  pop_size <- 100
  number_of_founders <- 10
  run_time <- 100
  morgan <- 1
  overlap <- 0.5

  pop1 <- create_population(pop_size, number_of_founders,
                            run_time, morgan, 42)

  pop2 <- create_population(pop_size, number_of_founders,
                            run_time, morgan, 24)
  pop2 <- increase_ancestor(pop2, number_of_founders)

  testthat::expect_true(verify_population(pop1))
  testthat::expect_true(verify_population(pop2))

  female_1 <- create_iso_female(pop1, n = 1,
                                run_time = 2000, seed = 1)
  female_2 <- create_iso_female(pop2, n = 1,
                                run_time = 2000, seed = 2)



  testthat::expect_true(verify_individual(female_1[[1]]))
  testthat::expect_true(verify_individual(female_2[[1]]))


  females <- create_iso_female(pop1, n = 2, run_time = 2000)

  vy <- create_population_from_individuals(females,
                                           pop_size, 2000,
                                           morgan,
                                           seed = 666)

  testthat::expect_equal(length(vy), pop_size)
  testthat::expect_true(verify_population(vy))

  vy <- create_population_from_individuals(list(female_1[[1]],
                                                female_2[[1]]),
                                                pop_size,
                                                2000,
                                                morgan,
                                                seed = 666)

  testthat::expect_equal(length(vy), pop_size)
  testthat::expect_true(verify_population(vy))

  plot_chromosome(female_1[[1]]$chromosome1, 0, 1)
})

test_that("cpp classes", {
  test_fish_functions()

  a <- matrix(c(0.1,1, 2, 2), nrow = 2)
  b <- matrix(c(0,1,1,-1), nrow = 2)
  indiv <- list(chromosome1 = a, chromosome2 = a)
  class(indiv) <- "individual"

  # chromosome 1
  testthat::expect_output(v <- verify_individual(indiv),
                           "Chromosome doesn't start at 0")

  indiv <- list(chromosome1 = b, chromosome2 = a)
  class(indiv) <- "individual"

  # chromosome 2
  testthat::expect_output(v <- verify_individual(indiv),
                          "Chromosome doesn't start at 0")

  a <- matrix(c(0.0, 1, 2, 2), nrow = 2)
  b <- matrix(c(0,1,1,-1), nrow = 2)
  indiv <- list(chromosome1 = b, chromosome2 = a)
  class(indiv) <- "individual"
  testthat::expect_output(v <- verify_individual(indiv),
                          "Chromosome doesn't end with -1")
  indiv$chromosome2 <-  indiv$chromosome1
  indiv$chromosome1 <- a
  testthat::expect_output(v <- verify_individual(indiv),
                          "Chromosome doesn't end with -1")


  a <- matrix(c(0.0, 1, 0.5, 29192875037,  1, -1), ncol = 2)

  indiv$chromosome1 <- a
  testthat::expect_output(v <- verify_individual(indiv),
                          "Memory error recorded in chromosome")

  a <- matrix(c(0.0, 1, 0.5, -92875037,  1, -1), ncol = 2)
  indiv$chromosome2 <- a
  indiv$chromosome1 <- b
  testthat::expect_output(v <- verify_individual(indiv),
                          "Memory error recorded in chromosome")

})


test_that("create_population_from_individuals", {

  pop_size <- 100
  number_of_founders <- 20
  run_time <- 100
  morgan <- 1


  pop1 <- create_population(pop_size, number_of_founders,
                            run_time, morgan, 42)

  pop2 <- create_population(pop_size, number_of_founders,
                            run_time, morgan, 24)
  pop2 <- increase_ancestor(pop2, number_of_founders)

  testthat::expect_true(verify_population(pop1))
  testthat::expect_true(verify_population(pop2))

  isofemale_1 <- create_iso_female(source_pop = pop1,
                                   n = 1,
                                   inbreeding_pop_size = 100,
                                   run_time = 2000,
                                   morgan = 1)

  isofemale_2 <- create_iso_female(source_pop = pop2,
                                   n = 1,
                                   inbreeding_pop_size = 100,
                                   run_time = 2000,
                                   morgan = 1)

  testthat::expect_true(verify_individual(isofemale_1[[1]]))
  testthat::expect_true(verify_individual(isofemale_2[[1]]))

  input <- list(isofemale_1[[1]],
                isofemale_2[[1]])

  mixed_population <- create_population_from_individuals(list(isofemale_1[[1]],
                                                              isofemale_2[[1]]),
                                                         pop_size = 100,
                                                         total_runtime = 100,
                                                         morgan = 1,
                                                         seed = 44)

  testthat::expect_true(verify_population(mixed_population))

  mixed_population <- create_population_from_individuals(input,
                                                         pop_size = 100,
                                                         total_runtime = 100,
                                                         morgan = 1,
                                                         seed = 43)

  testthat::expect_true(verify_population(mixed_population))

  isofemales <- create_iso_female(source_pop = pop1,
                                  n = 5,
                                  inbreeding_pop_size = 100,
                                  run_time = 2000,
                                  morgan = 1)

  testthat::expect_true(verify_individual(isofemales[[1]]))
  testthat::expect_true(verify_individual(isofemales[[2]]))
  testthat::expect_true(verify_individual(isofemales[[3]]))
  testthat::expect_true(verify_individual(isofemales[[4]]))
  testthat::expect_true(verify_individual(isofemales[[5]]))

  mixed_population_2 <- create_population_from_individuals(isofemales,
                                                           pop_size = 100,
                                                           total_runtime = 2000,
                                                           morgan = 1,
                                                           seed = 42)

  testthat::expect_true(verify_population(mixed_population_2))
})