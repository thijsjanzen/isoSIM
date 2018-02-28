create_population <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   write_to_file) {

  #call C_function
  pop <- isoSIM::create_population_cpp(pop_size,
                           number_of_founders,
                           total_runtime,
                           morgan,
                           seed,
                           write_to_file)

  popstruct <- isoSIM::create_pop_class(pop$population)
  return(popstruct)
}

create_two_populations <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        overlap,
                                        write_to_file) {

  pops <- create_two_populations_cpp(pop_size,
                                 number_of_founders,
                                 total_runtime,
                                 morgan,
                                 seed,
                                 overlap,
                                 write_to_file)

  pop1 <- isoSIM::create_pop_class(pops$population_1)
  pop2 <- isoSIM::create_pop_class(pops$population_2)

  output <- list("Population_1" = pop1,
                 "Population_2" = pop2)
  return(output)
}

create_two_populations_migration <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   migration,
                                   write_to_file) {

  pops <- create_two_populations_migration_cpp(pop_size,
                                     number_of_founders,
                                     total_runtime,
                                     morgan,
                                     seed,
                                     migration,
                                     write_to_file)

  pop1 <- isoSIM::create_pop_class(pops$population_1)
  pop2 <- isoSIM::create_pop_class(pops$population_2)

  output <- list("Population_1" = pop1,
                 "Population_2" = pop2)
  return(output)
}


create_population_from_individuals <- function(indiv_1,
                                               indiv_2,
                                               pop_size,
                                               total_runtime,
                                               morgan,
                                               seed,
                                               write_to_file) {
  indiv <- c()
  for (i in seq_along(indiv_1$chromosome1[, 1])) {
    indiv <- c(indiv, indiv_1$chromosome1[i, ])
  }
  for (i in seq_along(indiv_1$chromosome2[, 1])) {
    indiv <- c(indiv, indiv_1$chromosome2[i, ])
  }

  for (i in seq_along(indiv_2$chromosome1[, 1])) {
    indiv <- c(indiv, indiv_2$chromosome1[i, ])
  }
  for (i in seq_along(indiv_2$chromosome2[, 1])) {
    indiv <- c(indiv, indiv_2$chromosome2[i, ])
  }

  inbred_pop <- create_isofemale_line_cpp(indiv, pop_size,
                                  total_runtime, morgan,
                                  seed)

  inbred_population <- isoSIM::create_pop_class(inbred_pop$population)

  return(inbred_population)
}