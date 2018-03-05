create_population <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   write_to_file) {

  set.seed(seed)
  #call C_function
  pop <- isoSIM::create_population_cpp(pop_size,
                           number_of_founders,
                           total_runtime,
                           morgan,
                           write_to_file)

  popstruct <- isoSIM::create_pop_class(pop$population)
  return(popstruct)
}

simulate_from_population <- function(file_name,
                                     total_runtime,
                                     morgan,
                                     number_of_markers,
                                     seed) {
  set.seed(seed)
  simulate_from_population_cpp(file_name,
                               total_runtime,
                               morgan,
                               number_of_markers)
}





create_two_populations <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        overlap,
                                        write_to_file) {

  set.seed(seed)
  pops <- create_two_populations_cpp(pop_size,
                                 number_of_founders,
                                 total_runtime,
                                 morgan,
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


create_population_from_individuals <- function(individuals,
                                               pop_size,
                                               total_runtime,
                                               morgan,
                                               seed,
                                               write_to_file) {
  indiv <- c()

  for (j in seq_along(individuals)) {
    if (class(individuals[[j]]) != "individual") {
      stop("Input individuals not found\n are you sure you provided a list() ?")
    }
    for (i in seq_along(individuals[[j]]$chromosome1[, 1])) {
      indiv <- c(indiv, individuals[[j]]$chromosome1[i, ])
    }
    for (i in seq_along(individuals[[j]]$chromosome2[, 1])) {
      indiv <- c(indiv, individuals[[j]]$chromosome2[i, ])
    }
  }

  set.seed(seed)
  inbred_pop <- create_isofemale_line_cpp(indiv, pop_size,
                                  total_runtime, morgan)

  inbred_population <- isoSIM::create_pop_class(inbred_pop$population)

  return(inbred_population)
}