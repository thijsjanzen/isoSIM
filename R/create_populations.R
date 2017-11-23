create_full_population <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   write_to_file) {

  #call C_function
  pop <- isoSIM::create_population(pop_size,
                           number_of_founders,
                           total_runtime,
                           morgan,
                           seed,
                           write_to_file)

  popstruct <- isoSIM::create_pop_class(pop$population)
  return(popstruct)
}

create_two_full_populations <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        overlap,
                                        write_to_file) {

  pops <- create_two_populations(pop_size,
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