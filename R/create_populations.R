
create_full_population <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   writeToFile) {

  #call C_function
  pop <- create_population(pop_size,
                           number_of_founders,
                           total_runtime,
                           morgan,
                           seed,
                           writeToFile)

  popStruct <- create_pop_class(pop$population);
  return(list("Population" = popStruct));
}

create_two_full_populations <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        overlap,
                                        writeToFile) {

  pops <- create_two_populations(pop_size,
                                 number_of_founders,
                                 total_runtime, 
                                 morgan,
                                 seed,
                                 overlap,
                                 writeToFile)

  pop1 <- create_pop_class(pops$population_1)
  pop2 <- create_pop_class(pops$population_2)

  output <- list("Population_1" = pop1,
                 "Population_2" = pop2)
  return(output)
}

