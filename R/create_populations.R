create_population <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   progress_bar = TRUE,
                                   track_junctions = FALSE) {

  #call C_function
  set.seed(seed)
  pop <- create_population_cpp(pop_size,
                           number_of_founders,
                           total_runtime,
                           morgan, progress_bar,
                           track_junctions)

  popstruct <- isoSIM::create_pop_class(pop$population)
  if(track_junctions == FALSE) return(popstruct)

  output <- list("population" = popstruct,
                 "junctions" = pop$junctions)
  return(output)
}

create_two_populations_migration <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed,
                                   migration,
                                   progress_bar = TRUE) {

  set.seed(seed)
  pops <- create_two_populations_migration_cpp(pop_size,
                                     number_of_founders,
                                     total_runtime,
                                     morgan,
                                     migration,
                                     progress_bar)

  pop1 <- isoSIM::create_pop_class(pops$population_1)
  pop2 <- isoSIM::create_pop_class(pops$population_2)

  output <- list("Population_1" = pop1,
                 "Population_2" = pop2)
  return(output)
}


create_population_from_individuals <- function(individuals,
                                               pop_size = 100,
                                               total_runtime = 2000,
                                               morgan,
                                               seed,
                                               progress_bar = TRUE) {

  pop_for_cpp <- population_to_vector(individuals)

  set.seed(seed)
  inbred_pop <- create_isofemale_line_cpp(pop_for_cpp, pop_size,
                                  total_runtime, morgan, progress_bar)

  inbred_population <- isoSIM::create_pop_class(inbred_pop$population)

  return(inbred_population)
}


save_population <- function(population, file_name, compression = TRUE) {

  if(class(population) != "population") {
    stop("Not providing a population structure")
  }

  saveRDS(population, file = file_name, compress = compression)
}

load_population <- function(file_name) {
  readRDS(file_name)
}