create_population <- function(pop_size,
                                   number_of_founders,
                                   total_runtime,
                                   morgan,
                                   seed) {

  set.seed(seed)
  #call C_function
  pop <- isoSIM::create_population_cpp(pop_size,
                           number_of_founders,
                           total_runtime,
                           morgan)

  popstruct <- isoSIM::create_pop_class(pop$population)
  return(popstruct)
}

create_two_populations <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        overlap) {

  set.seed(seed)
  pops <- create_two_populations_cpp(pop_size,
                                 number_of_founders,
                                 total_runtime,
                                 morgan,
                                 overlap)

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
                                   migration) {

  set.seed(seed)
  pops <- create_two_populations_migration_cpp(pop_size,
                                     number_of_founders,
                                     total_runtime,
                                     morgan,
                                     migration)

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
                                               seed) {

  pop_for_cpp <- c()
  for (i in seq_along(individuals)) {
    x <- individuals[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- individuals[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }

  set.seed(seed)
  inbred_pop <- create_isofemale_line_cpp(pop_for_cpp, pop_size,
                                  total_runtime, morgan)

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



