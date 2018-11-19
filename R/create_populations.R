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

  pop1 <- create_pop_class(pops$population_1)
  pop2 <- create_pop_class(pops$population_2)

  output <- list("Population_1" = pop1,
                 "Population_2" = pop2)
  return(output)
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