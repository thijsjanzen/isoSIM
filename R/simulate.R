simulate <- function(input_population = NA,
                     pop_size,
                     number_of_founders,
                     total_runtime,
                     morgan,
                     seed,
                     select_matrix = NA,
                     progress_bar = TRUE,
                     track_junctions = FALSE,
                     track_frequency = FALSE) {

  if(!is.list(input_population)) {
    if(!is.matrix(select_matrix)) {
      return(create_population(pop_size, number_of_founders, total_runtime,
                             morgan, seed, progress_bar, track_junctions))
    }
    if(is.matrix(select_matrix)) {
      return(create_population_selection(pop_size, number_of_founders,
                                         total_runtime, morgan,
                                         select_matrix, seed,
                                         track_frequency, progress_bar))
    }
  }

  if(is.list(input_population)) {
    if(!is.matrix(select_matrix)) {
      return(create_population_from_individuals(input_population, pop_size, total_runtime,
                                              morgan, seed, progress_bar))
    }
    if(is.matrix(select_matrix)) {
      return(select_population(input_population, select_matrix, pop_size,
                               total_runtime, morgan, seed,
                               track_frequency, progress_bar))
    }
  }
}