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

  selected_pop <- simulate_cpp(input_population, select_matrix,
               pop_size, number_of_founders, total_runtime,
               morgan, progress_bar, track_frequency, track_junctions)

  selected_popstruct <- create_pop_class(selected_pop$population)

  initial_freq_tibble <- create_tibble_from_freq_mat(selected_pop$initial_frequencies,
                                                     select)
  final_freq_tibble <- create_tibble_from_freq_mat(selected_pop$final_frequencies,
                                                   select)



  if(track_frequency == FALSE && track_junctions == FALSE) {
    output <- list("population" = selected_popstruct,
                 "initial_frequency" = initial_freq_tibble,
                 "final_frequency" = final_freq_tibble,
                 "junctions" = selected_pop$junctions)
  }

  if(track_frequency == FALSE && track_junctions == TRUE) {
    output <- list("population" = selected_popstruct,
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if(track_frequency == TRUE && track_junctions == FALSE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(selected_pop$frequencies,
                                                                 select),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if(track_frequency == TRUE && track_junctions == TRUE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(selected_pop$frequencies,
                                                                 select),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}

simulate_2 <- function(input_population = NA,
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