
create_population_selection_region <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        select_matrix,
                                        seed,
                                        track_frequency,
                                        progress_bar = TRUE) {

  if (sum(is.na(select_matrix))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  if (dim(select_matrix)[[2]] != 5) {
    stop("Incorrect dimensions of select_matrix, are you sure you provided all fitnesses?\n")
  }

  # we assume that track_frequency is a vector with three numbers:
  # start, end, number of markers
  # we append these to the select_matrix
  markers <- seq(track_frequency[1],
                 track_frequency[2],
                 length.out = track_frequency[3])

  to_add <- cbind(markers, 1, 1, 1, -1)
  select_matrix <- rbind(select_matrix, to_add)
  vx <- which(duplicated(select_matrix[,1]))
  # remove duplicate entries
  select_matrix <- select_matrix[-vx,]

  track_frequency <- TRUE

  set.seed(seed)
  pop <- create_population_selection_cpp(select_matrix,
                                         pop_size,
                                         number_of_founders,
                                         total_runtime,
                                         morgan,
                                         progress_bar,
                                         track_frequency)
  popstruct <- create_pop_class(pop$population)

  initial_freq_tibble <- create_tibble_from_freq_mat(pop$initial_frequencies,
                                                     select_matrix)
  final_freq_tibble <- create_tibble_from_freq_mat(pop$final_frequencies,
                                                   select_matrix)


  output <- list("population" = popstruct,
                 "initial_frequency" = initial_freq_tibble,
                 "final_frequency" = final_freq_tibble)

  if(track_frequency == TRUE) {

    output <- list("population" = popstruct,
                   "frequencies" = create_tibble_from_freq_table(pop$frequencies,
                                                                 select_matrix),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  return(output)
}

select_population_region <- function(source_pop,
                              select_matrix,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed,
                              track_frequency,
                              progress_bar = TRUE) {

  # first we have to convert source_pop to vector...
  pop_for_cpp <- population_to_vector(source_pop)
  #then select_matrix to vector
  #select <- as.vector(t(select_matrix))
  select <- select_matrix
  if (sum(is.na(select))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  if (dim(select_matrix)[[2]] != 5) {
    stop("Incorrect dimensions of select_matrix, are you sure you provided all fitnesses?\n")
  }

  # we assume that track_frequency is a vector with three numbers:
  # start, end, number of markers
  # we append these to the select_matrix
  markers <- seq(track_frequency[1],
                 track_frequency[2],
                 length.out = track_frequency[3])

  to_add <- cbind(markers, 1, 1, 1, -1)
  select_matrix <- rbind(select_matrix, to_add)
  vx <- which(duplicated(select_matrix[,1]))
  # remove duplicate entries
  select_matrix <- select_matrix[-vx,]

  track_frequency <- TRUE

  set.seed(seed)
  selected_pop <- select_population_cpp(pop_for_cpp,
                                        select,
                                        pop_size,
                                        total_runtime,
                                        morgan,
                                        progress_bar,
                                        track_frequency)

  selected_popstruct <- create_pop_class(selected_pop$population)

  initial_freq_tibble <- create_tibble_from_freq_mat(selected_pop$initial_frequencies,
                                                     select_matrix)
  final_freq_tibble <- create_tibble_from_freq_mat(selected_pop$final_frequencies,
                                                   select_matrix)

  output <- list("population" = selected_popstruct,
                 "initial_frequency" = initial_freq_tibble,
                 "final_frequency" = final_freq_tibble)

  if(track_frequency == TRUE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(selected_pop$frequencies,
                                                                 select_matrix),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }
  return(output)
}