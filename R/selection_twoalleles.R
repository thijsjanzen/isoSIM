create_population_selection_twoalleles <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        select_matrix,
                                        seed,
                                        track_frequency = FALSE,
                                        progress_bar = TRUE) {

  if (sum(is.na(select_matrix))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  if(track_frequency == TRUE) {
    if(length(select_matrix[,1]) > 1) {
      stop("Can not track the frequency of more than one marker\n")
    }
  }

  set.seed(seed)
  pop <- create_population_selection_twoalleles_cpp(select_matrix,
                                                 pop_size,
                                                 number_of_founders,
                                                 total_runtime,
                                                 morgan,
                                                 progress_bar,
                                                 track_frequency)
  popstruct <- create_pop_class(pop$population)

  output <- list("population" = popstruct)

  if(track_frequency == TRUE) {
    time <- 0:(length(pop$frequencies[,1])-1)
    freq_tibble <- tibble::as.tibble(cbind(time,pop$frequencies))
    colnames(freq_tibble) <- c("time", 0:(length(pop$frequencies[1,])-1))

    freq_tibble <- tidyr::gather(freq_tibble,
                                 key = "ancestor",
                                 value = "frequency",
                                 -1)

    output <- list("population" = popstruct,
                   "frequencies" = freq_tibble)
  }

  return(output)
}

select_population_twoalleles <- function(source_pop,
                              select_matrix,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed,
                              track_frequency = FALSE,
                              progress_bar = TRUE) {

  # first we have to convert source_pop to vector...
  pop_for_cpp <- population_to_vector(source_pop)
  #then select_matrix to vector
  #select <- as.vector(t(select_matrix))
  select <- select_matrix
  if (sum(is.na(select))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  if(track_frequency == TRUE) {
    if(length(select_matrix[,1]) > 1) {
      stop("Can not track the frequency of more than one marker\n")
    }
  }

  set.seed(seed)
  selected_pop <- select_population_twoalleles_cpp(pop_for_cpp,
                                                select,
                                                pop_size,
                                                total_runtime,
                                                morgan,
                                                progress_bar,
                                                track_frequency)

  selected_popstruct <- create_pop_class(selected_pop$population)

  output <- list("population" = selected_popstruct)

  if(track_frequency == TRUE) {
    time <- 0:(length(selected_pop$frequencies[, 1]) - 1)
    freq_tibble <- tibble::as.tibble(cbind(time,
                                           selected_pop$frequencies))
    colnames(freq_tibble) <- c("time", 0:(
                                    length(selected_pop$frequencies[1, ]) - 1
                                        )
                               )

    freq_tibble <- tidyr::gather(freq_tibble,
                                 key = "ancestor",
                                 value = "frequency",
                                 -1)

    output <- list("population" = selected_popstruct,
                   "frequencies" = freq_tibble)
  }

  return(output)
}