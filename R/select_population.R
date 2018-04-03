population_to_vector <- function(source_pop) {
  pop_for_cpp <- c()
  for (i in seq_along(source_pop)) {
    x <- source_pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- source_pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }
  return(pop_for_cpp)
}

calculate_allele_frequencies <- function(source_pop,
                                         step_size,
                                         progress_bar = TRUE) {
  pop_for_cpp <- population_to_vector(source_pop)

  frequency_table <- calculate_allele_spectrum_cpp(pop_for_cpp,
                                                   step_size,
                                                   progress_bar)

  output <- tibble::as.tibble(frequency_table)
  colnames(output) <- c("location", "ancestor", "frequency")

  return(output)
}


create_population_selection <- function(pop_size,
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

  if (track_frequency == TRUE) {
    if (length(select_matrix[,1]) > 1) {
      stop("Can not track the frequency of more than one marker\n")
    }
  }

  if (dim(select_matrix)[[2]] != 5) {
    stop("Incorrect dimensions of select_matrix, are you sure you provided all fitnesses?\n")
  }

  set.seed(seed)
  pop <- create_population_selection_cpp(select_matrix,
                                                 pop_size,
                                                 number_of_founders,
                                                 total_runtime,
                                                 morgan,
                                                 progress_bar,
                                                 track_frequency)
  popstruct <- create_pop_class(pop$population)

  output <- list("population" = popstruct)

  if(track_frequency == TRUE) {
    found_markers <- c()
    for(i in 1:(dim(pop$frequencies)[[3]])) {
      local_mat <- pop$frequencies[,,i]
      time <- 0:(length(local_mat[,1])-1)
      marker_indicator <- rep(select_matrix[i, 1], length(time))
      freq_tibble <- tibble::as.tibble(cbind(time, marker_indicator, local_mat))
      colnames(freq_tibble) <- c("time", "location", 0:(length(local_mat[1,])-1))

      freq_tibble <- tidyr::gather(freq_tibble,
                                   key = "ancestor",
                                   value = "frequency",
                                   -c(1,2))

      found_markers <- rbind(found_markers, freq_tibble)
    }

    output <- list("population" = popstruct,
                   "frequencies" = found_markers)
  }

  return(output)
}

select_population <- function(source_pop,
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

  if (dim(select_matrix)[[2]] != 5) {
    stop("Incorrect dimensions of select_matrix, are you sure you provided all fitnesses?\n")
  }

  set.seed(seed)
  selected_pop <- select_population_cpp(pop_for_cpp,
                                        select,
                                        pop_size,
                                        total_runtime,
                                        morgan,
                                        progress_bar,
                                        track_frequency)

  selected_popstruct <- create_pop_class(selected_pop$population)

  output <- list("population" = selected_popstruct)

  if(track_frequency == TRUE) {
    found_markers <- c()
    for(i in 1:(dim(selected_pop$frequencies)[[3]])) {
      local_mat <- selected_pop$frequencies[,,i]
      time <- 0:(length(local_mat[,1])-1)
      marker_indicator <- rep(select_matrix[i, 1], length(time))
      freq_tibble <- tibble::as.tibble(cbind(time, marker_indicator, local_mat))
      colnames(freq_tibble) <- c("time", "location", 0:(length(local_mat[1,])-1))

      freq_tibble <- tidyr::gather(freq_tibble,
                                   key = "ancestor",
                                   value = "frequency",
                                   -c(1,2))

      found_markers <- rbind(found_markers, freq_tibble)
    }

    output <- list("population" = selected_popstruct,
                   "frequencies" = found_markers)
  }
  return(output)
}