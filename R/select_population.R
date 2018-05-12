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

create_tibble_from_freq_table <- function(frequencies, select_matrix) {
  input_list <- list()
  for(i in 1:(dim(frequencies)[[3]])) {
    local_mat <- frequencies[ , , i]
    input_list[[i]] <- local_mat
  }

  to_apply <- function(local_mat, i) {
    time <- 0:(length(local_mat[[i]][,1])-1)
    marker_indicator <- rep(select_matrix[i, 1], length(time))
    freq_tibble <- tibble::as.tibble( cbind(time,
                                            marker_indicator,
                                            local_mat[[i]]))
    colnames(freq_tibble) <- c("time",
                               "location",
                               0:(length(local_mat[[i]][1,]) - 1))

    freq_tibble <- tidyr::gather(freq_tibble,
                                 key = "ancestor",
                                 value = "frequency",
                                 -c(1,2))
    return(as.data.frame(freq_tibble))
  }

  interm_list <- lapply(seq_along(input_list), to_apply, local_mat = input_list)
  vy <- tibble::as.tibble(dplyr::bind_rows(interm_list))

  return(vy)
}

create_tibble_from_freq_mat <- function(frequencies, select_matrix) {
  found_markers <- c()
  for(i in 1:(dim(frequencies)[[1]])) {
    local_mat <- frequencies[i,]
    time <- 0
    marker_indicator <- select_matrix[i, 1]
    freq_tibble <- c(time, marker_indicator, local_mat)

    found_markers <- rbind(found_markers, freq_tibble)
  }
  colnames(found_markers) <- c("time",
                               "location",
                               0:(length(frequencies[1,])-1))
  found_markers <- tibble::as.tibble(found_markers)
  found_markers <- tidyr::gather(found_markers,
                                 key = "ancestor",
                                 value = "frequency",
                                  -c(1,2))

  return(found_markers)
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

  if (dim(select_matrix)[[2]] != 5) {
    stop("Incorrect dimensions of select_matrix,
         are you sure you provided all fitnesses?\n")
  }

  if(length(track_frequency) == 3)  {
    markers <- seq(track_frequency[1],
                   track_frequency[2],
                   length.out = track_frequency[3])

    to_add <- cbind(markers, 1, 1, 1, -1)
    select_matrix <- rbind(select_matrix, to_add)
    vx <- which(duplicated(select_matrix[,1]))
    # remove duplicate entries
    if(length(vx) > 0) {
      select_matrix <- select_matrix[-vx,]
    }

    track_frequency <- TRUE
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

  if (dim(select_matrix)[[2]] != 5) {
    stop("Incorrect dimensions of select_matrix, are you sure you provided all fitnesses?\n")
  }

  if(length(track_frequency) == 3)  {
    markers <- seq(track_frequency[1],
                   track_frequency[2],
                   length.out = track_frequency[3])

    to_add <- cbind(markers, 1, 1, 1, -1)
    select <- rbind(select, to_add)
    vx <- which(duplicated(select[,1]))
    # remove duplicate entries
    if(length(vx) > 0) {
      select <- select[-vx,]
    }

    track_frequency <- TRUE
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

  initial_freq_tibble <- create_tibble_from_freq_mat(selected_pop$initial_frequencies,
                                                     select)
  final_freq_tibble <- create_tibble_from_freq_mat(selected_pop$final_frequencies,
                                                   select)

  output <- list("population" = selected_popstruct,
                 "initial_frequency" = initial_freq_tibble,
                 "final_frequency" = final_freq_tibble)

  if(track_frequency == TRUE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(selected_pop$frequencies,
                                                                 select),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }
  return(output)
}

increase_ancestor <- function(population, increment = 20) {
  increase_indiv <- function(indiv) {

    # -1 indicates the end, no increment there!
    positive <- which(indiv$chromosome1[,2] > -1)
    indiv$chromosome1[positive, 2] <- indiv$chromosome1[positive, 2] + increment

    # -1 indicates the end, no increment there!
    positive <- which(indiv$chromosome2[,2] > -1)
    indiv$chromosome2[positive, 2] <- indiv$chromosome2[positive, 2] + increment
    return(indiv)
  }

  pop_2 <- lapply(population, increase_indiv)
  class(pop_2) <- "population"
  return(pop_2)
}