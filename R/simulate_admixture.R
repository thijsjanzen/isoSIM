simulate_admixture <- function(input_population = NA,
                     pop_size = 100,
                     number_of_founders = 2,
                     initial_frequencies = NA,
                     total_runtime = 100,
                     morgan = 1,
                     seed,
                     select_matrix = NA,
                     progress_bar = TRUE,
                     track_junctions = FALSE,
                     track_frequency = FALSE,
                     multiplicative_selection = TRUE) {

  select <- select_matrix

  if(is.list(input_population)) {
    input_population <- population_to_vector(input_population)
  } else {
    input_population <- c(-1e6, -1e6)
  }

  if(sum(is.na(initial_frequencies))) {
    initial_frequencies = rep(1/number_of_founders, number_of_founders)
  }

  if(sum(initial_frequencies) != 1) {
    initial_frequencies = initial_frequencies / sum(initial_frequencies)
    cat("starting frequencies were normalized to 1\n")
  }

  no_selection <- FALSE
  if(is.na(select)) no_selection <- TRUE

  if(is.matrix(select)) {
    if (sum(is.na(select))) {
      stop("Can't start, there are NA values in the selection matrix!\n")
    }

    if (dim(select_matrix)[[2]] != 5) {
      stop("Incorrect dimensions of select_matrix,
           are you sure you provided all fitnesses?\n")
    }
  } else {
    if(is.na(select)) {
      select <- matrix(-1, nrow=2,ncol=2)
    }
  }

  markers <- c(-1,-1)

  if(is.matrix(select) & no_selection == FALSE) {
    markers <- (select[,1])
  }

  if(length(track_frequency) == 3)  {
    markers <- seq(track_frequency[1],
                   track_frequency[2],
                   length.out = track_frequency[3])

    if(is.matrix(select) && no_selection == FALSE) {
      markers <- c(markers, select[,1])
      markers <- sort(markers)
      markers <- unique(markers)
    }

    track_frequency <- TRUE
  }

  set.seed(seed)

  selected_pop <- simulate_cpp( input_population,
                                select,
                                pop_size,
                                number_of_founders,
                                initial_frequencies,
                                total_runtime,
                                morgan,
                                progress_bar,
                                track_frequency,
                                markers,
                                track_junctions,
                                multiplicative_selection)

  selected_popstruct <- create_pop_class(selected_pop$population)

  initial_freq_tibble <- create_tibble_from_freq_mat(selected_pop$initial_frequencies,
                                                     markers)

  final_freq_tibble   <- create_tibble_from_freq_mat(selected_pop$final_frequencies,
                                                     markers)

  output <- list()
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
                                                                 markers),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble)
  }

  if(track_frequency == TRUE && track_junctions == TRUE) {

    output <- list("population" = selected_popstruct,
                   "frequencies" = create_tibble_from_freq_table(selected_pop$frequencies,
                                                                 markers),
                   "initial_frequency" = initial_freq_tibble,
                   "final_frequency" = final_freq_tibble,
                   "junctions" = selected_pop$junctions)
  }

  return(output)
}

create_admixed_individuals <- function(num_individuals,
                           population_size,
                           number_of_founders,
                           size_in_morgan) {

  pop <- create_pop_admixed_cpp(num_individuals,
                                number_of_founders,
                                population_size,
                                size_in_morgan)

  popstruct <- create_pop_class(pop$population)

  output <- list("population" = popstruct)
  return(output)
}



