simulate_admixture <- function(input_population = NA,
                     pop_size = 100,
                     number_of_founders = 2,
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


  if(is.matrix(select)) {
    if (sum(is.na(select))) {
      stop("Can't start, there are NA values in the selection matrix!\n")
    }

    if (dim(select_matrix)[[2]] != 5) {
      stop("Incorrect dimensions of select_matrix,
           are you sure you provided all fitnesses?\n")
    }
  }

  if(length(track_frequency) == 3)  {
   markers <- seq(track_frequency[1],
                   track_frequency[2],
                   length.out = track_frequency[3])

    to_add <- cbind(markers, -1, -1, -1, -1)
    if(is.matrix(select)) {
      select <- rbind(select, to_add)
    } else {
      select <- to_add
    }
    vx <- which(duplicated(select[,1]))
    # remove duplicate entries
    if(length(vx) > 0) {
      select <- select[-vx,]
    }

    track_frequency <- TRUE
  } else {
    select <- matrix(-1e6, -1, -1, -1, -1)
  }


  set.seed(seed)

  selected_pop <- simulate_cpp( input_population,
                                select,
                                pop_size,
                                number_of_founders,
                                total_runtime,
                                morgan,
                                progress_bar,
                                track_frequency,
                                track_junctions,
                                multiplicative_selection)

  selected_popstruct <- create_pop_class(selected_pop$population)

  initial_freq_tibble <- create_tibble_from_freq_mat(selected_pop$initial_frequencies,
                                                     select)
  final_freq_tibble   <- create_tibble_from_freq_mat(selected_pop$final_frequencies,
                                                     select)

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