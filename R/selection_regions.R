
select_population <- function(source_pop,
                              select_matrix,
                              selection,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed,
                              progress_bar = TRUE) {

  # first we have to convert source_pop to vector...
  pop_for_cpp <- population_to_vector(source_pop)
  #then select_matrix to vector
  #select <- as.vector(t(select_matrix))
  select <- select_matrix
  if (sum(is.na(select))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  set.seed(seed)
  selected_pop <- select_population_cpp(pop_for_cpp,
                                        select,
                                        selection,
                                        pop_size,
                                        total_runtime,
                                        morgan,
                                        progress_bar)

  selected_pop <- create_pop_class(selected_pop$population)

  return(selected_pop)
}

create_population_selection <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        select_matrix,
                                        selection,
                                        seed,
                                        progress_bar = TRUE) {

  if (sum(is.na(select_matrix))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  if (sum(select_matrix[,2] < select_matrix[,1])) {
    stop("Can't start, select matrix incorrect format\n")
  }

  set.seed(seed)
  pop <- create_population_selection_cpp(pop_size,
                                         number_of_founders,
                                         total_runtime,
                                         morgan,
                                         select_matrix,
                                         selection,
                                         progress_bar)

  popstruct <- create_pop_class(pop$population)
  return(popstruct)
}