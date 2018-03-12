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

select_population <- function(source_pop,
                              select_matrix,
                              selection,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed) {

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
                                        morgan)

  selected_pop <- create_pop_class(selected_pop$population)

  return(selected_pop)
}

calculate_allele_frequencies <- function(source_pop,
                                         step_size) {
  pop_for_cpp <- population_to_vector(source_pop)

  frequency_table <- calculate_allele_spectrum_cpp(pop_for_cpp,
                                                   step_size)

  output <- tibble::as.tibble(frequency_table)
  colnames(output) <- c("location", "ancestor", "frequency")

  return(output)
}

create_population_selection <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        select_matrix,
                                        selection,
                                        seed) {

  if (sum(is.na(select_matrix))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  if (sum(select_matrix[,2] < select_matrix[,1])) {
    stop("Can't start, select matrix incorrect format\n")
  }



  set.seed(seed)
  pop <- create_population_selection_cpp(pop_size, number_of_founders, total_runtime, morgan,
                                                  select_matrix, selection)
  popstruct <- create_pop_class(pop$population)
  return(popstruct)
}

create_population_selection_markers <- function(pop_size,
                                        number_of_founders,
                                        total_runtime,
                                        morgan,
                                        select_matrix,
                                        seed) {

  if (sum(is.na(select_matrix))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  set.seed(seed)
  pop <- create_population_selection_markers_cpp(select_matrix, pop_size,
                                       number_of_founders, total_runtime, morgan)
  popstruct <- create_pop_class(pop$population)
  return(popstruct)
}

select_population_markers <- function(source_pop,
                              select_matrix,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed) {

  # first we have to convert source_pop to vector...
  pop_for_cpp <- population_to_vector(source_pop)
  #then select_matrix to vector
  #select <- as.vector(t(select_matrix))
  select <- select_matrix
  if (sum(is.na(select))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  set.seed(seed)
  selected_pop <- select_population_markers_cpp(pop_for_cpp,
                                        select,
                                        pop_size,
                                        total_runtime,
                                        morgan)

  selected_pop <- create_pop_class(selected_pop$population)

  return(selected_pop)
}


