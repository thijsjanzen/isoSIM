select_population <- function(source_pop,
                              select_matrix,
                              selection,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed,
                              write_to_file) {

  # first we have to convert source_pop to vector...
  pop_for_cpp <- c()
  for (i in seq_along(source_pop)) {
    x <- source_pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- source_pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }
  #then select_matrix to vector
  select <- as.vector(t(select_matrix))
  if (sum(is.na(select))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  selected_pop <- isoSIM::select_population_cpp(pop_for_cpp,
                                        select,
                                        selection,
                                        pop_size,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        write_to_file)

  selected_pop <- isoSIM::create_pop_class(selected_pop$population)

  return(selected_pop)
}

calculate_allele_frequencies <- function(source_pop,
                                         number_of_founders,
                                         step_size) {
  pop_for_cpp <- c()
  for (i in seq_along(source_pop)) {
    x <- source_pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- source_pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }

  frequency_table <- isoSIM::calculate_allele_spectrum_cpp(pop_for_cpp,
                                                           number_of_founders,
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
                                        seed,
                                        write_to_file = FALSE) {

  if (sum(is.na(select_matrix))) {
    stop("Can't start, there are NA values in the selection matrix!\n")
  }

  pop <- isoSIM::create_population_simulation_cpp(pop_size, number_of_founders, total_runtime, morgan,
                                                  select_matrix, selection, seed, write_to_file)
  popstruct <- isoSIM::create_pop_class(pop$population)
  return(popstruct)
}