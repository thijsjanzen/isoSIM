select_population <- function(source_pop,
                              selectMatrix,
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
  
  
  
  #then selectMatrix to vector
  select <- as.vector(t(selectMatrix))
  
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
  
  output <- as.data.frame("location" = frequency_table[,1],
                          "ancestor" = frequency_table[,2],
                          "frequency" = frequency_table[,3])
  
  return(output)
}