select_population <- function(source_pop,
                              selectMatrix,
                              selection,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed,
                              write_to_file) {
  
  # first we have to convert source_pop to vector...
  cat("R: converting pop\n")
  
  #for(i in seq_along(source_pop)) {
  #  source_pop[[i]]$chromosome1 <- t(source_pop[[i]]$chromosome1)
  #  source_pop[[i]]$chromosome2 <- t(source_pop[[i]]$chromosome2)
  #}
  
  #individuals <- unlist(t(source_pop))
  pop_for_cpp <- c()
  for (i in seq_along(source_pop)) {
    x <- source_pop[[i]]$chromosome1
    chrom1 <- as.vector(t(x))
    x <- source_pop[[i]]$chromosome2
    chrom2 <- as.vector(t(x))
    pop_for_cpp <- c(pop_for_cpp, chrom1, chrom2)
  }
  
  
  
  cat("R: converting select\n")
  #then selectMatrix to vector
  select <- as.vector(t(selectMatrix))
  
  
  
  cat("R: passing on to CPP\n")
  
  Sys.sleep(1)
  
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