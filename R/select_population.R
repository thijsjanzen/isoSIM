select_population <- function(source_pop,
                              selectMatrix,
                              s,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed,
                              write_to_file) {
  
  # first we have to convert source_pop to vector...
  cat("R: converting pop\n")
  
  for(i in seq_along(source_pop)) {
    source_pop[[i]]$chromosome1 <- t(source_pop[[i]]$chromosome1)
    source_pop[[i]]$chromosome2 <- t(source_pop[[i]]$chromosome2)
  }
  
  individuals <- unlist(t(source_pop))
  
  
  cat("R: converting select\n")
  #then selectMatrix to vector
  select <- c();
  for(i in 1:length(selectMatrix[,1])) {
    select <- c(select, selectMatrix[i,])
  }
  
  cat("R: passing on to CPP\n")
  selected_pop <- select_population_cpp(individuals,
                                        select,
                                        s,
                                        pop_size,
                                        total_runtime,
                                        morgan,
                                        seed,
                                        write_to_file)
  
  selected_pop <- isoSIM::create_pop_class(selected_pop$population)
  
  return(selected_pop)
}