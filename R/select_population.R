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
  individuals <- c()
  for(j in 1:length(source_pop)) {
    for (i in seq_along(source_pop[[j]]$chromosome1[, 1])) {
      individuals <- c(individuals, source_pop[[j]]$chromosome1[i, ])
    }
    for (i in seq_along(source_pop[[j]]$chromosome2[, 1])) {
      individuals <- c(individuals, source_pop[[j]]$chromosome2[i, ])
    }
  }
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