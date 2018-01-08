select_population <- function(source_pop,
                              selectMatrix,
                              s,
                              pop_size,
                              total_runtime,
                              morgan,
                              seed) {
  
  # first we have to convert source_pop to vector...
  individuals <- c()
  for(j in 1:length(source_pop)) {
    for (i in seq_along(source_pop[[j]]$chromosome1[, 1])) {
      individuals <- c(individuals, source_pop[[1]]$chromosome1[i, ])
    }
    for (i in seq_along(source_pop[[j]]$chromosome2[, 1])) {
      individuals <- c(individuals, source_pop[[j]]$chromosome2[i, ])
    }
  }

  #then selectMatrix to vector
  select <- c();
  for(i in 1:length(selectMatrix[,1])) {
    select <- c(select, selectMatrix[i,])
  }
  
  selected_pop <- select_population_cpp(individuals,
                                        select,
                                        s,
                                        pop_size,
                                        total_runtime,
                                        morgan,
                                        seed)
  
  selected_pop <- isoSIM::create_pop_class(selected_pop$population)
  
  return(selected_pop)
}