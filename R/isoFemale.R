create_iso_female <- function(source_pop,
                             n = 1,
                             inbreeding_pop_size = 100,
                             run_time = 1000,
                             morgan = 1,
                             seed = 42) {

  if(n == 1) {
    indices <- sample(1:length(source_pop), 2, replace = FALSE)
    parents = list(source_pop[[indices[1]]], source_pop[[indices[2]]])
    inbred_population <- create_population_from_individuals(parents,
                                                            pop_size,
                                                            run_time,
                                                            morgan,
                                                            seed)
    output_females <- list()
    output_females[[1]] <- inbred_population[[sample(1:length(inbred_population), 1)]]
    class(output_females[[1]]) <- "individual"
    return(output_females)
  }



  # first we select the individuals that will be the parents of the isofemales
  indices <- sample(seq_along(source_pop), n * 2, replace = FALSE)

  #for each isofemale
  #pick the female, then simulate until run_time
  #then pick a random individual (assuming the population is fixed now)

  iso_females <- source_pop[indices]
  output_females <- list()
  for (i in 1:n) {
    parents <- list(iso_females[[i]], iso_females[[i + n]])

    inbred_population <- create_population_from_individuals(parents,
                                             pop_size,
                                             run_time,
                                             morgan,
                                             seed + i)
    output_females[[i]] <- inbred_population[[sample(1:length(inbred_population), 1)]]
    class(output_females[[i]]) <- "individual"
  }
  return(output_females)
}