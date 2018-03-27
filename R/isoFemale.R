create_iso_female <- function(source_pop,
                             n = 1,
                             inbreeding_pop_size = 100,
                             run_time = 2000,
                             morgan = 1,
                             seed = 42,
                             progress_bar = TRUE) {

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
                                                            inbreeding_pop_size,
                                             run_time,
                                             morgan,
                                             seed + i)
    output_females[[i]] <- inbred_population[[sample(seq_along(inbred_population), 1)]]
    class(output_females[[i]]) <- "individual"
  }
  return(output_females)
}