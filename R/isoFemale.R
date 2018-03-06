create_iso_female_line <- function(parents,
                                 pop_size = 100,
                                 run_time = 100,
                                 morgan = 1,
                                 seed = 42) {

  inbred_population <- create_population_from_individuals(parents,
                                                          pop_size,
                                                          run_time,
                                                          morgan,
                                                          seed)

  output <- inbred_population[[sample(1:length(inbred_population), 1)]]

  return(output)
}

create_iso_female <- function(source_pop,
                             n = 1,
                             inbreeding_pop_size = 100,
                             run_time = 1000,
                             morgan = 1) {

  # first we select the individuals that will be the parents of the isofemales
  indices <- sample(seq_along(source_pop), n * 2, replace = FALSE)

  #for each isofemale
  #pick the female, then simulate until run_time
  #then pick a random individual (assuming the population is fixed now)

  iso_females <- source_pop[indices]
  output_females <- c()
  for (i in seq_len(n)) {
    parents <- list(iso_females[[i]], iso_females[[i + n]])

    set.seed(i)
    vx <- create_iso_female_line(parents,
                                 pop_size = inbreeding_pop_size,
                                 run_time = run_time,
                                 morgan = morgan)
    output_females[[i]] <- vx
    class(output_females[[i]]) <- "individual"
  }
  return(output_females)
}