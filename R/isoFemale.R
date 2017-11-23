create_iso_female_line <- function(parents,
                                 pop_size = 100,
                                 run_time = 100,
                                 morgan = 1,
                                 seed = 42) {

  #we have to convert the parents back to a vector
  indiv <- c()
  for (i in 1:length(parents[[1]]$chromosome1[, 1])) {
    indiv <- c(indiv, parents[[1]]$chromosome1[i, ])
  }
  for (i in 1:length(parents[[1]]$chromosome2[, 1])) {
    indiv <- c(indiv, parents[[1]]$chromosome2[i, ])
  }

  for (i in 1:length(parents[[2]]$chromosome1[, 1])) {
    indiv <- c(indiv, parents[[2]]$chromosome1[i, ])
  }
  for (i in 1:length(parents[[2]]$chromosome2[, 1])) {
    indiv <- c(indiv, parents[[2]]$chromosome2[i, ])
  }

  inbred_pop <- create_femaleLine(indiv, pop_size,
                                  run_time, morgan,
                                  seed);

  inbred_population <- isoSIM::create_pop_class(inbred_pop$population)

  output <- inbred_population[[sample(1:length(inbred_population), 1)]]

  return(output)
}


create_iso_female <- function(pop,
                             n = 1,
                             run_time = 1000,
                             morgan = 1) {

  # first we select the individuals that will be the parents of the isofemales
  indices <- sample(1:length(pop), n * 2, replace = FALSE)

  #for each isofemale
  #pick the female, then simulate until run_time
  #then pick a random individual (assuming the population is fixed now)

  iso_females <- pop[indices]
  output_females <- c()
  for (i in 1:n) {
    parents <- list(iso_females[[i]], iso_females[[i + 1]])

    vx <- create_iso_female_line(parents,
                                 pop_size = length(pop),
                                 run_time = run_time,
                                 morgan = morgan,
                                 seed = i)
    output_females[[i]] <- vx
  }
  class(output_females) <- "population"
  return(output_females)
}