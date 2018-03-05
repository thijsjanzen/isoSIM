create_iso_female_line <- function(parents,
                                 pop_size = 100,
                                 run_time = 100,
                                 morgan = 1,
                                 seed = 42) {

  #we have to convert the parents back to a vector
  indiv <- c()
  for (i in seq_along(parents[[1]]$chromosome1[, 1])) {
    indiv <- c(indiv, parents[[1]]$chromosome1[i, ])
  }
  for (i in seq_along(parents[[1]]$chromosome2[, 1])) {
    indiv <- c(indiv, parents[[1]]$chromosome2[i, ])
  }

  for (i in seq_along(parents[[2]]$chromosome1[, 1])) {
    indiv <- c(indiv, parents[[2]]$chromosome1[i, ])
  }
  for (i in seq_along(parents[[2]]$chromosome2[, 1])) {
    indiv <- c(indiv, parents[[2]]$chromosome2[i, ])
  }

  set.seed(seed)
  inbred_pop <- create_isofemale_line_cpp(indiv, pop_size,
                                  run_time, morgan)

  inbred_population <- isoSIM::create_pop_class(inbred_pop$population)

  if (length(inbred_population) < 1) {
    stop("creating isofemale failed\n")
  }

  if (length(inbred_population) == 1) {
    return(inbred_population)
  }

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