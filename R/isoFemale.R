create_isoFemaleLine <- function(isoFemale, pop_size = 100, 
                                 run_time = 100, morgan = 1, SEED = 42) {
  
  #we have to convert the individual back to a vector
  indiv <- c()
  for(i in 1:length(isoFemale$chromosome1[,1])) {
    indiv <- c(indiv,isoFemale$chromosome1[i,])
  }
  for(i in 1:length(isoFemale$chromosome2[,1])) {
    indiv <- c(indiv,isoFemale$chromosome2[i,])
  }
  
  inbred_pop <- create_femaleLine(indiv, pop_size, 
                              run_time, morgan, 
                              SEED);
  
  inbred_population <- create_pop_class(inbred_pop$population)
  
  output <- inbred_population[[sample(1:length(inbred_population),1)]]
  
  return(output);
}


create_isoFemale <- function(pop, n = 1, simulate = FALSE, 
                             run_time = 1000, morgan = 1) {
  if(simulate == FALSE) {
    #just pick n individuals, and make them homozygous
    indices <- sample(1:length(pop), n)
    isoFemales <- pop[indices]
    
    for(i in 1:n) {
      isoFemales[i]$chromosome2 = isoFemales[i]$chromosome1
    }
    class(isoFemales) <- "population"
    return(isoFemales)
  }
  
  if(simulate == TRUE) {
    #for each isofemale
    #pick the female, then simulate until run_time
    #then pick a random individual (assuming the population is fixe now)
    
    indices <- sample(1:length(pop$Population), n)
    isoFemales <- pop$Population[indices]
    output_females <- c()
    for(i in 1:n) {
      vx <- create_isoFemaleLine(isoFemales[[i]], pop_size = 100, 
                                 run_time = run_time, 
                                 morgan = morgan, 
                                 SEED = i)
      output_females[[i]] <- vx
    }
    class(output_females) <- "population"
    return(output_females)
  }
}