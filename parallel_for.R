library(isoSIM)
library(foreach)

library(doParallel)

library(isoSIM)
library(tidyverse)

assess_match_R <- function(chrom, start, end, ancestor) {
  
  if(is.null(chrom)) {
    cat("error! chromosome empty\n")
  }
  
  block <- c()
  for(i in 1:length(chrom[, 1])) {
     if(chrom[i, 1][[1]] > end) break;
    
    if(chrom[i, 1][[1]] > start) {
      block <- rbind(block, chrom[i, ])
    } else {
      if(i + 1 < length(chrom[, 1])) {
        if(chrom[i+1, 1][[1]] > start) {
          block <- rbind(block, chrom[i+1, ])
        }
      }
    }
  }
  
  if(length(block[,1]) == 1) {
    if(block[1, 2] == ancestor) {
      return(1.0)
    } else {
      return(0.0)
    }
  }
  
  match = 0.0
  if(is.null(block)) return(0.0)
  
  for(i in 1:length(block[,1])) {
    if(block[i, 2] == ancestor) {
      local_right <- end
      if(i + 1 < length(block[,1])) {
        local_right <- block[i+1, 1]
      }
      local_left <- block[i, 1]
      if(local_left < start) local_left <- start
      
      match <- match + (local_right - local_left)
    }
  }
  match <- match * 1.0 / (end - start)
  return(match)
}




calculate_in_R <- function(pop, 
                           number_of_founders,
                           step_size) {
  
  assess_match_R <- function(chrom, start, end, ancestor) {
    
    if(is.null(chrom)) {
      cat("error! chromosome empty\n")
    }
    
    block <- c()
    for(i in 1:length(chrom[, 1])) {
      if(chrom[i, 1][[1]] > end) break;
      
      if(chrom[i, 1][[1]] > start) {
        block <- rbind(block, chrom[i, ])
      } else {
        if(i + 1 < length(chrom[, 1])) {
          if(chrom[i+1, 1][[1]] > start) {
            block <- rbind(block, chrom[i+1, ])
          }
        }
      }
    }
    
    if(length(block[,1]) == 1) {
      if(block[1, 2] == ancestor) {
        return(1.0)
      } else {
        return(0.0)
      }
    }
    
    match = 0.0
    if(is.null(block)) return(0.0)
    
    for(i in 1:length(block[,1])) {
      if(block[i, 2] == ancestor) {
        local_right <- end
        if(i + 1 < length(block[,1])) {
          local_right <- block[i+1, 1]
        }
        local_left <- block[i, 1]
        if(local_left < start) local_left <- start
        
        match <- match + (local_right - local_left)
      }
    }
    match <- match * 1.0 / (end - start)
    return(match)
  }
  
  grid <- seq(step_size, 1, by = step_size)
  output <- matrix(ncol = 3, nrow = 1 + number_of_founders * length(grid))
  left <- 0;
  pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
  for(j in seq_along(grid)) {
    step <- grid[j]
    right <- step
    
    for(a in 1:number_of_founders) {
      toAdd <- c(right, a, 0)
      for(i in seq_along(pop)) {
          a1 <- assess_match_R(pop[[i]]$chromosome1, left, right, a)
          a2 <- assess_match_R(pop[[i]]$chromosome2, left, right, a)
          toAdd[3] <- toAdd[3] + a1 + a2  
      }
      index <- (a-1) * length(grid) + j
      output[index,] <- toAdd
    }
    left <- right
    setTxtProgressBar(pb, j)
  }
  return(output)
}


calculate_in_R_parallel <- function(pop, 
                           number_of_founders,
                           step_size) {
  
  assess_match_R <- function(chrom, start, end, ancestor) {
    
    if(is.null(chrom)) {
      cat("error! chromosome empty\n")
    }
    
    block <- c()
    for(i in 1:length(chrom[, 1])) {
      if(chrom[i, 1][[1]] > end) break;
      
      if(chrom[i, 1][[1]] > start) {
        block <- rbind(block, chrom[i, ])
      } else {
        if(i + 1 < length(chrom[, 1])) {
          if(chrom[i+1, 1][[1]] > start) {
            block <- rbind(block, chrom[i+1, ])
          }
        }
      }
    }
    
    if(length(block[,1]) == 1) {
      if(block[1, 2] == ancestor) {
        return(1.0)
      } else {
        return(0.0)
      }
    }
    
    match = 0.0
    if(is.null(block)) return(0.0)
    
    for(i in 1:length(block[,1])) {
      if(block[i, 2] == ancestor) {
        local_right <- end
        if(i + 1 < length(block[,1])) {
          local_right <- block[i+1, 1]
        }
        local_left <- block[i, 1]
        if(local_left < start) local_left <- start
        
        match <- match + (local_right - local_left)
      }
    }
    match <- match * 1.0 / (end - start)
    return(match)
  }
  
  grid <- seq(step_size, 1, by = step_size)
  output <- matrix(ncol = 3, nrow = 1 + number_of_founders * length(grid))
  left <- 0;
  pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
  
  cl <- makeCluster(6)
  registerDoParallel(cl)
  
  
  foreach(j = 1:length(grid)) %dopar%
  {
    step <- grid[j]
    right <- step
    
    for(a in 1:number_of_founders) {
      toAdd <- c(right, a, 0)
      for(i in seq_along(pop)) {
        a1 <- assess_match_R(pop[[i]]$chromosome1, left, right, a)
        a2 <- assess_match_R(pop[[i]]$chromosome2, left, right, a)
        toAdd[3] <- toAdd[3] + a1 + a2  
      }
      index <- (a-1) * length(grid) + j
      output[index,] <- toAdd
    }
    left <- right
    setTxtProgressBar(pb, j)
  }
  
  stopCluster(cl)
  return(output)
}


number_founders <- 20
sourcepop =  isoSIM::create_full_population(pop_size = 100, 
                                            number_of_founders = number_founders,
                                            total_runtime = 1, 
                                            morgan = 1, 
                                            seed = 123, 
                                            write_to_file = FALSE)

selectMatrix = matrix(ncol=3, nrow = 1)

under_selection <- 1


selectMatrix[1,] = c(0.4, 0.401, under_selection)


selected_pop <- isoSIM::select_population(sourcepop, selectMatrix,
                                          selection = 0.01,
                                          pop_size = 100,
                                          total_runtime = 200,
                                          morgan = 1,
                                          seed = 12345,
                                          write_to_file = FALSE)

#system.time(
#  freq_output <- calculate_allele_frequencies(selected_pop, 
                                              number_of_founders = number_founders,
                                              step_size = 0.01)
#)

#system.time(
#  freq_output2 <- calculate_in_R(selected_pop, number_founders, 0.01)
  
#  freq_output3 <- calculate_in_R_parallel(selected_pop, number_founders, 0.01)
#)

require(microbenchmark)  
  
  vx <- microbenchmark(
    "C++" = calculate_allele_frequencies(selected_pop, 
                                                 number_of_founders = number_founders,
                                                 step_size = 0.001),
    #"R_serial" = calculate_in_R(selected_pop, number_founders, 0.001),
    "R_parallel" = calculate_in_R_parallel(selected_pop, number_founders, 0.001),
    times = 3
  )

  autoplot(vx)




