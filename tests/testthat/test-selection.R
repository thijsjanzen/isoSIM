context("selection")

test_that("select on population", {
  
  sourcepop =  isoSIM::create_full_population(pop_size = 100, 
                                    number_of_founders = 10,
                                    total_runtime = 1000, 
                                    morgan = 1, 
                                    seed = 123, 
                                    write_to_file = FALSE)
  
  selectMatrix = matrix(ncol=3, nrow = 2)
  selectMatrix[1,] = c(0.05, 0.1, 0)
  selectMatrix[2,] = c(0.15, 0.5, 1)
  
  
#  require("debug")
 # mtrace(select_population_cpp)
  
  selected_pop <- isoSIM::select_population(sourcepop, selectMatrix,
                                    selection = 0.1,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1234,
                                    write_to_file = FALSE)
  
#  mtrace(select_population_cp, FALSE)
  
  
  testthat::expect_equal(length(selected_pop), 100)
  
  if(1 == 2) {
    source_pop <- sourcepop
    selection = 0.1
    pop_size = 100
    total_runtime = 100
    morgan = 1
    seed = 1234
    write_to_file = FALSE
  
    
    v1 <- pop_for_cpp
    selectM <- select
    population_size <- pop_size
    run_time <- total_runtime
    writeToFile <- write_to_file
    #selectM, s, population_size, run_time, morgan, seed, writeToFile  
  }
  
  
  
})