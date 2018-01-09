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
  
  selected_pop <- isoSIM::select_population(sourcepop, selectMatrix,
                                    selection = 0.1,
                                    pop_size = 100,
                                    total_runtime = 100,
                                    morgan = 1,
                                    seed = 1234,
                                    write_to_file = FALSE)
  

  testthat::expect_equal(length(selected_pop), 100)
  
})

test_that("allele frequencies", {
  
  sourcepop =  isoSIM::create_full_population(pop_size = 100, 
                                              number_of_founders = 10,
                                              total_runtime = 1000, 
                                              morgan = 1, 
                                              seed = 123, 
                                              write_to_file = FALSE)
  
  selectMatrix = matrix(ncol=3, nrow = 2)
  selectMatrix[1,] = c(0.05, 0.1, 0)
  selectMatrix[2,] = c(0.15, 0.5, 1)
  
  selected_pop <- isoSIM::select_population(sourcepop, selectMatrix,
                                            selection = 0.1,
                                            pop_size = 100,
                                            total_runtime = 100,
                                            morgan = 1,
                                            seed = 1234,
                                            write_to_file = FALSE)
  
  freq_output <- calculate_allele_frequencies(selected_pop, 
                                                   number_of_founders = 10,
                                                   step_size = 0.01)
  
}
  



