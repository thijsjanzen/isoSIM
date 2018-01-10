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
                                              number_of_founders = 2,
                                              total_runtime = 100, 
                                              morgan = 1, 
                                              seed = 123, 
                                              write_to_file = FALSE)
  
  selectMatrix = matrix(ncol=3, nrow = 2)
  selectMatrix[1,] = c(0.0, 0.5, 0)
  selectMatrix[2,] = c(0.5, 1.0, 1)
  
  selected_pop <- isoSIM::select_population(sourcepop, selectMatrix,
                                            selection = 5,
                                            pop_size = 100,
                                            total_runtime = 10000,
                                            morgan = 1,
                                            seed = 1234,
                                            write_to_file = FALSE)
  
  freq_output <- calculate_allele_frequencies(selected_pop, 
                                                   number_of_founders = 2,
                                                   step_size = 0.01)
  
  require(ggplot2)
  ggplot(freq_output, aes(x = location, y = frequency, col = as.factor(ancestor))) +
    geom_line() + 
    ylim(c(0,1))
  
  a <- subset(freq_output, freq_output$location < 0.5)
  b <- a %>%
        group_by(as.factor(ancestor)) %>%
        summarise("mean_freq" = mean(frequency))
  
  a <- subset(freq_output, freq_output$location > 0.5)
  b <- a %>%
    group_by(as.factor(ancestor)) %>%
    summarise("mean_freq" = mean(frequency))
  
  
  sourcepop =  isoSIM::create_full_population(pop_size = 10000, 
                                              number_of_founders = 20,
                                              total_runtime = 1, 
                                              morgan = 1, 
                                              seed = 123, 
                                              write_to_file = FALSE)
  freq_output <- calculate_allele_frequencies(sourcepop, 
                                              number_of_founders = 20,
                                              step_size = 0.01)
  b <- freq_output %>%
    group_by(as.factor(ancestor)) %>%
    summarise("mean_freq" = mean(frequency))
  
  testthat::expect_equal(mean(b$mean_freq), 1/20, tolerance = 0.01)
  
  sourcepop =  isoSIM::create_full_population(pop_size = 1000, 
                                              number_of_founders = 5,
                                              total_runtime = 1000, 
                                              morgan = 1, 
                                              seed = 123, 
                                              write_to_file = FALSE)
  freq_output <- calculate_allele_frequencies(sourcepop, 
                                              number_of_founders = 5,
                                              step_size = 0.01)
  b <- freq_output %>%
    group_by(as.factor(ancestor)) %>%
    summarise("mean_freq" = mean(frequency))
  
  testthat::expect_equal(mean(b$mean_freq), 1/5, tolerance = 0.01)
  
  
  
  
})
  



