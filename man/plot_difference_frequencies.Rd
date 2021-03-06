\name{plot_difference_frequencies}
\alias{plot_difference_frequencies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot the change in frequency between the start and end of a simulation
}
\description{
 This function plots the change in frequency of one or multiple ancestors after performing a simulation.
}
\usage{
plot_difference_frequencies(results, picked_ancestor = "ALL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
    An object which is the result of \code{select_population} or \code{create_population_selection}, being a list with four properties: \code{population}, \code{frequencies}, \code{initial_frequencies} and \code{final frequencies}
}

  \item{picked_ancestor}{
    Default is "ALL", where different colors indicate different ancestors. Alternatively, for clarity, the user can specify a specific ancestral allele, and only that allele is plotted
}
}
\value{
  a ggplot2 object
}

\author{
  Thijs Janzen
}

\examples{
  s <- 0.1
  select_matrix <- matrix(nrow = 1, ncol = 5)
  select_matrix[1, ] <- c(0.25, 1.0, 1 + 0.5 * s, 1 + s, 0)

  track_freq <- c(0.2, 0.3, 100)

  selected_pop <- create_population_selection(pop_size = 1000,
                                                    number_of_founders = 10,
                                                    total_runtime = 11,
                                                    morgan = 1,
                                                    select_matrix,
                                                    seed = 1234,
                                                    track_frequency = track_freq)
  plot_difference_frequencies(selected_pop,
                  picked_ancestor = "ALL")
}

