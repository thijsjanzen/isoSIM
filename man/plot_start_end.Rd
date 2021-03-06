\name{plot_start_end}
\alias{plot_start_end}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot both the starting frequencies and the final frequencies in one plot
}
\description{
 This function plots the distribution of both the starting and the final frequencies in one plot
}
\usage{
plot_start_end(results, picked_ancestor = "ALL")
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
  s <- 0.01
  select_matrix <- matrix(nrow = 1, ncol = 5)
  select_matrix[1, ] <- c(0.25, 1.0, 1 + 0.5 * s, 1 + s, 0)

  track_freq <- c(0.2, 0.3, 100)

  selected_pop <- isoSIM::create_population_selection(pop_size = 1000,
                                                    number_of_founders = 10,
                                                    total_runtime = 11,
                                                    morgan = 1,
                                                    select_matrix,
                                                    seed = 1234,
                                                    track_frequency = track_freq)
  plot_start_end(selected_pop,
                  picked_ancestor = "ALL")
}

