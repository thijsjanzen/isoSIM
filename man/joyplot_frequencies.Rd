\name{joyplot_frequencies}
\alias{joyplot_frequencies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 make a joy plot of the distribution of allele frequencies within a region
}
\description{
 This function plots the distribution of allele frequencies within a region over time, making use of a 'joyplot'
}
\usage{
joyplot_frequencies(frequencies, time_points, picked_ancestor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frequencies}{
    A tibble containing four columns: \code{time}, \code{location}, \code{ancestor}, \code{frequency}. Typically one of the items returned by \code{create_population_selection} or \code{select_population} when the user specifies \code{track_frequency}.
}
  \item{time_points}{
    A sequence of time points for which the user wants to create the joyplot
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

  markers <- seq(from = 0.2, to = 0.3, length.out = 100)

  selected_pop <- create_population_selection(pop_size = 1000,
                                              number_of_founders = 10,
                                              total_runtime = 11,
                                              morgan = 1,
                                              select_matrix,
                                              seed = 1234)
  joyplot_frequencies(selected_pop$frequencies,
                    0:11, picked_ancestor = "ALL")
}

