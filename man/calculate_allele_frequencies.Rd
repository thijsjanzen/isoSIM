\name{calculate_allele_frequencies}
\alias{calculate_allele_frequencies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate allele frequencies
}
\description{
  Calculate for a number of regularly spaced markers the relative frequency of each ancestor in the population.
}
\usage{
calculate_allele_frequencies(source_pop, number_of_founders, step_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_pop}{
  Population for which to estimate allele frequencies
}
  \item{number_of_founders}{
  Number of founders in the population
}
  \item{step_size}{
  Distance between markers, a smaller step_size provides a more complete overview, but is also more computationally heavy.
}
}
\details{
  Markers are equidistantly spaced, with a distance of \code{step_size} in between them.
}
\value{
  A tibble containing the allele frequencies
}

\author{
  Thijs Janzen
}

\examples{
  number_founders = 20
  wildpop =  create_population(pop_size = 1000,
                             number_of_founders = number_founders,
                             total_runtime = 10,
                             morgan = 1,
                             seed = 666,
                             write_to_file = FALSE)

 freq_output <- calculate_allele_frequencies(wildpop,
                                            number_of_founders = number_founders,
                                            step_size = 0.001)

  require(ggplot2)
  ggplot(freq_output, aes(x=location, y = frequency, col = as.factor(ancestor))) +
    geom_line()
}