\name{calculate_marker_frequency}
\alias{calculate_marker_frequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate allele frequencies at a specific marker location
}
\description{
  Calculate the relative frequency of each ancestor in the population at a specific marker location
}
\usage{
calculate_marker_frequency(pop, location)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
  Population for which to estimate allele frequencies at the given marker
}
  \item{location}{
  A vector or scalar of location(s) along the chromosome for which allele frequencies are to be calculated. Locations are in Morgan.
}
}

\value{
  A tibble containing the frequency of each present ancestor at the provided location. Ancestors with frequency = 0 are dropped out of the table. The tibble contains three columns: location, ancestor and frequency.
}

\author{
  Thijs Janzen
}

\examples{
\dontrun{
  number_founders = 20
  wildpop =  create_population(pop_size = 1000,
                             number_of_founders = number_founders,
                             total_runtime = 10,
                             morgan = 1,
                             seed = 666)

  calculate_marker_frequency(wildpop, 0.5)
  frequencies <- calculate_marker_frequency(wildpop, seq(0.4, 0.5, by = 0.001))
  ggplot(frequencies, aes(x = location, y = frequency, col = ancestor)) +
      geom_step()

}
}