\name{create_population}
\alias{create_population}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a population
}
\description{
  Creates a population of individuals
}
\usage{
create_population(pop_size, number_of_founders,
                       total_runtime, morgan, seed, write_to_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop_size}{
  Total number of diploid individuals
}
  \item{number_of_founders}{
  Number of unique ancestors
}
  \item{total_runtime}{
  Number of generations allowed for inbreeding
}
  \item{morgan}{
  Length of the chromosome in Morgan (e.g. the number of crossovers during meiosis)
}
  \item{seed}{
  Seed of the pseudo-random number generator
}
  \item{write_to_file}{
  Boolean indicating whether the population should be stored into a text file as well.
}
}
\value{
  A population structure
}
\author{
 Thijs Janzen
}

\examples{
wildpop =  create_population(pop_size = 100,
                             number_of_founders = 10,
                             total_runtime = 1000,
                             morgan = 1,
                             seed = 42,
                             write_to_file = FALSE)
}