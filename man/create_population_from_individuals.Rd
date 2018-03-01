\name{create_population_from_individuals}
\alias{create_population_from_individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create a population from isofemale individuals
}
\description{
  Create a population from isofemale individuals
}
\usage{
create_population_from_individuals(indiv_1, indiv_2, pop_size,
              total_runtime, morgan, seed, write_to_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individuals}{
  List of isofemale individuals to be used to create a population from.
}
  \item{pop_size}{
  Number of individuals in the population
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
 A population object
}

\author{
  Thijs Janzen
}

\examples{
  two_populations <- create_two_populations(pop_size = 100,
                                          number_of_founders = 20,
                                          total_runtime = 5,
                                          morgan = 1,
                                          seed = 42,
                                          overlap = 0.25,
                                          write_to_file = FALSE)

  isofemale_1 <- create_iso_female(source_pop = two_populations$Population_1,
                               n = 1,
                               inbreeding_pop_size = 100,
                               run_time = 1000000,
                               morgan = 1)

 isofemale_2 <- create_iso_female(source_pop = two_populations$Population_2,
                               n = 1,
                               inbreeding_pop_size = 100,
                               run_time = 1000000,
                               morgan = 1)

  mixed_population <- create_population_from_individuals(c(isofemale_1[[1]],
                                                         isofemale_2[[1]]),
                                                       pop_size = 100, total_runtime = 100,
                                                       morgan = 1, seed = 42, write_to_file = FALSE)


}