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
create_population_from_individuals(individuals, pop_size,
              total_runtime, morgan, seed, progress_bar)
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
\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}

}

\value{
 A population object
}

\author{
  Thijs Janzen
}

\examples{

\dontrun{

  two_populations <- create_two_populations(pop_size = 100,
                                          number_of_founders = 20,
                                          total_runtime = 2,
                                          morgan = 1,
                                          seed = 42,
                                          overlap = 0.25,
                                          progress_bar = FALSE)

  isofemale_1 <- create_iso_female(source_pop = two_populations$Population_1,
                               n = 1,
                               inbreeding_pop_size = 100,
                               run_time = 1000000,
                               morgan = 1,
                               progress_bar = FALSE)

 isofemale_2 <- create_iso_female(source_pop = two_populations$Population_2,
                               n = 1,
                               inbreeding_pop_size = 100,
                               run_time = 1000000,
                               morgan = 1,
                               progress_bar = FALSE)

  mixed_population <- create_population_from_individuals(list(isofemale_1[[1]],
                                                         isofemale_2[[1]]),
                                                       pop_size = 100, total_runtime = 10,
                                                       morgan = 1, seed = 42,
                                                       progress_bar = FALSE)
}

}