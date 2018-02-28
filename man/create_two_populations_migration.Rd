\name{create_two_populations_migration}
\alias{create_two_populations_migration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create two populations, linked by migration
}
\description{
  This function creates two independent populations linked through migration.
}
\usage{
create_two_populations_migration(pop_size, number_of_founders,
            total_runtime, morgan, seed, migration, write_to_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop_size}{
    Size of the population, the two populations are assumed to be equal in size.
}
  \item{number_of_founders}{
    The number of unique ancestors per population.
}
  \item{total_runtime}{
    Each population is allowed to inbreed until the total_runtime is reached.
}
  \item{morgan}{
    Size of the chromosome in Morgan.
}
  \item{seed}{
    The random seed.
}
  \item{migration}{
    Probability of having an individual in the next generation being the result of a migration event from the other population, rather than the result of a reproduction event in the focal population.
}
 \item{write_to_file}{
    Boolean flag indicating whether output should be written to file.
}
}

\value{
  An object containing two population objects: \code{Population_1} and \code{Population_2}.
}

\author{
 Thijs Janzen
}

\examples{
 two_populations <- create_two_populations_migration(pop_size = 1000,
                                          number_of_founders = 20,
                                          total_runtime = 100,
                                          morgan = 1,
                                          seed = 42,
                                          migration = 0.25,
                                          write_to_file = FALSE)
}

