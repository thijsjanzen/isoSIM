\name{create_two_populations}
\alias{create_two_populations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Create two populations, possibly sharing ancestry
}
\description{
  This function creates two independent populations with a user provided degree of shared ancestry.
}
\usage{
create_two_populations(pop_size, number_of_founders,
                      total_runtime, morgan, seed, overlap)
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
  \item{overlap}{
    Degree of shared ancestry between the populations, reflected in the probability of sharing an ancestor.
}
\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}

}

\value{
 An object containing two population objects: \code{Population_1} and \code{Population_2}.
}

\author{
 Thijs Janzen
}

\examples{
two_populations <- create_two_populations(pop_size = 100,
                                          number_of_founders = 20,
                                          total_runtime = 10,
                                          morgan = 1,
                                          seed = 42,
                                          overlap = 0.25)

}
