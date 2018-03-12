\name{select_population_markers}
\alias{select_population_markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate a population, given selection for a specific ancestor on user specified markers
}
\description{
  A population is forward simulated, given selection for specific ancestors on genetic markers
}
\usage{
select_population_markers(source_pop, select_matrix,
                  pop_size, total_runtime, morgan, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_pop}{
  Source population with which to start.
}
  \item{select_matrix}{
  Selection matrix indicating the areas which are under selection. Each row in the matrix contains three entries:
    \code{location}{
      location of the marker under selection (in Morgan)
    }
    \code{ancestor} {
      if there is genomic material present of this ancestor at \code{location}, the individual received positive selection.
    }
    \code{selection} {
      The amount of positive selection received when there is a matching ancestor at the site.
    }
  }
  \item{pop_size}{
  Size of the population
}
  \item{total_runtime}{
  Total number of generations
}
  \item{morgan}{
  Length of the chromosome in Morgan (e.g. the number of crossovers during meiosis)
}
  \item{seed}{
  Seed of the pseudo-random number generator
}

}

\value{
 A population structure
}
\author{
 Thijs Janzen
}

\examples{
\dontrun{
  number_founders = 20
  wildpop =  create_population(pop_size = 100,
                             number_of_founders = number_founders,
                             total_runtime = 10,
                             morgan = 1,
                             seed = 666)

  selection_matrix <- matrix(ncol = 3, nrow = 1)
  selection_matrix[1, ] <- c(0.5, 0, 1)

  selected_pop <- select_population(wildpop, selection_matrix,
                                  pop_size = 100,
                                  total_runtime = 10,
                                  morgan = 1,
                                  seed = 12345)
  }
}