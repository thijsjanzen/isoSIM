\name{select_population}
\alias{select_population}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate a population, given selection for a specific ancestor
}
\description{
  A population is forward simulation, given selection on genomic sections for specific ancestors.
}
\usage{
select_population(source_pop, select_matrix,
                  pop_size, total_runtime, morgan, seed,
                  track_frequency = FALSE, progress_bar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_pop}{
  Source population with which to start.
}
  \item{select_matrix}{
  Selection matrix indicating the markers which are under selection. Each row in the matrix contains three entries:
    \code{location}{
      location of the marker under selection (in Morgan)
    }
    \code{ancestor} {
      if there is genomic material present of this ancestor at \code{location}, the individual received positive selection.
    \code{selection} {
      The amount of positive selection received when there is a matching ancestor at the site.
    }
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
\item{track_frequency}{
  If TRUE, every generation the frequency of all ancestors at the marker under selection is tracked. Currently only supports one marker under selection. See also \code{details}
}

\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}

}

\value{
A list with: \code{population} a population object, and \code{frequencies} a tibble containing the ancestral frequencies - this tibble is only included if \code{track_frequency} = TRUE.
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

selection_matrix = matrix(ncol=3, nrow = 2)

selection_matrix[1,] = c(0.2, 0.3, 1)
selection_matrix[2,] = c(0.6, 0.7, 2)

selected_pop <- select_population(wildpop, selection_matrix,
                                  pop_size = 100,
                                  total_runtime = 10,
                                  morgan = 1,
                                  seed = 12345)
  }
}