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
                  track_frequency = FALSE, progress_bar = TRUE,
                  multiplicative_selection = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_pop}{
  Source population with which to start.
}
 \item{select_matrix}{
  Selection matrix indicating the markers which are under selection. Each row in the matrix contains five entries:
    \code{location}{
      location of the marker under selection (in Morgan)
    }
    \code{fitness of wildtype (aa)}
    \code{fitness of heterozygote (aA)}
    \code{fitness of homozygote mutant (AA)}
    \code{Ancestral type that representes the mutant allele A}
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
  If TRUE, every generation the frequency of all ancestors at the marker under selection is tracked. Currently only supports one marker under selection. Alternatively, the user can provide a region for which a fixed number of markers is tracked. The user does this by providing a vector with three entries: \code{start}, \code{end}, \code{number_of_markers}.
}

\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}

\item{multiplicative_selection}{
  Default: TRUE. If TRUE, fitness is calculated for multiple markers by multiplying fitness values for each marker. If FALSE, fitness is calculated by adding fitness values for each marker.
}

}

\value{
A list with: \code{population} a population object, and three tibbles with allele frequencies: \code{frequencies} (if \code{track_frequency} = TRUE), \code{initial_frequencies} and \code{final_frequencies}. Each tibble contains four columns, \code{time}, \code{location}, \code{ancestor} and \code{frequency}, which indicates the number of generations, the location along the chromosome of the marker, the ancestral allele at that location in that generation, and finally, the frequency of that allele.
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

selection_matrix = matrix(ncol=5, nrow = 1)

selection_matrix[1,] = c(0.5, 1.0, 1.2, 1.0, 1)

selected_pop <- select_population(wildpop,
                                  selection_matrix,
                                  pop_size = 100,
                                  total_runtime = 10,
                                  morgan = 1,
                                  seed = 12345)
  }
}