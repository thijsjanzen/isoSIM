\name{select_population_twoalleles}
\alias{select_population_twoalleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate a population, given selection for a specific ancestor
}
\description{
  Forward simulation of a given population, where there is selection on a specific marker. Wildtype alleles are indicated by 'a', mutant alleles with a selective benefit with 'A'.
}
\usage{
select_population_twoalleles(source_pop, select_matrix,
                  pop_size, total_runtime, morgan, seed,
                  track_frequency = FALSE, progress_bar = TRUE)
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
\details{
 Fitness is applied here analogous to a single locus, multiple allele system where here, alleles are represented by the unique ancestors in the source population. Allele 'a' represents any ancestor not under selection at the given marker, whilst allele 'A' represents the ancestor providing a selective benefit.
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

selection_matrix = matrix(ncol=5, nrow = 1)
s <- 0.1

selection_matrix[1,] = c(0.5, 1.0, 1 + 0.5 * s, 1 + s, 0)

selected_pop <- select_population(wildpop, selection_matrix,
                                  pop_size = 100,
                                  total_runtime = 10,
                                  morgan = 1,
                                  seed = 12345)
  }
}