\name{create_population_selection_twoalleles}
\alias{create_population_selection_twoalleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create a population taking into account selection
}
\description{
  Create a population while specific genomic areas are under selection to be of a specific ancestor.
}
\usage{
create_population_selection_twoalleles(pop_size, number_of_founders,
            total_runtime, morgan, select_matrix,
            seed, track_frequency = FALSE, progress_bar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop_size}{
  Number of individuals in the population
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
s <- 0.1
selection_matrix <- matrix(ncol=5, nrow = 1)
selection_matrix[1,] <- c(0.5, 1.0, 1.0 + 0.5 * s, 1.0 + s, 0)

selected_pop <- create_population_selection(pop_size = 100,
                                            number_of_founders = 10,
                                            total_runtime = 10,
                                            morgan = 1,
                                            selection_matrix,
                                            seed = 12345,
                                            progress_bar = TRUE)
}