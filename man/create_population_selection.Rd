\name{create_population_selection}
\alias{create_population_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create a population taking into account selection
}
\description{
  Create a population while specific genomic areas are under selection to be of a specific ancestor.
}
\usage{
create_population_selection(pop_size, number_of_founders,
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
  If TRUE, every generation the frequency of all ancestors at the marker under selection is tracked. Currently only supports one marker under selection. \item{track_frequency}{
  If TRUE, every generation the frequency of all ancestors at the marker under selection is tracked. Currently only supports one marker under selection. Alternatively, the user can provide a region for which a fixed number of markers is tracked. The user does this by providing a vector with three entries: \code{start}, \code{end}, \code{number_of_markers}.
}
}

\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}

}

\value{
A list with: \code{population} a population object, and three tibbles with allele frequencies: \code{frequencies} (if \code{track_frequency} = TRUE), \code{initial_frequencies} and \code{final_frequencies}. Each tibble contains four columns, \code{time}, \code{location}, \code{ancestor} and \code{frequency}, which indicates the number of generations, the location along the chromosome of the marker, the ancestral allele at that location in that generation, and finally, the frequency of that allele.
}

\author{
  Thijs Janzen
}

\examples{
selection_matrix <- matrix(ncol=5, nrow = 1)
selection_matrix[1,] <- c(0.5, 1.0, 1.0, 1.2, 1)

selected_pop <- create_population_selection(pop_size = 100,
                                            number_of_founders = 10,
                                            total_runtime = 10,
                                            morgan = 1,
                                            selection_matrix,
                                            seed = 12345,
                                            progress_bar = TRUE)
}