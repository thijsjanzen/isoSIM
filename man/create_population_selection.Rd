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
            total_runtime, morgan, select_matrix, selection, seed, write_to_file)
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
  Selection matrix indicating the areas which are under selection. Each row in the matrix contains three entries:
    \code{start}{
      start of the genomic section under selection (in Morgan)
    }
    \code{end}{
      end of the genomic section under selection (in Morgan)
    }
    \code{ancestor} {
      if there is genomic material present of this ancestor between \code{start} and \code{end}, the individual received positive selection.
    }
}
  \item{selection}{
   Amount of selection. Selection is calculated as (1 + selection) * (end - start), hence, selection is dependent on the length of the genomic stretch as well.
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
selection_matrix = matrix(ncol=3, nrow = 2)
selection_matrix[1,] = c(0.2, 0.3, 1)
selection_matrix[2,] = c(0.6, 0.7, 2)

selection_coefficient <- 0.2

selected_pop <- create_population_selection(pop_size = 100,
                                            number_of_founders = 10,
                                            total_runtime = 10,
                                            morgan = 1,
                                            selection_matrix,
                                            selection = selection_coefficient,
                                            seed = 12345,
                                            write_to_file = FALSE)
}