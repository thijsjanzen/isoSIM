\name{create_population_selection_markers}
\alias{create_population_selection_markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create a population taking into account selection on specific markers
}
\description{
  Create a population while specific markers (e.g. SNPs) are under selection to be of a specific ancestor.
}
\usage{
create_population_selection_markers(pop_size, number_of_founders,
            total_runtime, morgan, select_matrix, seed)
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
  selection_matrix <- matrix(ncol = 3, nrow = 1)
  selection_matrix[1, ] <- c(0.5, 0, 1)

 selected_pop <- create_population_selection_markers(selection_matrix,
                                                      pop_size = 100,
                                                      number_of_founders = 10,
                                                      total_runtime = 10,
                                                      morgan = 1,
                                                      seed = 1234)
}