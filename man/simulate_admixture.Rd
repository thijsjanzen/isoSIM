\name{simulate_admixture}
\alias{simulate_admixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Individual based simulation of the breakdown of contiguous ancestry blocks.
}
\description{
  Individual based simulation of the breakdown of contiguous ancestry blocks, with or without selection. Simulations can be started from scratch, or from a predefined input population.
}
\usage{
simulate_admixture(input_population = NA, pop_size,
                     number_of_founders,
                     total_runtime,
                     morgan,
                     seed,
                     select_matrix = NA,
                     progress_bar = TRUE,
                     track_junctions = FALSE,
                     track_frequency = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_population}{
    Potential earlier simulated population used as starting point for the simulation. If not provided by the user, the simulation starts from scratch.
}
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

\item{seed}{
  Seed of the pseudo-random number generator
}

 \item{select_matrix}{
  Selection matrix indicating the markers which are under selection. If not provided by the user, the simulation proceeds neutrally. If provided, each row in the matrix should contain five entries:
    \code{location}{
      location of the marker under selection (in Morgan)
    }
    \code{fitness of wildtype (aa)}
    \code{fitness of heterozygote (aA)}
    \code{fitness of homozygote mutant (AA)}
    \code{Ancestral type that representes the mutant allele A}
    }

\item{progress_bar}{
   Displays a progress_bar if TRUE. Default value is TRUE
}

\item{track_junctions}{
   Track the average number of junctions over time if TRUE
}

\item{track_frequency}{
  If TRUE, every generation the frequency of all ancestors at the marker under selection is tracked. Currently only supports one marker under selection. \code{track_frequency}{
  If TRUE, every generation the frequency of all ancestors at the marker under selection is tracked. Currently only supports one marker under selection. Alternatively, the user can provide a region for which a fixed number of markers is tracked. The user does this by providing a vector with three entries: \code{start}, \code{end}, \code{number_of_markers}.
}
}

}

\value{
A list with: \code{population} a population object, and three tibbles with allele frequencies (if track_frequency = TRUE): \code{frequencies} (if \code{track_frequency} = TRUE), \code{initial_frequencies} and \code{final_frequencies}. Each tibble contains four columns, \code{time}, \code{location}, \code{ancestor} and \code{frequency}, which indicates the number of generations, the location along the chromosome of the marker, the ancestral allele at that location in that generation, and finally, the frequency of that allele.
}

\author{
  Thijs Janzen
}

\examples{

}