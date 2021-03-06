\name{calculate_fst}
\alias{calculate_fst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate FST
}
\description{
 The FST value between two populations is calculated, given a number of markers. Markers are superimposed upon the (known) ancestry along the chromosome for all sampled individuals. Markers can be chosen to be regularly spaced, or randomly distributed.
}
\usage{
calculate_fst(pop1, pop2,
        sampled_individuals, number_of_markers = 100, random_markers = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop1}{
  Population object
}
  \item{pop2}{
  Population object
}

  \item{sampled_individuals}{
  Number of individuals to base the FST upon. Individuals are randomly drawn from each population, a lower number speeds up calculations.
}
  \item{number_of_markers}{
  Number of markers along the chromosome used to calculate FST metrics.
}
  \item{random_markers}{
  If TRUE, markers are randomly spaced along the chromosome, if FALSE, markers are equidistantly spaced along the chromosome.
}
}
\details{
  Uses the function \code{wc} from the package \code{hierfstat} to calculate the FST. The function \code{wc} computes the Weir and Cockerham F statistic.
}
\value{
FST value
}

\author{
Thijs Janzen
}


\examples{
\dontrun{
two_populations <- create_two_populations(pop_size = 1000,
                                          number_of_founders = 20,
                                          total_runtime = 100,
                                          morgan = 1,
                                          seed = 42,
                                          overlap = 0.25,
                                          write_to_file = FALSE)

  FST <- calculate_fst(two_populations$Population_1,
                 two_populations$Population_2,
                 sampled_individuals = 10,
                 number_of_markers = 100,
                 random_markers = TRUE)
  }
}