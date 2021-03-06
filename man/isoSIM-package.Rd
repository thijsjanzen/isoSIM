\name{isoSIM-package}
\alias{isoSIM-package}
\alias{isoSIM}
\docType{package}
\title{
 This package simulates the generation of an isofemale line
}
\description{
   Version 0.57 - Added function to generate admixed individuals \cr
   Version 0.56 - Added starting frequencies to 'simulate_admixture' \cr
   Version 0.55 - extended 'calculate_marker_frequency' to handle a vector of locations \cr
   Version 0.55 - increased accuracy of choosing a random position for recombination, this should prevent the rare bug fixed in version 0.54 \cr
   Version 0.54 - Fixed a MAJOR bug regarding recombination: in rare cases, a crossover position could be picked on an existing junction, due to the limited number of digits in uniform() \cr
   Version 0.54 - Improved plot_difference_frequencies to handle modified input \cr
   Version 0.53 - Added multiplicative_selection \cr
   Version 0.52 - Added plot_difference_frequencies \cr
   Version 0.51 - Added tajima's d calculation \cr
   Version 0.50 - added simulated_admixture until \cr
   Version 0.49 - Added 'simulate' to cpp \cr
   Version 0.48 - Added a general 'simulate' function \cr
   Version 0.47 - Changed the effect of migration \cr
   Version 0.46 - Added joyplot & increase_ancestor \cr
   Version 0.45 - Removed create_two_populations \cr
   Version 0.44 - Added tracking regions \cr
   Version 0.43 - Fixed bugs in select_population \cr
   Version 0.42 - Added initial and final frequency tables \cr
   Version 0.41 - Added multiple marker support \cr
   Version 0.40 - Collapsed selection functions \cr
   Version 0.39 - Added support for non-additive selection \cr
   Version 0.38 - Added track frequencies \cr
   Version 0.37 - Removed selection on regions \cr
   Version 0.36 - Added progress_bar option \cr
   Version 0.35 - Added calculate_marker_frequency \cr
   Version 0.34 - Added selection_markers \cr
   Version 0.33 - Fixed bugs in selection \cr
   Version 0.32 - Moved Fish.h code to Fish.cpp \cr
   Version 0.31 - Changed random number generator to R based \cr
   Version 0.30 - Added Recombination = 1 code \cr
   Version 0.29 - Changed internal junction representation: removed .left \cr
   Version 0.28 - Reverted to Agner Fog Random number generation \cr
   Version 0.27 - Speed up return types \cr
   Version 0.26 - Added class verification code \cr
   Version 0.25 - Squashed plotting bug \cr
   Version 0.24 - Removed Output.cpp \cr
   Version 0.23 - Removed number_of_founders from calc_allele_spectrum \cr
   Version 0.22 - Added save and load functions \cr
   Version 0.21 - Changed random-seed management \cr
   Version 0.20 - Removed superfluous code \cr
   Version 0.19 - Removed number_of_founders from Fst and LD code \cr
   Version 0.18 - Start of tracking changes \cr
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Thijs Janzen (thijsjanzen@gmail.com)
}
