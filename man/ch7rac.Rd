\name{ch7rac}
\docType{data}
\alias{ch7rac}
\title{RAC Data for the Vectorized Root System of a Young Oak Tree}
\description{For each observation date, this dataset obtained with Data Analysis of Root Tracings summarizes the length of each root constituting the vectorized root system of a young oak tree (\emph{Quercus robur} L.) growing in a vertical rhizobox.    
}
\usage{ch7rac}
\format{A data frame containing the following columns:
\itemize{
\item \bold{Root}: the identification number of each root constituting the vectorized root system.
\item \bold{Mother}: the parent root identification number (a first-order root having no parent root, its value in this column is set to -1).
\item \bold{Ord}: the branching order (equal to 1 for the first-order root, 2 for a secondary root (or a first-order lateral root), 3 for a tertiary root (or a second-order lateral root), etc).
\item \bold{DBase}: the distance between the branching point to the parent root base.
\item \bold{DApp}: the date of emergence.
\item \bold{Lengths}: the root length at each observation date (in pixels). The number of column "Lengths" is equal to the number of observation dates.}
}
\source{\url{www6.paca.inra.fr/psh/Modeles-et-logiciels/Dart}
}
\references{
Pages L., Serra V. (1994) Growth and branching of the taproot of young oak trees - a dynamic study, \emph{Journal of Experimental Botany}, \bold{45}, 1327--1334.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.
}
\keyword{datasets}