\name{ch7lie}
\docType{data}
\alias{ch7lie}
\title{LIE Data for the Vectorized Root System of a Young Oak Tree}
\description{
  This data set obtained with Data Analysis of Root Tracings summarizes the attributes of each point used to construct the links constituting the interpreted root system of a young oak tree (\emph{Quercus robur} L.) growing in a vertical rhizobox.    
}
\usage{ch7lie}
\format{A data frame containing the following columns:
\itemize{
\item \bold{Num}: the identification number of each point (i.e. the junctions between two consecutive links) used to construct the vectorized root system.
\item \bold{Date}: the identification number coding for the date at which a point was seen for the first time.
\item \bold{Bran}: the branching state (\sQuote{true} denotes the start of a branching link while \sQuote{false} denotes the start of a continuing link).
\item \bold{Apic}: the apical state (\sQuote{true} denotes a point ending a root while \sQuote{false} is related to a point starting a new link).
\item \bold{Prec}: the identification number of the previous point used to construct a link.
\item \bold{Suiv}: the identification number of the following point used to contruct a link.
\item \bold{X}: the X coordinate of each point used to construct the vectorized root system.
\item \bold{Y}: the Y coordinate of each point used to construct the vectorized root system.}
}
\source{\url{www6.paca.inra.fr/psh/Modeles-et-logiciels/Dart}
}
\references{
Pages L., Serra V. (1994) Growth and branching of the taproot of young oak trees - a dynamic study, \emph{Journal of Experimental Botany}, \bold{45}, 1327--1334.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.
}
\keyword{datasets}