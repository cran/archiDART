\name{architect}
\alias{architect}
\title{Computing Root System Architecture Parameters}
\description{Calculation of common root system architecture (RSA) parameters from Data Analysis of Root Tracings (DART) output files.      
}
\usage{
architect(inputrac, inputtps, res=NULL, unitlength="px", rootdiv=1)}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as inputtps.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as inputrac.}
  \item{res}{A numeric value specifying the image resolution (dots per inch, dpi) of the pictures used in DART for the vectorization of the root systems. This argument must only be specified if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that will be used by the function to express the results in the final data frame. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimeters or \dQuote{mm} for millimeters.}
  \item{rootdiv}{Either a numeric value or a vector of numeric values defining a class number or the class limits that will be used by the function to divide the primary root in zones. Using this argument, \code{architect} will also compute the distribution of the second-order lateral root number, length and density on the primary root.}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex RSAs from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  After vectorization of a root system, a final data set can be exported as a table file (rac file) containing several attributes for each individual root, notably the root length at each observation date.  A second table file (tps file) containing the root system age at each observation date can also be exported from DART.  These two data sets can finally be used to calculate derived root system architecture parameters.  In this context, \code{architect} is an R function allowing the calculation of common RSA parameters at each observation date for each DART rac file: the total root system length, the length of the primary root, the total number of lateral roots, the total lateral root length, the length and the number of lateral roots for each branching order, the mean length of lateral roots for each branching order, the density of second-order lateral roots on the primary root, the primary and the lateral root growth rates, as well as the growth rate of the entire root system.  When needed, the distribution of the second-order lateral root number, length and density on the primary root can also be calculated (via \code{rootdiv}).  

The function \code{architect} will successively locate and read the rac and tps files created by DART and calculate RSA parameters.  First, \code{architect} will check if all arguments have been written in the suitable mode.  If \code{res} is specified, the function will check if the numeric value required by this argument is positive.  The function will also automatically check if the numeric elements in \code{rootdiv} are superior or equal to zero and sorted by increasing values.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputrac} and \code{inputtps} to locate and read the raw rac and tps files created by DART respectively.  To run \code{architect} efficiently, DART files must have been saved with their appropriate extensions (.rac or .tps).  Third, the function will create vectors and matrices for RSA parameters calculation.  If \code{inputtps} contains a single tps file, it will be used by \code{architect} for each rac file located in \code{inputrac}.  If \code{inputps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of rac files in \code{inputrac} and corresponding rac and tps files must have the same name.  If an error occurs at this step, a message will be returned by the function.  

Growth rates are calculated as follows:
\itemize{
\item \emph{Primary root growth rate}: for the first observation date, it is calculated as the ratio of the primary root length to the root system age. For other observation dates (t), it is calculated as the difference between the primary root length at time t and t-1 divided by the difference between the root system age at time t and t-1.
\item \emph{Second- and third-order lateral root growth rates}: for the first observation date, it is calculated as the ratio of the second- or third-order lateral root length to the root system age. For other observation dates (t), it is calculated as the difference between the second- or third-order lateral root length at time t and t-1 divided by the difference between the root system age at time t and t-1.
\item \emph{Growth rate of the root system}: for the first observation date, it is calculated as the ratio of the total root system length to the root system age. For other observation dates (t), it is calculated as the difference between the total root system length at time t and t-1 divided by the difference between the root system age at time t and t-1.}
}
\value{
Returns a data frame summarizing all RSA parameters calculated for each rac file located in \code{inputrac}. 

The data frame contains the following columns: \bold{FileName} (the name of the rac files), \bold{Time} (the observation date expressed in the same unit as in the tps file(s) (=unittime)), \bold{TRL} (the total length of the root system expressed in \code{unitlength}), \bold{GRTR} (the growth rate of the root system expressed in \code{unitlength}/unittime), \bold{PRL} (the primary root length expressed in \code{unitlength}), \bold{GRPR} (the growth rate of the primary root expressed in \code{unitlength}/unittime), \bold{TNLR} (the total number of lateral roots), \bold{TLRL} (the total length of lateral roots expressed in \code{unitlength}), \bold{NxLR} (the total number of lateral roots with a branching order equal to x), \bold{LxLR} (the total length of lateral roots with a branching order equal to x expressed in \code{unitlength}), \bold{MLxLR} (the mean length of lateral roots with a branching order equal to x expressed in \code{unitlength}), \bold{GRxL} (the growth rate of lateral roots with a branching order equal to x expressed in \code{unitlength}/unittime), \bold{D2LR} (the density of second-order lateral roots on the primary root expressed in root(s)/\code{unitlength}).

\bold{If \code{rootdiv} is a numeric value}, the following columns will be added to the data frame: \bold{N2LR.LayerX} (the total number of second-order lateral roots in layer X), \bold{L2LR.LayerX} (the total length of second-order lateral roots in layer X expressed in \code{unitlength}), \bold{D2LR.LayerX} (the density of second-order lateral roots on the primary root in the layer X expressed in root(s)/\code{unitlength}).

\bold{If \code{rootdiv} is a vector of numeric values}, the following columns will be added to the data frame: \bold{N2LR.XtoY} (the total number of second-order lateral roots in a layer delimited by X and Y on the primary root), \bold{L2LR.XtoY} (the total length of second-order lateral roots in a layer delimited by X and Y on the primary root expressed in \code{unitlength}), \bold{D2LR.XtoY} (the density of second-order lateral roots in a layer delimited by X and Y on the primary root expressed in root(s)/\code{unitlength}).
}
\references{
Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\examples{
## Importation of example data sets
data(ch7rac)
data(ch7tps)

## Creation of a folder (DART_architect) inside a temporary directory to contain the rac and tps
## example data sets created by DART
dirDART <- paste(tempdir(), "/DART_architect", sep="")
dir.create(dirDART)

## Definition of the number of rac and tps files that will be saved inside dirDART
n <- 2

## Creation of rac and tps files inside dirDART
for (i in 1:n) {
write.table(ch7rac, file=paste(dirDART,"/ch7_",i,".rac", sep=""), row.names=FALSE)
write.table(ch7tps, file=paste(dirDART,"/ch7_",i,".tps", sep=""), row.names=FALSE)}

## Use of architect to compute RSA parameters from DART output files
## rootdiv=1
ch7_output_2 <- architect(inputrac=dirDART, inputtps=dirDART, res=75, unitlength="cm")

## Use of the rootdiv argument to compute the distribution of the second-order lateral root number,
## length and density in the [0-15[, [15-30[, [30-45[ and [45-60[ first centimeters of the primary
## root
ch7_output_3 <- architect(inputrac=dirDART, inputtps=dirDART, res=75, unitlength="cm",
rootdiv=c(0,15,30,45,60))

unlink(dirDART, recursive=TRUE)
}