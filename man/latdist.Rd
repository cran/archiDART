\name{latdist}
\alias{latdist}
\title{Computing Lateral Root Distributions on the Mother Roots}
\description{Calculation of local lateral root lengths and densities, as well as distances between neighbouring lateral roots according to the distance between the branching point and the parent root base (DBase) using Data Analysis of Root Tracings (DART) output files.  
}
\usage{
latdist(inputrac, output=c("lrd","dtp"), res=NULL, unitlength="px",
        int.length=NULL, interpol=NULL)
}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART.}
  \item{output}{A character string specifying the type of output that must be returned by the function. Two values are acceptable for this argument: \dQuote{lrd} and \dQuote{dtp}. If \code{output="lrd"} (default value), the function will compute lateral root lengths and densities according to DBase for each mother root of each vectorized root system. If \code{output="dtp"}, the function will compute the distances between neighbouring lateral roots according to DBase for each mother root of each analysed root system (see details).}
  \item{res}{A numeric value specifying the image resolution (dots per inch, dpi) of the pictures used in DART for the vectorization of the root systems. This argument must only be specified if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that will be used by the function to express the results in the final list. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimeters or \dQuote{mm} for millimeters.}
  \item{int.length}{A numeric value specifying the interval length (expressed in \code{unitlength}) that will be used by the function to calculate local lateral root lengths and densities. This value must only be specified when \code{output="lrd"} (see details).}
  \item{interpol}{A numeric value specifying the number of points that must be used by the function on each mother root to calculate local lateral root lengths and densities. By default, \code{interpol=NULL} means that the function will calculate local lateral root lengths and densities at each branching point of each mother root. A numeric value could be specified when \code{output="lrd"} (see details).}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex root system architectures from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  After vectorization of a root system, a final data set can be exported as a table file (rac file) containing several attributes for each individual root, notably the root length at each observation date and the distance between each branching point to the parent root base (DBase).  This R function allows the calculation of either local lateral root lengths and densities or distances between neighbouring lateral roots according to DBase on each mother root consisting a vectorized root system.  

The function \code{latdist} will successively locate and read the rac files created by DART and calculate RSA parameters.  First, \code{latdist} will check if all arguments have been written in the suitable mode.  If \code{res}, \code{interpol} and \code{int.length} are specified, the function will check if the numeric values required by these arguments are positive.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputrac} to locate and read the raw rac files created by DART.  To run \code{latdist} efficiently, the DART files must have been saved with their appropriate extensions (.rac).  Third, depending on the \code{output} value, the function will compute RSA parameters.  If \code{output="lrd"}, the function will use \code{interpol} to select the appropriate DBase values on which the calculation of RSA parameters should be performed.  Then, a lateral root density and a total lateral root length will be calculated on intervals possessing a length equal to \code{int.length} and centred on each DBase value.  If \code{output="dtp"}, only the distance between neighbouring lateral roots on each mother root will be computed by the function.  
}
\value{
Returns a list including:

\item{root}{A list of data frames. Each element of the list is named as its corresponding rac file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system in the corresponding rac file), \bold{Ord} (the branching order), \bold{LatRootNum} (the number of lateral roots), \bold{FinalRootLength} (the root length at the last observation date expressed in \code{unitlength}), \bold{LatRootDensity} (the lateral root density calculated as LatRootNum/FinalRootLength and expressed in root(s)/\code{unitlength}).}

\item{res}{A list of secondary lists. Each element of the primary list (res) is named as its corresponding rac file. The secondary lists contain as much elements as roots consisting a vectorized root system. Each element of the secondary lists could be either a \code{NULL} value if a root do not have lateral roots or if no DBase value is comprised in [\code{int.length}/2, FinalRootLength-(\code{int.length}/2)], or a data frame made of two or three columns depending on the \code{output} value. 

\bold{If \code{output="lrd"}}: \bold{DBase} (the distance between the branching point and the parent root base expressed in \code{unitlength}), \bold{LRD} (Local lateral root density expressed in root(s)/\code{unitlength}), \bold{LRL} (Local lateral root length expressed in \code{unitlength}/\code{unitlength}).

\bold{If \code{output="dtp"}}: \bold{DBase} (the distance between the branching point and the parent root base expressed in \code{unitlength}), \bold{DTP} (the distance to the previous lateral root expressed in \code{unitlength}).}
}
\references{
Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.
}
\author{Benjamin M. Delory, Loic Pages
}
\examples{
## Importation of an example data set
data(ch7rac)

## Creation of a folder (DART_latrootdensity) inside a temporary directory to contain the rac
## example data sets created by DART
dirDART <- paste(tempdir(), "/DART_latrootdensity", sep="")
dir.create(dirDART)

## Definition of the number of rac files that will be saved inside dirDART
n <- 1

## Creation of rac files inside dirDART
for (i in 1:n) {write.table(ch7rac, file=paste(dirDART,"/ch7_",i,".rac", sep=""), row.names=FALSE)}

## Use of latdist to compute RSA parameters from DART output files
## Plotting the results for the primary root
ch7_output_4 <- latdist(inputrac=dirDART, output="lrd", res=75, unitlength="cm",
int.length=5)
plot(ch7_output_4$res$ch7_1[[1]]$DBase, ch7_output_4$res$ch7_1[[1]]$LRD, pch=16, 
xlab="DBase (cm)",ylab="Lateral root density (roots/cm)", main="LRD-interpol=NULL-int.length=5",
las=1, bty="l", xaxp=c(0,90,9))
plot(ch7_output_4$res$ch7_1[[1]]$DBase, ch7_output_4$res$ch7_1[[1]]$LRL, pch=16, 
xlab="DBase (cm)",ylab="Lateral root length (cm/cm)", main="LRL-interpol=NULL-int.length=5",
las=1, bty="l", xaxp=c(0,90,9))

ch7_output_5 <- latdist(inputrac=dirDART, output="lrd", res=75, unitlength="cm",
int.length=5, interpol=1000)
plot(ch7_output_5$res$ch7_1[[1]]$DBase, ch7_output_5$res$ch7_1[[1]]$LRD, pch=16, 
xlab="DBase (cm)", ylab="Lateral root density (roots/cm)", main="LRD-interpol=1000-int.length=5", 
las=1, bty="l", xaxp=c(0,90,9))
plot(ch7_output_5$res$ch7_1[[1]]$DBase, ch7_output_5$res$ch7_1[[1]]$LRL, pch=16, 
xlab="DBase (cm)", ylab="Lateral root length (cm/cm)", main="LRL-interpol=1000-int.length=5", 
las=1, bty="l", xaxp=c(0,90,9))

ch7_output_6 <- latdist(inputrac=dirDART, output="dtp", res=75, unitlength="cm")
plot(ch7_output_6$res$ch7_1[[1]]$DBase, ch7_output_6$res$ch7_1[[1]]$DTP, pch=16, 
xlab="DBase (cm)", ylab="Distance to the previous root (cm)", main="DTP", las=1, bty="l",
xaxp=c(0,90,9))

unlink(dirDART, recursive=TRUE)
}