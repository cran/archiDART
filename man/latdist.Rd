\name{latdist}
\alias{latdist}
\title{Computing Lateral Root Length and Density Distribution}
\description{Calculation of local lateral root lengths and densities, as well as distances between neighbouring lateral roots according to the distance between the branching point and the parent root base (DBase) using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). More information can be found in \cite{Delory et al (2015)}.  
}
\usage{
latdist(inputrac=NULL, inputrsml=NULL, output=c("lrd","dtp"), res=NULL, unitlength="px",
        int.length=NULL, interpol=NULL, rsml.date=NULL, rsml.connect=FALSE)
}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as \code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the rsml files. This argument could be the same as \code{inputrac}.}
  \item{output}{A character string specifying the type of output that must be returned by the function. Two values are acceptable for this argument: \dQuote{lrd} and \dQuote{dtp}. If \code{output="lrd"} (default value), the function will compute lateral root lengths and densities according to DBase for each mother root of each vectorized root system. If \code{output="dtp"}, the function will compute the distances between neighbouring lateral roots according to DBase for each mother root of each analysed root system (see details).}
  \item{res}{\bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. This argument must be specified only if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that must be used by the function to express the results in the final list. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{int.length}{A numeric value specifying the interval length (expressed in \code{unitlength}) that must be used by the function to calculate local lateral root lengths and densities. This value must be specified only when \code{output="lrd"} (see details).}
  \item{interpol}{A numeric value specifying the number of points that must be used by the function on each mother root to calculate local lateral root lengths and densities. By default, \code{interpol=NULL} means that the function will calculate local lateral root lengths and densities at each branching point of each mother root. A numeric value can be specified when \code{output="lrd"} (see details).}
  \item{rsml.date}{A numeric value that must be specified when \code{is.null(inputrsml)=FALSE}. This argument specifies the root system age (the unit of time is defined by the user) and will be used only for RSA data encoded with the Root System Markup Language.}
  \item{rsml.connect}{A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. This argument will be used only for RSA data encoded with the Root System Markup Language. If \code{TRUE}, the lateral roots will be connected to their parent root if no connection information can be found in the rsml file(s) (see details). If \code{FALSE}, no connection will be made between daughter and mother roots if no connection information can be found in the rsml file(s).}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex root system architectures from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  After vectorization of a root system, a final dataset can be exported as a table file (rac file) containing several attributes for each individual root, notably the root length at each observation date and the distance between each branching point to the parent root base (DBase).  This R function allows the calculation of either local lateral root lengths and densities or distances between neighbouring lateral roots according to DBase on each mother root constituting a vectorized root system.  

The \code{latdist} function will successively locate and read the rac files created by DART and calculate RSA parameters.  First, \code{latdist} will check if all arguments have been written in the suitable mode.  If \code{res}, \code{interpol} and \code{int.length} are specified, the function will check if the numeric values required by these arguments are positive.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputrac} to locate and read the raw rac files created by DART.  To run \code{latdist} efficiently, the DART files must have been saved with their appropriate extensions (.rac).  Third, depending on the \code{output} value, the function will compute RSA parameters.  If \code{output="lrd"}, the function will use \code{interpol} to select the appropriate DBase values for which the calculation of RSA parameters should be performed.  Then, a lateral root density and a total lateral root length will be calculated on intervals possessing a length equal to \code{int.length} and centred on each DBase value.  If \code{output="dtp"}, only the distances between neighbouring lateral roots on each mother root will be computed by the function.

If RSA data have been encoded with the Root System Markup Language (\cite{Lobet et al, 2015}), the path to the folder containing the rsml files can be specified using the \code{inputrsml} argument.  To run \code{latdist} efficiently, rsml files must have been saved with their appropriate extensions (.rsml).  If \code{is.null(inputrsml)=FALSE}, \code{latdist} will check if \code{rsml.date} and \code{rsml.connect} are specified and written in the suitable mode.  After reading the rsml files located in \code{inputrsml}, the RSA data associated with one root system in each rsml file will be converted into a data frame possessing the same structure as the rac files created by DART.  The \code{latdist} function will then use these data frames to compute lateral root length and density distribution.  As some rsml files do not provide connection information between daughter and mother roots, the \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root will be connected to the nearest point located on its mother root.  If \code{rsml.connect=FALSE} or \code{rsml.connect=TRUE}, the distance between a branching point and the parent root base (DBase) will be approximated by the distance between the point located on the mother root that is closest to the point starting a lateral root and the parent root base.  Using rsml files, only RSA data associated with roots possessing a branching order lower or equal to 5 will be used for the computation of RSA traits.
}
\value{
Returns a list including:

\item{root}{A list of data frames. Each element of the list is named as its corresponding rac/rsml file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system in the corresponding rac/rsml file), \bold{Ord} (the branching order), \bold{LatRootNum} (the number of lateral roots), \bold{FinalRootLength} (the root length at the last observation date expressed in \code{unitlength}), \bold{LatRootDensity} (the lateral root density calculated as LatRootNum/FinalRootLength and expressed in root(s)/\code{unitlength}).}

\item{res}{A list of secondary lists. Each element of the primary list (res) is named as its corresponding rac/rsml file. The secondary lists contain as much elements as roots constituting a vectorized root system. Each element of the secondary lists could be either a \code{NULL} value if a root do not have lateral roots or if no DBase value is comprised in [\code{int.length}/2, FinalRootLength-(\code{int.length}/2)], or a data frame made of two or three columns depending on the \code{output} value. 

\bold{If \code{output="lrd"}}: \bold{DBase} (the distance between the branching point and the parent root base expressed in \code{unitlength}), \bold{LRD} (Local lateral root density expressed in root(s)/\code{unitlength}), \bold{LRL} (Local lateral root length expressed in \code{unitlength}/\code{unitlength}).

\bold{If \code{output="dtp"}}: \bold{DBase} (the distance between the branching point and the parent root base expressed in \code{unitlength}), \bold{DTP} (the distance to the previous lateral root expressed in \code{unitlength}).}
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2015) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.
}
\author{Benjamin M. Delory, Guillaume Lobet, Loic Pages
}
\examples{
## Importation of an example dataset
data(ch7rac)

## Creation of a folder (DART_latrootdensity) inside a temporary directory to contain the rac
## example datasets created by DART
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