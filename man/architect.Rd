\name{architect}
\alias{architect}
\title{Computing Traits Describing The Global Root System Architecture}
\description{Calculation of common root system architecture (RSA) traits from Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). More information can be found in \cite{Delory et al (2016)}.      
}
\usage{
architect(inputrac=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px",
          rsml.date=NULL, rsml.connect=FALSE, rootdiv=1)}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as \code{inputtps}/\code{inputrsml}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as \code{inputrac}/\code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the rsml files. This argument could be the same as \code{inputrac}/\code{inputtps}.}
  \item{res}{\bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. This argument must be specified only if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that must be used by the function to express the results in the final data frame. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{rsml.date}{A numeric value that must be specified when \code{is.null(inputrsml)=FALSE}. This argument specifies the root system age (the unit of time is defined by the user) and will be used only for RSA data encoded with the Root System Markup Language.}
  \item{rsml.connect}{A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. This argument will be used only for RSA data encoded with the Root System Markup Language. If \code{TRUE}, the lateral roots will be connected to their parent root if no connection information can be found in the rsml file(s) (see details). If \code{FALSE}, no connection will be made between daughter and mother roots if no connection information can be found in the rsml file(s).}
  \item{rootdiv}{Either a numeric value or a vector of numeric values defining a class number or the class limits that must be used by the function to divide the first-order root in zones. Using this argument, \code{architect} will also compute the distribution of the secondary root number, length and density on the first-order root.}
}
\details{

\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex RSA from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  After vectorization of a root system, a final dataset can be exported as a table file (rac file) containing several attributes for each individual root, notably the root length at each observation date.  A second table file (tps file) containing the root system age at each observation date can also be exported from DART.  These two datasets can finally be used to calculate derived RSA traits.  In this context, \code{architect} is an R function allowing the calculation of common RSA traits at each observation date for each DART rac file: the total root system length, the length of the first-order root, the total number of lateral roots, the total lateral root length, the length and the number of lateral roots for each branching order, the mean lateral root length for each branching order, the density of secondary roots on the first-order root, the first-order and the lateral root growth rates, as well as the growth rate of the entire root system.  When needed, the distribution of the secondary root number, length and density on the first-order root can also be calculated (via \code{rootdiv}).  

The \code{architect} function will successively locate and read the rac and tps files created by DART and calculate RSA traits.  First, \code{architect} will check if all arguments have been written in the suitable mode.  If \code{res} is specified, the function will check if the numeric value required by this argument is positive.  The function will also automatically check if the numeric elements in \code{rootdiv} are superior or equal to zero and sorted by increasing values.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputrac} and \code{inputtps} to locate and read the raw rac and tps files created by DART, respectively.  To run \code{architect} efficiently, DART files must have been saved with their appropriate extensions (.rac or .tps).  Third, the function will create vectors and matrices for RSA traits calculation.  If \code{inputtps} contains a single tps file, this will be used by \code{architect} for each rac file located in \code{inputrac}.  If \code{inputtps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of rac files in \code{inputrac} and corresponding rac and tps files must have the same name.  If an error occurs at this step, a message will be returned by the function.

If RSA data have been encoded with the Root System Markup Language (\cite{Lobet et al, 2015}), the path to the folder containing the rsml files can be specified using the \code{inputrsml} argument.  To run \code{architect} efficiently, rsml files must have been saved with their appropriate extensions (.rsml).  If \code{is.null(inputrsml)=FALSE}, \code{architect} will check if \code{rsml.date} and \code{rsml.connect} are specified and written in the suitable mode.  After reading the rsml files located in \code{inputrsml}, the RSA data associated with one root system in each rsml file will be converted into a data frame possessing the same structure as the rac files created by DART.  The \code{architect} function will then use these data frames to calculate RSA traits describing the global architecture of each root system.  The argument \code{rsml.date} will be used to create a tps file for each root system located in \code{inputrsml}.  As some rsml files do not provide connection information between daughter and mother roots, the \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root will be connected to the nearest point located on its mother root.  If \code{rsml.connect=FALSE} or \code{rsml.connect=TRUE}, the distance between a branching point and the parent root base (DBase) will be approximated by the distance between the point located on the mother root that is closest to the point starting a lateral root and the parent root base.  Using rsml files, only RSA data associated with roots possessing a branching order lower or equal to 5 will be used for the computation of RSA traits.  


Growth rates are calculated as follows:
\itemize{
\item \emph{First-order root growth rate}: for the first observation date, it is calculated as the ratio of the first-order root length to the root system age. For other observation dates (t), it is calculated as the difference between the first-order root length at time t and t-1 divided by the difference between the root system age at time t and t-1.
\item \emph{Lateral root growth rates}: for the first observation date, it is calculated as the ratio of the lateral root length to the root system age. For other observation dates (t), it is calculated as the difference between the lateral root length at time t and t-1 divided by the difference between the root system age at time t and t-1.
\item \emph{Growth rate of the root system}: for the first observation date, it is calculated as the ratio of the total root system length to the root system age. For other observation dates (t), it is calculated as the difference between the total root system length at time t and t-1 divided by the difference between the root system age at time t and t-1.}
}
\value{
Returns a data frame summarizing all RSA traits calculated for each rac/rsml file located in \code{inputrac}/\code{inputrsml}. 

The data frame contains the following columns: \bold{FileName} (the name of the rac/rsml files), \bold{Time} (the root system age expressed in the same unit (=unittime) as in the tps file(s) (for rac files) or in the same unit as \code{rsml.date} (for rsml files)), \bold{TRL} (the total root system length expressed in \code{unitlength}), \bold{GRTR} (the growth rate of the root system expressed in \code{unitlength}/unittime), \bold{L1R} (the first-order root length expressed in \code{unitlength}), \bold{GR1R} (the growth rate of the first-order root expressed in \code{unitlength}/unittime), \bold{TNLR} (the total number of lateral roots), \bold{TLRL} (the total lateral root length expressed in \code{unitlength}), \bold{NxLR} (the total number of lateral roots with a branching order equal to x), \bold{LxLR} (the total length of lateral roots with a branching order equal to x expressed in \code{unitlength}), \bold{MLxLR} (the mean length of lateral roots with a branching order equal to x expressed in \code{unitlength}), \bold{GRxL} (the growth rate of lateral roots with a branching order equal to x expressed in \code{unitlength}/unittime), \bold{D2LR} (the density of secondary roots on the first-order root expressed in root(s)/\code{unitlength}).

\bold{If \code{rootdiv} is a numeric value > 1}, the following columns will be added to the data frame: \bold{N2LR.LayerX} (the secondary root number in layer X), \bold{L2LR.LayerX} (the total secondary root length in layer X expressed in \code{unitlength}), \bold{D2LR.LayerX} (the secondary root density on the first-order root in the layer X expressed in root(s)/\code{unitlength}).

\bold{If \code{rootdiv} is a vector of numeric values}, the following columns will be added to the data frame: \bold{N2LR.XtoY} (the secondary root number in a layer delimited by X and Y on the first-order root), \bold{L2LR.XtoY} (the total secondary root length in a layer delimited by X and Y on the first-order root expressed in \code{unitlength}), \bold{D2LR.XtoY} (the secondary root density in a layer delimited by X and Y on the first-order root expressed in root(s)/\code{unitlength}).
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625. 
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Guillaume Lobet, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\examples{
## Importation of example datasets
data(ch7rac)
data(ch7tps)

## Creation of a folder (DART_architect) inside a temporary directory to contain the rac and tps
## example datasets created by DART
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

## Use of the rootdiv argument to compute the distribution of the secondary root number, length
## and density in the [0-15[, [15-30[, [30-45[ and [45-60[ first centimeters of the first-order
## root
ch7_output_3 <- architect(inputrac=dirDART, inputtps=dirDART, res=75, unitlength="cm",
rootdiv=c(0,15,30,45,60))

unlink(dirDART, recursive=TRUE)
}