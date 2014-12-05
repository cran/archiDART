\name{archidraw}
\alias{archidraw}
\title{X-Y Plotting of Vectorized Root Systems}
\description{X-Y plotting of vectorized root systems for selected observation dates using Data Analysis of Root Tracings (DART) output files.  
}
\usage{
archidraw(inputlie, res=NULL, unitlength="px", unitangle="d", rotation=0, numdate=NULL,
          finalscale=NULL, coldate=par("col"), ltydate=par("lty"), lwddate=par("lwd"), 
          main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, ...)
}
\arguments{
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART.}
  \item{res}{A numeric value specifying the image resolution (dots per inch, dpi) of the pictures used in DART for the vectorization of the root systems. This argument must only be specified if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that will be used by the function to scale the X and Y axes of the graphical outputs. The value acceptable for this argument could either be \dQuote{px} for pixels, \dQuote{cm} for centimeters or \dQuote{mm} for millimeters.}
  \item{unitangle}{A character string specifying the unit that will be used by the function to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that will be used by the function to perform a clockwise rotation of the root system.}
  \item{numdate}{A numeric value or a vector of numeric values specifying the identification number of all observation dates at which the root system architecture must be plotted (see details).}
  \item{finalscale}{A logical value that must be specified when \code{is.null(numdate)=FALSE}. If \code{TRUE}, the root system architecture at the selected observation dates will be plotted with the same X and Y-axis scales to that used to draw the root system architecture at the last observation date. If \code{FALSE}, the root system architecture at each observation date will be plotted using different X and Y-axis ranges of coordinates (\code{xlim} and \code{ylim}).}
  \item{coldate}{A vector specifying the colour(s) that will be used to draw the root system architecture at the selected observation dates (see details).}
  \item{ltydate}{A vector specifying the line type(s) that will be used to draw the root system architecture at the selected observation dates (see details).}
  \item{lwddate}{A vector specifying the line width(s) that will be used to draw the root system architecture at the selected observation dates (see details).}
  \item{main}{A character string giving a main title for the plot. When \code{main=NULL}, the default title for each plot is the name of the corresponding lie file.}
  \item{xlab}{A character string giving a label for the X axis. When \code{xlab=NULL}, the default label for the X axis is X (\code{unitlength}).}
  \item{ylab}{A character string giving a label for the Y axis. When \code{ylab=NULL}, the default label for the Y axis is Y (\code{unitlength}).}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot using the same unit as \code{unitlength}.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot using the same unit as \code{unitlength}.}
  \item{...}{Additional graphical parameters (see details).}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex root system architectures from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  At the end of the vectorization process, a final data set can be exported as a table file (lie file) containing several attributes for each point used to construct the links constituting the vectorized root system.  Using these lie files created by DART, \code{archidraw} allows the X-Y plotting of vectorized root systems for selected observation dates.

The function \code{archidraw} will locate and read the lie files created by DART and plot the vectorization results at selected observation dates for each analysed root system.  First, \code{archidraw} will check if all arguments have been written in the suitable mode.  If \code{res} and \code{rotation} are specified, the function will check if the numeric values required by these arguments are positive.  If \code{is.null(numdate)=FALSE}, the function will also automatically check if the numeric values in \code{numdate} are positive and sorted by increasing values.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputlie} to locate and read the raw lie files created by DART.  To run \code{archidraw} efficiently, DART files must have been saved with their appropriate extensions (.lie).  Finally, \code{archidraw} will plot each vectorized root system located in \code{inputlie} at the observation dates specified by \code{numdate}.  By default (when \code{numdate=NULL}), only the root system architecture of the last observation date will be plotted.  If \code{is.null(numdate)=FALSE}, \code{archidraw} will only plot the root system architecture for the selected observation dates.  If the number of elements specified in \code{coldate}/\code{ltydate}/\code{lwddate} is lower than \code{max(Date)} in the lie files (if \code{numdate=NULL}) or \code{max(numdate)} (if \code{is.null(numdate)=FALSE}), \code{archidraw} will replicate elements stored in \code{coldate}/\code{ltydate}/\code{lwddate} to draw the vectorization results.  However when the number of elements specified in \code{coldate}/\code{ltydate}/\code{lwddate} is greater than \code{max(Date)} in the lie files (if \code{numdate=NULL}) or \code{max(numdate)} (if \code{is.null(numdate)=FALSE}), only the first corresponding elements in \code{coldate}/\code{ltydate}/\code{lwddate} will be used by \code{archidraw} to plot the root system architecture for each file located in \code{inputlie}.  

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In that case, the vectorized root system will depict a rotation compared to the natural plant root system.  To make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} that will be used by the function to rotate the vectorized root system clockwise before plotting.  

Additional graphical parameters can also be used to personalize the graphical outputs (see \code{\link{par}}).  
}
\value{Returns a X-Y plot for each vectorized root system located in \code{inputlie}.
}
\references{
Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\seealso{\code{\link{par}}.
}
\examples{
## Importation of an example data set
data(ch7lie)

## Creation of a folder (called DART_archidraw) inside a temporary directory to contain lie 
## example datasets created by DART
dirDART <- paste(tempdir(), "/DART_archidraw", sep="")
dir.create(dirDART)

## Definition of the number of lie files that will be saved inside dirDART 
n <- 1

## Creation of lie files inside dirDART
for (i in 1:n) {
write.table(ch7lie, file=paste(dirDART,"/ch7_",i,".lie", sep=""), row.names=FALSE)}

## Use of archidraw to plot the root system architecture in black for the last observation date
archidraw(inputlie=dirDART, res=75, unitlength="cm", lwddate=2, las=1, bty="l", asp=1,
xaxp=c(0,30,3), yaxp=c(0,90,9))


## Use of archidraw to plot the root system architecture at the last observation date 
## Use a specific colour for each link related to a single observation date
archidraw(inputlie=dirDART, res=75, unitlength="cm", coldate=rainbow(31), lwddate=2, las=1,
bty="l", asp=1, xaxp=c(0,30,3), yaxp=c(0,90,9))

## Use of archidraw to plot the root system architecture at selected observation dates 
## Use the same X and Y-axis scales to that used for the last observation date and a specific colour
## for each observation date
archidraw(inputlie=dirDART, res=75, unitlength="cm", numdate=c(15,30), finalscale=TRUE, 
coldate=rainbow(31), lwddate=2, las=1, bty="l", asp=1, xaxp=c(0,30,3), yaxp=c(0,90,9))

unlink(dirDART, recursive=TRUE)
}