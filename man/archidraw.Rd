\name{archidraw}
\alias{archidraw}
\title{X-Y Plotting of Vectorized Root Systems}
\description{X-Y plotting of vectorized root systems for selected observation dates using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML).  
}
\usage{
archidraw(inputlie=NULL, inputrsml=NULL, res=NULL, unitlength="px", unitangle="d", 
          rotation=0, rsml.connect=FALSE, numdate=NULL, finalscale=NULL,
          coldate=par("col"), ltydate=par("lty"), lwddate=par("lwd"), 
          main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, ...)
}
\arguments{
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the rsml files. This argument could be the same as \code{inputlie}.}
  \item{res}{A numeric value specifying the image resolution (dots per inch, dpi) of the pictures used for the vectorization of the root systems. This argument must be specified only if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that must be used by the function to scale the X and Y axes of the graphical outputs. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{unitangle}{A character string specifying the unit that must be used by the function to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used by the function to perform a clockwise rotation of the root system.}
  \item{rsml.connect}{A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. This argument will be used only for RSA data encoded with the Root System Markup Language. If \code{TRUE}, the lateral roots will be connected to their parent root if no connection information can be found in the rsml file(s) (see details). If \code{FALSE}, no connection will be made between daughter and mother roots if no connection information can be found in the rsml file(s).}
  \item{numdate}{A numeric value or a vector of numeric values specifying the identification number of all the observation dates at which the root system architecture must be plotted (see details).}
  \item{finalscale}{A logical value that must be specified when \code{is.null(numdate)=FALSE}. If \code{TRUE}, the root system architecture at the selected observation dates will be plotted with the same X and Y-axis scales to that used to draw the root system architecture at the last observation date. If \code{FALSE}, the root system architecture at each observation date will be plotted using different X and Y-axis ranges of coordinates (\code{xlim} and \code{ylim}).}
  \item{coldate}{A vector specifying the colour(s) that must be used to draw the root system architecture at the selected observation dates (see details).}
  \item{ltydate}{A vector specifying the line type(s) that must be used to draw the root system architecture at the selected observation dates (see details).}
  \item{lwddate}{A vector specifying the line width(s) that must be used to draw the root system architecture at the selected observation dates (see details).}
  \item{main}{A character string giving a main title for the plot. When \code{is.null(main)=TRUE}, the default title for each plot is the name of the corresponding lie/rsml file.}
  \item{xlab}{A character string giving a label for the X axis. When \code{is.null(xlab)=TRUE}, the default label for the X axis is X (\code{unitlength}).}
  \item{ylab}{A character string giving a label for the Y axis. When \code{is.null(ylab)=TRUE}, the default label for the Y axis is Y (\code{unitlength}).}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot using the same unit as \code{unitlength}.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot using the same unit as \code{unitlength}.}
  \item{...}{Additional graphical parameters (see details).}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex root system architectures from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  At the end of the vectorization process, a final dataset can be exported as a table file (lie file) containing several attributes for each point used to construct the links constituting the vectorized root system.  Using these lie files created by DART, \code{archidraw} allows the X-Y plotting of vectorized root systems for selected observation dates.

The \code{archidraw} function will locate and read the lie files created by DART and plot the vectorization results at selected observation dates for each analysed root system.  First, \code{archidraw} will check if all arguments have been written in the suitable mode.  If \code{res} and \code{rotation} are specified, the function will check if the numeric values required by these arguments are positive.  If \code{is.null(numdate)=FALSE}, the function will also automatically check if the numeric values in \code{numdate} are positive and sorted by increasing values. If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputlie} to locate and read the raw lie files created by DART.  To run \code{archidraw} efficiently, DART files must have been saved with their appropriate extensions (.lie).  Finally, \code{archidraw} will plot each vectorized root system located in \code{inputlie} at the observation dates specified by \code{numdate}.  By default (when \code{is.null(numdate)=TRUE}), only the root system architecture at the last observation date will be plotted.  If \code{is.null(numdate)=FALSE}, \code{archidraw} will plot only the root system architecture for the selected observation dates.  If the number of elements specified in \code{coldate}/\code{ltydate}/\code{lwddate} is lower than \code{max(Date)} in the lie files (if \code{is.null(numdate)=TRUE}) or \code{max(numdate)} (if \code{is.null(numdate)=FALSE}), \code{archidraw} will replicate the elements stored in \code{coldate}/\code{ltydate}/\code{lwddate} to draw the vectorization results.  When the number of elements specified in \code{coldate}/\code{ltydate}/\code{lwddate} is greater than \code{max(Date)} in the lie files (if \code{is.null(numdate)==TRUE}) or \code{max(numdate)} (if \code{is.null(numdate)=FALSE}), only the first corresponding elements in \code{coldate}/\code{ltydate}/\code{lwddate} will be used by \code{archidraw} to plot the root system architecture for each file located in \code{inputlie}.  

If RSA data have been encoded with the Root System Markup Language (\cite{Lobet et al, 2015}), the path to the folder containing the rsml files can be specified using the \code{inputrsml} argument.  To run \code{archidraw} efficiently, rsml files must have been saved with their appropriate extensions (.rsml).  If \code{is.null(inputrsml)=FALSE}, \code{archidraw} will check if \code{rsml.connect} is specified and written in the suitable mode.  After reading the rsml files located in \code{inputrsml}, the RSA data associated with one root system in each rsml file will be converted into a data frame possessing the same structure as the lie files created by DART.  The \code{archidraw} function will then use these data frames to plot the root system architectures encoded with the RSML format.   As some rsml files do not provide connection information between daughter and mother roots, the \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root will be connected to the nearest point located on its mother root.  Using rsml files, only RSA data associated with roots possessing a branching order lower or equal to 5 will be used for the computation of graphical outputs.

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In this case, the vectorized root system will depict a rotation compared with the natural plant root system.  In order to make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} that will be used by the function to rotate the vectorized root system clockwise before plotting.  

Additional graphical parameters can also be used to personalize the graphical outputs (see \code{\link{par}}).  
}
\value{Returns a X-Y plot for each vectorized root system located in \code{inputlie}/\code{inputrsml}.
}
\references{
Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, \bold{167}, 617--627.
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Guillaume Lobet, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\seealso{\code{\link{par}}.
}
\examples{
## Importation of an example dataset
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