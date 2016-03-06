\name{archigrow}
\alias{archigrow}
\title{Computing Growth Rates and X-Y Plotting of Vectorized Root Systems}
\description{Calculation of growth rate matrices and X-Y plotting of vectorized root systems for selected observation dates using Data Analysis of Root Tracings (DART) output files. In the final plots, the colour code used for each link constituting the vectorized root systems depends on their corresponding growth rate value. More information can be found in \cite{Delory et al (2016)}.  
}
\usage{
archigrow(inputlie, inputtps, res=NULL, unitlength="px", unittime, 
          unitangle="d", rotation=0, numdate=NULL, finalscale=NULL, 
          coldyn, GRscale=NULL, main=NULL, xlab=NULL, ylab=NULL,
          xlim=NULL, ylim=NULL, ...)}
\arguments{
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputtps}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as \code{inputlie}.}
  \item{res}{\bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. This argument must be specified only if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that must be used by the function both to calculate the growth rate matrices and to scale the X and Y axes of the graphical outputs. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{unittime}{A character string specifying the unit of time used to express the data in the tps files created by DART (all character strings are valid for this argument).}
  \item{unitangle}{A character string specifying the unit that must be used by the function to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used by the function to perform a clockwise rotation of the root system (see details).}
  \item{numdate}{A numeric value or a vector of numeric values specifying the identification number of all the observation dates at which the root system architecture must be plotted (see details).}
  \item{finalscale}{A logical value that must be specified when \code{is.null(numdate)=FALSE}. If \code{TRUE}, the root system architecture at the selected observation dates will be plotted with the same X and Y-axis scales to that used to draw the root system architecture at the last observation date. If \code{FALSE}, the root system architecture at each observation date will be plotted using different X and Y-axis ranges of coordinates (\code{xlim} and \code{ylim}).}
  \item{coldyn}{A vector specifying the colour(s) that must be used to draw the root system architecture at the selected observation dates. The function will colour each link constituting a vectorized root system depending on its corresponding growth rate value. To do so, \code{archigrow} interpolates the colour(s) contained in \code{coldyn} between the minimum and the maximum values of the computed standardized growth rate matrix (if \code{is.null(GRscale)=TRUE}) or between \code{min(GRscale)} and \code{max(GRscale)} (if \code{is.null(GRscale)=FALSE}).}
  \item{GRscale}{A numeric vector (length=2) specifying the minimum and the maximum growth rate values (expressed in \code{unitlength}/\code{unittime}) that must be used by the function to plot each vectorized root system.}
  \item{main}{A character string giving a main title for the plot. When \code{is.null(main)=TRUE}, the default title for each plot is the name of the corresponding lie file.}
  \item{xlab}{A character string giving a label for the X axis. When \code{is.null(xlab)=TRUE}, the default label for the X axis is X (\code{unitlength}).}
  \item{ylab}{A character string giving a label for the Y axis. When \code{is.null(ylab)=TRUE}, the default label for the Y axis is Y (\code{unitlength}).}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot using the same unit as \code{unitlength}.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot using the same unit as \code{unitlength}.}
  \item{...}{Additional graphical parameters (see details)}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex root system architectures from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  At the end of the vectorization process, a final dataset can be exported as a table file (lie file) containing several attributes for each point used to construct the links constituting the vectorized root system.  A second table file (tps file) containing the root system age at each observation date can also be exported from DART.  Using these files created by DART, \code{archigrow} allows the X-Y plotting of vectorized root systems for selected observation dates with a colour code depending on the growth rate value of each link constituting the vectorized root systems.   

The \code{archigrow} function will successively locate and read the lie and tps files created by DART, compute a growth rate matrix and plot the vectorization results at selected observation dates for each analysed root system.  First, \code{archigrow} will check if all arguments have been written in the suitable mode.  If \code{res} and \code{rotation} are specified, the function will check if the numeric values required by these arguments are positive.  If \code{is.null(numdate)=FALSE}, the function will also automatically check if the numeric values in \code{numdate} are positive and sorted by increasing values.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputlie} and \code{inputtps} to locate and read the raw lie and tps files created by DART.  To run \code{archigrow} efficiently, DART files must have been saved with their appropriate extensions (.lie or .tps).  If \code{inputtps} contains a single tps file, it will be used by \code{archigrow} for each lie file located in \code{inputlie}.  If \code{inputps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of lie files in \code{inputlie} and corresponding lie and tps files must have the same name.  If an error occurs at this step, a message will be returned by the function.  Third, \code{archigrow} will compute a growth rate matrix containing for each root constituting the vectorized root system its corresponding growth rate value at each observation date.  Finally, \code{archigrow} will plot each vectorized root system located in \code{inputlie} at the observation dates specified by \code{numdate}. By default (when \code{is.null(numdate)=TRUE}), only the root system architecture at the last observation date will be plotted.  If \code{is.null(numdate)=FALSE}, \code{archigrow} will plot only the root system architecture for the selected observation dates.  

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In this case, the vectorized root system will depict a rotation compared with the natural plant root system.  To make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} that will be used by the function to rotate the vectorized root system clockwise before plotting and root system architecture parameters calculation.  

Additional graphical parameters can also be used to personalize the graphical outputs (see \code{\link{par}}).  By default, a legend showing the colour gradient used to draw each vectorized root system is shown at the right side of the final plot.  

Growth rates are calculated as follows:  
\itemize{
\item For the first observation date, it is calculated as the ratio of the root length to the root system age.  
\item For other observation dates (t), it is calculated as the difference between the root length at time t and t-1 divided by the difference between the root system age at time t and t-1.}
}
\value{
Returns the computed growth rate matrices in a list of data frames. Each element of the list is named as its corresponding lie file and contains the following columns: \bold{Root} (the identification number of each root constituting a vectorized root system), \bold{GR.DateX} (the calculated growth rate value of a single root at the observation date X expressed in \code{unitlength}/\code{unittime}).

Returns a X-Y plot for each vectorized root system located in \code{inputlie}.
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\seealso{\code{\link{par}}.
}
\examples{
## Importation of example datasets
data(ch7lie)
data(ch7tps)

## Creation of a folder (DART_archigrow) inside a temporary directory to contain the lie and tps 
## example datasets created by DART
dirDART <- paste(tempdir(), "/DART_archigrow", sep="")
dir.create(dirDART)

## Definition of the number of lie and tps files that will be saved inside dirDART 
n <- 1

## Creation of lie and tps files inside dirDART
for (i in 1:n) {
write.table(ch7lie, file=paste(dirDART,"/ch7_",i,".lie", sep=""), row.names=FALSE)
write.table(ch7tps, file=paste(dirDART,"/ch7_",i,".tps", sep=""), row.names=FALSE)}

## Use of archigrow to plot the root system architecture for selected observation dates
## Use the same X and Y-axis scales to that used for the last observation date
ch7_output_1 <- archigrow(inputlie=dirDART, inputtps=dirDART, res=75, unitlength="cm", 
unittime="day", coldyn=c("blue", "orange", "red"), numdate=c(15,30), finalscale=TRUE, 
las=1, bty="l", asp=1, xaxp=c(0,30,3), lwd=2)

unlink(dirDART, recursive=TRUE)
}