\name{trajectory}
\alias{trajectory}
\title{Computing Root Growth Directions and Trajectories}
\description{Calculation of root system architecture (RSA) parameters describing the directions and the trajectories of the root system growth using Data Analysis of Root Tracings (DART) output files.    
}
\usage{
trajectory(inputrac, inputlie, inputtps, res=NULL, unitlength="px", 
           unitangle="d", rotation=0, l.brangle, l.varangle, l.tipangle)}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as inputlie/inputtps.}
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as inputrac/inputtps.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as inputrac/inputlie.}
  \item{res}{A numeric value specifying the image resolution (dots per inch, dpi) of the pictures used in DART for the vectorization of the root systems. This argument must only be specified if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that will be used by the function to express the results. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimeters or \dQuote{mm} for millimeters.}
  \item{unitangle}{A character string specifying the unit that will be used by the function to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used by the function to perform a clockwise rotation of the root system (see details).}
  \item{l.brangle}{A numeric value specifying the minimum root length that will be used by the function to calculate the branching angle of a daughter root on its corresponding mother root (see details).}
  \item{l.varangle}{A numeric value specifying the distance that will be used by the function to place equidistantly spaced interpolated points along a single root to calculate statistical parameters (mean and standard deviation) describing the variation of the angle formed by the successive links with a horizontal direction vector (see details).}
  \item{l.tipangle}{A numeric value specifying the distance between the root tip and an interpolated point located along a single root. Each interpolated point will be used by the function to calculate a direction vector following the root growth direction near the apex prior to the calculation of the inclination of the root tip relative to horizontal at each observation date.}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex RSAs from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  After vectorization of a root system, three final data sets can be exported from DART as a table files: (1) a rac file containing several attributes for each individual root, notably the root length at each observation date, (2) a lie file containing several attributes for each point used to construct the links constituting the vectorized root system, and (3) a tps file containing the root system age at each observation date.  This R function allows the calculation of the branching angle of each daughter root on its corresponding mother root, statistical parameters describing the trajectory of each root consisting a vectorized root system, and the root tip angles relative to horizontal at each observation date.  

The function \code{trajectory} will successively locate and read the rac, lie and tps files created by DART and calculate RSA parameters.  First, \code{trajectory} will check if all arguments have been written in the suitable mode.  If \code{res}, \code{rotation}, \code{l.brangle}, \code{l.varangle} and \code{l.tipangle} are specified, the function will check if the numeric values required by these arguments are positive.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputrac}, \code{inputlie} and \code{inputtps} to locate and read the raw rac, lie and tps files created by DART respectively.  To run \code{trajectory} efficiently, the DART files must have been saved with their appropriate extensions (.rac, .lie or .tps).  If \code{inputtps} contains a single tps file, it will be used by \code{archigrow} for each lie/rac file located in \code{inputlie}/\code{inputrac}.  Moreover, the number of lie files in \code{inputlie} must be equal to the number of rac files in \code{inputrac} and corresponding rac and lie files must have the same name.  If \code{inputps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of rac/lie files in \code{inputrac}/\code{inputlie} and corresponding rac/lie and tps files must have the same name.  If an error occurs at this step, a message will be returned by the function. Third, the function will compute RSA parameters.

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In that case, the vectorized root system will depict a rotation compared to the natural plant root system.  To make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} that will be used by the function to rotate the vectorized root system clockwise before RSA parameters calculation.    

Hereafter, we will consider that the normal vector that is orthogonal to a direction vector (a,b) will have the following coordinates: (b,-a).  

A daughter root will be considered to grow at the left side of the mother root if the scalar product between a direction vector (=vector A) going from the branching point (X0) to the following point (Xd) on the daughter root and a vector (=vector B) normal to a direction vector (=vector C) going from the branching point (X0) to the following point on the mother root (Xm) is positive.  If the scalar product between A and B is negative, the daughter root will be considered to grow at the right side of the mother root.  If the scalar product between A and B is equal to zero, the calculation of the scalar product between A and B will be performed again using a vector A going from the branching point (X0) to a point located on the daughter root at Xd<-Xd+1 until the calculated scalar product is different from zero.  If A is always normal to B, a random lateral root growth orientation (either left or right) will be defined.  

The branching angle of a daughter root on its mother root is approximated by the angle between two direction vectors going from the branching point to an interpolated point located at a distance \code{l.brangle} from the branching point on the mother or on the daughter root respectively.  The distance \code{l.brangle} between a branching point and an interpolated point is measured along the mother and the daughter roots.  A primary root having no mother root, a \code{NA} value will be returned for Branching.Angle.  If a root has a final length lower than \code{l.brangle}, no branching angle can be calculated and the function will return a \code{NA} value for Branching.Angle.  

The algorithm used to calculate the variation of the angle formed by a horizontal direction vector and the successive links of a single root comprises three main steps.  First, the function will use \code{l.varangle} to interpolate equidistantly spaced points along each single root consisting a vectorized root system.  Second, the angle formed between each successive link of a single root and a horizontal direction vector is calculated.  Then, the mean and the standard deviation of that angle are calculated for each single root.  If a root has a final length lower than 2*\code{l.varangle}, no angle variation can be calculated along that root and the function will return a \code{NA} value for Mean.AngleVar and SD.AngleVar.  If a root has a final length lower than 3*\code{l.varangle}, only the mean of the angle variation along that root can be calculated and the function will return a \code{NA} value for SD.AngleVar.  
}
\value{
Returns a list including:

\item{root}{A list of data frames. Each element of the list is named as its corresponding rac/lie file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system), \bold{Mother} (the parent root identification number), \bold{Ord} (the branching order), \bold{DBase} (the distance between the branching point to the parent root base expressed in \code{unitlength}), \bold{DApp} (the date of emergence expressed in the same unit as the root system age in the tps file(s)), \bold{FinalRootLength} (the root length at the last observation date expressed in \code{unitlength}), \bold{Orientation} (the direction of the lateral root growth), \bold{Branching.Angle} (the branching angle of a daughter root on its mother root expressed in \code{unitangle}), \bold{Mean.AngleVar} (the mean of the variation of the angle calculated between a horizontal direction vector and each successive link constructed using equidistantly spaced interpolated points along each root of a vectorized root system; the calculated mean is expressed in \code{unitangle}), \bold{SD.AngleVar} (the standard deviation of the variation of the angle calculated between a horizontal direction vector and each successive link constructed using equidistantly spaced interpolated points along each root of a vectorized root system; the calculated standard deviation is expressed in \code{unitangle}).}

\item{tip}{A list of data frames. Each element of the list is named as its corresponding rac/lie file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system), \bold{Ang.DateX} (the calculated root tip angle relative to horizontal at the observation date X expressed in \code{unitangle}).}
}
\references{
Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, \bold{326}, 261--273.
}
\author{Benjamin M. Delory, Loic Pages
}
\examples{
## Importation of example data sets
data(ch7rac)
data(ch7lie)
data(ch7tps)

## Creation of a folder (DART_trajectory) inside a temporary directory to contain the rac and lie
## example datasets created by DART
dirDART <- paste(tempdir(), "/DART_trajectory", sep="")
dir.create(dirDART)

## Definition of the number of rac, tps and lie files that will be saved inside dirDART
n <- 1

## Creation of rac, tps and lie files inside dirDART
for (i in 1:n) {
write.table(ch7rac, file=paste(dirDART,"/ch7_",i,".rac", sep=""), row.names=FALSE)
write.table(ch7lie, file=paste(dirDART,"/ch7_",i,".lie", sep=""), row.names=FALSE)
write.table(ch7tps, file=paste(dirDART,"/ch7_",i,".tps", sep=""), row.names=FALSE)}

## Use of trajectory to compute RSA parameters from DART output files
ch7_output_7 <- trajectory(inputrac=dirDART, inputlie=dirDART, inputtps=dirDART, res=75, 
unitlength="cm", l.brangle=1, l.varangle=1, l.tipangle=0.5)

## Plotting the results for ch7_1
## Distribution of the calculated insertion angles
hist(ch7_output_7$root$ch7_1$Branching.Angle, breaks=seq(from=0, to=180, by=5),
main="ch7_1-Branching angle distribution", xlab="Branching angle (d)", las=1,
xaxp=c(0,150,15), xlim=c(0,150))

## Evolution of the root tip angle for the primary root (black) and the two longest lateral roots
## (green and red)
date<-c(1:31)
plot(date, ch7_output_7$tip$ch7_1[1,2:ncol(ch7_output_7$tip$ch7_1)], type="l", lwd=2, bty="l",
las=1, ylim=c(0,180), ylab="Root tip angle (d)", xlab="Time (Num)", col="black",
main="Root tip angle", yaxp=c(0,180,18))
lines(date, ch7_output_7$tip$ch7_1[206,2:ncol(ch7_output_7$tip$ch7_1)], lwd=2, col="green")
lines(date, ch7_output_7$tip$ch7_1[221,2:ncol(ch7_output_7$tip$ch7_1)], lwd=2, col="red")
abline(h=90)

unlink(dirDART, recursive=TRUE)
}