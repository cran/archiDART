\name{trajectory}
\alias{trajectory}
\title{Computing Root Growth Directions and Trajectories}
\description{Calculation of root system architecture (RSA) parameters describing the directions and the trajectories of the root system growth using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). More information can be found in \cite{Delory et al (2016)}.    
}
\usage{
trajectory(inputrac=NULL, inputlie=NULL, inputtps=NULL, inputrsml=NULL, res=NULL,
          unitlength="px", unitangle="d", rotation=0, l.brangle, l.curv, l.tipangle,
          rsml.date=NULL)}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as \code{inputlie}/\code{inputtps}/\code{inputrsml}.}
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputrac}/\code{inputtps}/\code{inputrsml}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as \code{inputrac}/\code{inputlie}/\code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the rsml files. This argument could be the same as \code{inputrac}/\code{inputlie}/\code{inputtps}.}
  \item{res}{\bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. This argument must be specified only if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that must be used by the function to express the results. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{unitangle}{A character string specifying the unit that must be used by the function to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used by the function to perform a clockwise rotation of the root system (see details).}
  \item{l.brangle}{A numeric value specifying the minimum root length that must be used by the function to calculate the branching angle of a daughter root on its corresponding mother root (see details).}
  \item{l.curv}{A numeric value specifying the distance that must be used by the function to place equidistantly spaced interpolated points along each root before calculating the angles between the direction vectors of the successive links constituting a vectorized root (see details).}
  \item{l.tipangle}{A numeric value specifying the distance between the root tip and an interpolated point located along a root. Each interpolated point will be used by the function to calculate a direction vector following the root growth direction near the apex prior to the calculation of the root tip angle relative to a vertical direction vector at each observation date.}
  \item{rsml.date}{A numeric value that must be specified when \code{is.null(inputrsml)=FALSE}. This argument specifies the root system age (the unit of time is defined by the user) and will be used only for RSA data encoded with the Root System Markup Language.}
}
\details{
\cite{Le Bot et al (2010)} presented a free and open-access software (Data Analysis of Root Tracings - DART) allowing the analysis of complex RSA from captured images, particularly across time series.  Using this software, a user has to manually identify roots as a set of links.  After vectorization of a root system, three final datasets can be exported from DART as table files: (1) a rac file containing several attributes for each individual root, notably the root length at each observation date, (2) a lie file containing several attributes for each point used to construct the links constituting the vectorized root system, and (3) a tps file containing the root system age at each observation date.  This R function allows the calculation of the branching angle of each daughter root on its corresponding mother root, statistical parameters describing the curvature of each root constituting a vectorized root system, and the root tip angles relative to a vertical direction vector at each observation date.  

The \code{trajectory} function will successively locate and read the rac, lie and tps files created by DART and calculate RSA traits.  First, \code{trajectory} will check if all arguments have been written in the suitable mode.  If \code{res}, \code{rotation}, \code{l.brangle}, \code{l.curv} and \code{l.tipangle} are specified, the function will check if the numeric values required by these arguments are positive.  If an error occurs for one argument, the code stops and the corresponding error message will be returned by the function.  Second, the function will use \code{inputrac}, \code{inputlie} and \code{inputtps} to locate and read the raw rac, lie and tps files created by DART, respectively.  To run \code{trajectory} efficiently, the DART files must have been saved with their appropriate extensions (.rac, .lie or .tps).  If \code{inputtps} contains a single tps file, it will be used by \code{trajectory} for each lie/rac file located in \code{inputlie}/\code{inputrac}.  Moreover, the number of lie files in \code{inputlie} must be equal to the number of rac files in \code{inputrac} and corresponding rac and lie files must have the same name.  If \code{inputps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of rac/lie files in \code{inputrac}/\code{inputlie} and corresponding rac/lie and tps files must have the same name.  If an error occurs at this step, a message will be returned by the function. Third, the function will compute RSA traits.

If RSA data have been encoded with the Root System Markup Language (\cite{Lobet et al, 2015}), the path to the folder containing the rsml files can be specified using the \code{inputrsml} argument.  To run \code{trajectory} efficiently, rsml files must have been saved with their appropriate extensions (.rsml).  If \code{is.null(inputrsml)=FALSE}, \code{trajectory} will check if \code{rsml.date} is specified and written in the suitable mode.  After reading the rsml files located in \code{inputrsml}, the RSA data associated with one root system in each rsml file will be converted into two data frames possessing the same structure as the rac and lie files created by DART, respectively.  The \code{trajectory} function will then use these data frames to compute RSA traits describing root growth directions and trajectories.  The argument \code{rsml.date} will be used to create a tps file for each root system located in \code{inputrsml}.  By default, if no connection information between daughter and mother roots exists in a rsml file, the function will connect the lateral roots with their correponding mother root.  To do so, each point starting a lateral root will be connected to the nearest point located on its mother root.  The distance between a branching point and the parent root base (DBase) will be approximated by the distance between the point located on the mother root that is closest to the point starting a lateral root and the parent root base.  Using rsml files, only RSA data associated with roots possessing a branching order lower or equal to 5 will be used for the computation of RSA traits.

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In this case, the vectorized root system will depict a rotation compared with the natural plant root system.  To make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} that will be used by the function to rotate the vectorized root system clockwise before RSA parameters calculation.    

Hereafter, we will consider that the normal vector that is orthogonal to a direction vector (a,b) will have the following coordinates: (b,-a).  

A daughter root will be considered to grow at the left side of the mother root if the scalar product between a direction vector (=vector A) going from the branching point (X0) to the following point (Xd) on the daughter root and a vector (=vector B) normal to a direction vector (=vector C) going from the branching point (X0) to the following point on the mother root (Xm) is positive.  If the scalar product between A and B is negative, the daughter root will be considered to grow at the right side of the mother root.  If the scalar product between A and B is equal to zero, the calculation of the scalar product between A and B will be performed again using a vector A going from the branching point (X0) to a point located on the daughter root at Xd<-Xd+1 until the calculated scalar product is different from zero.  If A is always normal to B, a random lateral root growth orientation (either left or right) will be defined.  

The branching angle of a daughter root on its mother root is approximated by the angle between two direction vectors going from the branching point to an interpolated point located at a distance \code{l.brangle} from the branching point on the mother or on the daughter root, respectively.  The distance \code{l.brangle} between a branching point and an interpolated point is measured along the mother and the daughter roots.  A first-order root having no mother root, a \code{NA} value will be returned for Branching.Angle. If a root has a final length lower than \code{l.brangle}, no branching angle can be calculated and the function will return a \code{NA} value for Branching.Angle.  

The curvature of each root is evaluated by the mean and the standard deviation of the local angles formed by the direction vectors of the successive links constituting a vectorized root.  First, the function will use \code{l.curv} to interpolate equidistantly spaced points along each root constituting a vectorized root system.  Second, the angles between the direction vectors of the successive links of each root are calculated. Then, the mean and the standard deviation of the calculated angles are determined for each root.  If a root has a final length lower than 2*\code{l.curv}, the function will return a \code{NA} value for Mean.Curv and SD.Curv. If a root has a final length lower than 3*\code{l.curv}, the function will return a \code{NA} value for SD.Curv as no standard deviation can be calculated on a single angle value.  
}
\value{
Returns a list including:

\item{root}{A list of data frames. Each element of the list is named as its corresponding rac/lie/rsml file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system), \bold{Mother} (the parent root identification number), \bold{Ord} (the branching order), \bold{DBase} (the distance between the branching point to the parent root base expressed in \code{unitlength}), \bold{DApp} (the date of emergence expressed in the same unit as the root system age in the tps file(s) (for rac/lie files) or in rsml.date (for rsml files)), \bold{FinalRootLength} (the root length at the last observation date expressed in \code{unitlength}), \bold{Tortuosity} (the ratio of the final root length to the Euclidean distance between the branching point and the apex of the root), \bold{Orientation} (the direction of the lateral root growth), \bold{Branching.Angle} (the branching angle of a daughter root on its mother root expressed in \code{unitangle}), \bold{Mean.Curv} (the mean of the local angle values calculated between the direction vectors of the successive links constructed using equidistantly spaced interpolated points along each root of a vectorized root system; the calculated mean is expressed in \code{unitangle}), \bold{SD.Curv} (the standard deviation of the local angle values calculated between the direction vectors of the successive links constructed using equidistantly spaced interpolated points along each root of a vectorized root system; the calculated standard deviation is expressed in \code{unitangle}.}

\item{tip}{A list of data frames. Each element of the list is named as its corresponding rac/lie/rsml file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system), \bold{Ang.DateX} (the calculated root tip angle relative to a vertical direction vector at the observation date X expressed in \code{unitangle}).}
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.
}
\author{Benjamin M. Delory, Guillaume Lobet, Loic Pages
}
\examples{
## Importation of example datasets
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
unitlength="cm", l.brangle=1, l.curv=1, l.tipangle=0.5)

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