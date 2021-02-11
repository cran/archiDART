\name{plot.barcode}
\alias{plot.barcode}
\title{Plot the Persistence Barcode}
\description{Plot the persistence barcode of the topology of a root system. More information can be found in \cite{Li et al (2017)} and \cite{Delory et al (2018)}.}
\usage{
\method{plot}{barcode}(x, y = NULL, type = "p",  xlim = NULL, ylim = NULL,
        log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
        ann = par("ann"), axes = TRUE, frame.plot = axes,
        panel.first = NULL, panel.last = NULL, asp = NA,
        xgap.axis = NA, ygap.axis = NA, ...)
}
\arguments{
  \item{x}{A barcode object. Mandatory.}
  \item{y}{See \code{\link{plot.default}}}
  \item{type}{See \code{\link{plot.default}}}
  \item{xlim}{See \code{\link{plot.default}}}
  \item{ylim}{See \code{\link{plot.default}}}
  \item{log}{See \code{\link{plot.default}}}
  \item{main}{See \code{\link{plot.default}}}
  \item{sub}{See \code{\link{plot.default}}}
  \item{xlab}{A character string giving a label for the X axis.}
  \item{ylab}{A character string giving a label for the Y axis.}
  \item{ann}{See \code{\link{plot.default}}}
  \item{axes}{See \code{\link{plot.default}}}
  \item{frame.plot}{See \code{\link{plot.default}}}
  \item{panel.first}{See \code{\link{plot.default}}}
  \item{panel.last}{See \code{\link{plot.default}}}
  \item{asp}{See \code{\link{plot.default}}}
  \item{xgap.axis}{See \code{\link{plot.default}}}
  \item{ygap.axis}{See \code{\link{plot.default}}}
  \item{...}{Additional graphical parameters.}
}
\value{
Returns a plot.
}
\references{
Delory B.M., Li M., Topp C.N., Lobet G. (2018). archiDART v3.0: a new data analysis pipeline allowing the topological analysis of plant root systems, \emph{F1000Research}, 7:22, DOI: \doi{10.12688/f1000research.13541.1}

Li M., Duncan K., Topp C.N., Chitwood D.H. (2017) Persistent homology and the branching topologies of plants, \emph{American Journal of Botany}, DOI: 10.3732/ajb.1700046.
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\seealso{\code{\link{perhomology}}, \code{\link{par}}.
}
\examples{
## Locate folder with RSML file
path <- system.file("extdata", package="archiDART")

table<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age")

#Persistent homology
ph<-perhomology(table)

#Plot the barcode
plot(ph$'monocot-archisimple_1', las=1, xlab="Geodesic distance (cm)")
}