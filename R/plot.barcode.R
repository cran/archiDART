plot.barcode<-function(x, y = NULL, type = "p",  xlim = NULL, ylim = NULL,
                       log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                       ann = par("ann"), axes = TRUE, frame.plot = axes,
                       panel.first = NULL, panel.last = NULL, asp = NA,
                       xgap.axis = NA, ygap.axis = NA,
                       ...){
  if (is.null(xlab)==TRUE){xlab<-"Geodesic distance"}
  if (is.null(ylab)==TRUE){ylab<-expression(H[0])}
  xcoord<-c(min(x[,2:3]), max(x[,2:3]))
  ycoord<-c(0, nrow(x))
  plot(x=xcoord, y=ycoord, type="n", xlab=xlab, ylab=ylab, ...)
  segments(x0=x[,"birth"], y0=c(1:nrow(x)), x1=x[,"death"], y1=c(1:nrow(x)), ...)}