rsmlToDART <- function(rsml.path, final.date, connect){
  
  # Create the plant object
  rsml <- xmlToList(xmlParse(rsml.path))
  plants <- rsml$scene
  resolution<-rsml$metadata$resolution
  
  # Find unittime in 'property-definition'
  
  unittime1<-NULL
  
  if (is.character(final.date)==TRUE) {
    property<-rsml$metadata$'property-definition'
    
    if (is.null(property)==FALSE){
      for (i in 1:length(property)){if (property[[i]]$label==final.date){unittime1<-as.character(property[[i]]$unit)}}
      if (is.null(unittime1)==TRUE){stop(paste("No", final.date, "label in rsml metadata (property-definition)", sep=" "))}}}
  
  # Create LIE and RAC files for each root system 
  
  n <- 0 # Initialise number of LIE/RAC files (as 1 RSML file can contain more than 1 root system)
  lie.all<-list()
  rac.all<-list()
  tps.all<-list()
  
  for (r0 in plants){# For each plant
    
    n<-n+1 #Add one unit for each root system
    r <- 0 # Initialise the number of roots consisting a root system
    timeserie<-FALSE #Does the rsml file contain time series data? (Ex: the root system age for each node)
    
    #For each plant root system, create LIE, RAC and TPS files
    lie<-data.frame(Num=0,Date=0, Bran=0, Apic=0, Prec=0, Suiv=0, X=0, Y=0, Z=0, dist=0, root=0)
    rac<-data.frame(Root=0, Mother=0, Ord=0, DBase=0, DApp=0, Lengths1=0)
    
    for (r1 in r0){ # For each first-order root
      
      if (class(r1)=="list"){
        
        ns<-r1$geometry$polyline
        
        # Extract plant age when nodes appeared for LIE and TPS files
        
        age=NULL
        
        if (is.character(final.date)==TRUE & "functions" %in% names(r1)){
          age<-r1$functions
          for (i in 1:length(age)){if (final.date %in% age[[i]]$.attrs) {
            timeserie<-TRUE
            time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
        
        length1<-length(ns)
        r<-r+1 #Add one unit for each root
        
        if (r==1){#For the first first-order root only
        
        currentMother<-r
          
        #c(Num, Date, Bran, Apic, Prec, Suiv)
        if (length(age)==0) {lie[1:length1,1:6]<-c(c(1:length1),rep(1, length1),rep(0, length1), rep(0, length1), 0:(length1-1), 2:(length1+1))}
        if (length(age)!=0) {lie[1:length1,1:6]<-c(c(1:length1),time,rep(0, length1), rep(0, length1), 0:(length1-1), 2:(length1+1))}
        lie[1:length1,11]<-rep(r, length1)
        
        #c(X,Y,Z)
        lie[1:length1,7]<-sapply(ns, xnodes)
        lie[1:length1,8]<-sapply(ns, ynodes)
        if (length(ns[[1]])==3) {lie[1:length1,9]<-sapply(ns, znodes)} else {lie[1:length1,9]<-0}
 
        #c(dist)
        lie[1:length1, 10]<-c(0, cumsum(sqrt((diff(lie$X[1:length1]))^2+(diff(lie$Y[1:length1]))^2+(diff(lie$Z[1:length1]))^2)))
        
        #For the first point of the first order root
        if (length(age)==0) {lie[1,c(2:3)]<-c(0,1)} else {lie[1,3]<-1}
        start1<-as.numeric(lie$Num[1])
        
        #For the last point of the first order root
        lie[length1,c(4,6)]<-c(1,0)
        stop1<-as.numeric(lie$Num[length1])
        
        # Fill RAC file for the first order root
        #c(Root, Mother, Ord, DBase, DApp, Length)
        cumulDist<-sum(sqrt((diff(lie$X[1:length1]))^2+(diff(lie$Y[1:length1]))^2+(diff(lie$Z[1:length1]))^2))
        rac[r,1:6]<-c(r,-1,1,0,0,cumulDist)}
        
        else {#For the other first-order roots only
          
          currentMother<-r
          
          lie.lines<-nrow(lie)
          
          #c(Num, Date, Bran, Apic, Prec, Suiv)
          if (length(age)==0) {lie[(lie.lines+1):(lie.lines+length1),1:6]<-c(c((lie.lines+1):(lie.lines+length1)),rep(1, length1),rep(0, length1), rep(0, length1), c(0,(lie.lines+1):(lie.lines+length1-1)), (lie.lines+2):(lie.lines+length1+1))}
          if (length(age)!=0) {lie[(lie.lines+1):(lie.lines+length1),1:6]<-c(c((lie.lines+1):(lie.lines+length1)),time,rep(0, length1), rep(0, length1), c(0,(lie.lines+1):(lie.lines+length1-1)), (lie.lines+2):(lie.lines+length1+1))}
          lie[(lie.lines+1):(lie.lines+length1),11]<-rep(r, length1)
          
          #c(X,Y,Z)
          lie[(lie.lines+1):(lie.lines+length1),7]<-sapply(ns, xnodes)
          lie[(lie.lines+1):(lie.lines+length1),8]<-sapply(ns, ynodes)
          if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length1),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length1),9]<-0}

          #c(dist)
          lie[(lie.lines+1):(lie.lines+length1), 10]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length1)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length1)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length1)]))^2)))
          
          #For the first point of the first order root
          if (length(age)==0) {lie[lie.lines+1,c(2:3)]<-c(0,1)} else {lie[lie.lines+1,3]<-1}
          start1<-as.numeric(lie$Num[lie.lines+1])
          
          #For the last point of the first order root
          lie[(lie.lines+length1),c(4,6)]<-c(1,0)
          stop1<-as.numeric(lie$Num[lie.lines+length1])
          
          # Fill RAC file for the first order root
          #c(Root, Mother, Ord, DBase, DApp, Length)
          cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length1)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length1)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length1)]))^2))
          rac[r,1:6]<-c(r,-1,1,0,0,cumulDist)}
        
        #----------------------------------------------------------------------------------------------------
        
        # If there are lateral roots
        
        if ("root" %in% names(r1)){
          
          for (r2 in r1){# For each 2-order root
            
            if ("geometry" %in% names(r2)){
              
              r<-r+1
              ns <- r2$geometry$polyline
              length2<-length(ns)
              lie.lines<-nrow(lie)
              
              age=NULL
              if (is.character(final.date)==TRUE & "functions" %in% names(r2)){
                age<-r2$functions
                for (i in 1:length(age)){if (final.date %in% age[[i]]$.attrs) {time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
              
              #c(Num, Date, Bran, Apic, Prec, Suiv)
              if (length(age)==0) {lie[(lie.lines+1):(lie.lines+length2),1:6]<-c((lie.lines+1):(lie.lines+length2), rep(1, length2), rep(0, length2), rep(0, length2), lie.lines:(lie.lines+length2-1), (lie.lines+2):(lie.lines+length2+1))}
              if (length(age)!=0) {lie[(lie.lines+1):(lie.lines+length2),1:6]<-c((lie.lines+1):(lie.lines+length2), time, rep(0, length2), rep(0, length2), lie.lines:(lie.lines+length2-1), (lie.lines+2):(lie.lines+length2+1))}
              lie[(lie.lines+1):(lie.lines+length2),11]<-rep(r, length2)
              
              #c(X,Y,Z)
              lie[(lie.lines+1):(lie.lines+length2),7]<-sapply(ns, xnodes)
              lie[(lie.lines+1):(lie.lines+length2),8]<-sapply(ns, ynodes)
              if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length2),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length2),9]<-0}
              
              #c(dist)
              lie[(lie.lines+1):(lie.lines+length2), 10]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length2)]))^2)))
              
              # Search the closest point on the mother root (calculate DBase)
              scalx<-diff(lie$X[start1:stop1])*(lie$X[start1:(stop1-1)]-lie$X[lie.lines+1])
              scaly<-diff(lie$Y[start1:stop1])*(lie$Y[start1:(stop1-1)]-lie$Y[lie.lines+1])
              scalz<-diff(lie$Z[start1:stop1])*(lie$Z[start1:(stop1-1)]-lie$Z[lie.lines+1])
              d2<-diff(lie$X[start1:stop1])^2+diff(lie$Y[start1:stop1])^2+diff(lie$Z[start1:stop1])^2
              t<-(-(scalx+scaly+scalz)/d2)

              if (sum(t>=0 & t<=1)==0){
                
                index<-min(which(t<0))
                dbase<-lie$dist[start1+index-1]
                if (connect==TRUE){
                  lie[lie.lines+1,5]<-lie$Num[start1+index-1]
                  dist1<-sqrt((lie$X[start1+index-1]-lie$X[lie.lines+1])^2+(lie$Y[start1+index-1]-lie$Y[lie.lines+1])^2+(lie$Z[start1+index-1]-lie$Z[lie.lines+1])^2)}}
              
              else {
              
              t[t<0]<-NA
              t[t>1]<-NA
              xn<-diff(lie$X[start1:stop1])*t+lie$X[start1:(stop1-1)]
              yn<-diff(lie$Y[start1:stop1])*t+lie$Y[start1:(stop1-1)]
              zn<-diff(lie$Z[start1:stop1])*t+lie$Z[start1:(stop1-1)]
              dist1<-sqrt((xn-lie$X[lie.lines+1])^2+(yn-lie$Y[lie.lines+1])^2+(zn-lie$Z[lie.lines+1])^2)
              
              if (sum(is.na(dist1)==T)>0) {
                index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                dist1<-min(dist1, na.rm=T)}
              else {
                index<-as.numeric(match(min(dist1), dist1))
                dist1<-min(dist1)}
              
              dbase<-lie$dist[start1+index-1]+distance3D(x1=lie$X[start1+index-1], x2=xn[index], y1=lie$Y[start1+index-1], y2=yn[index], z1=lie$Z[start1+index-1], z2=zn[index])

              if (connect==TRUE){
              lie[nrow(lie)+1,1:11]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie$root[start1+index-1])
              lie<-lie[order(lie$root, lie$dist),]
              length1<-length1+1
              lie.lines<-lie.lines+1
              stop1<-stop1+1
              pos<-match(NA, lie$Num)
              lie$Prec<-match(lie$Prec, lie$Num)
              lie$Suiv<-match(lie$Suiv, lie$Num)
              lie$Prec[pos]<-pos-1
              lie$Suiv[pos]<-pos+1
              lie$Date[pos]<-lie$Date[lie$Suiv[pos]]
              lie$Prec[pos+1]<-pos
              lie$Suiv[pos-1]<-pos
              lie$Prec[is.na(lie$Prec)==TRUE]<-0
              lie$Suiv[is.na(lie$Suiv)==TRUE]<-0
              lie$Num<-match(lie$Num, lie$Num)
              lie[lie.lines+1,5]<-pos}}
              
              if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
              
              lie[lie.lines+1,3]<-1
              start2<-as.numeric(lie$Num[lie.lines+1])
              
              # Change Suiv and Apic values for the last point of a lateral root
              lie[lie.lines+length2,c(4,6)]<-c(1, 0)
              stop2<-as.numeric(lie$Num[lie.lines+length2])
              
              # Fill RAC file for the 2-order root
              
              if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length2)]))^2))}
              else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length2)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length2)]))^2))}
              
              rac[r, 1:6]<-c(max(rac$Root)+1, currentMother, 2, dbase, 1, cumulDist)
              
              currentMother2<-rac$Root[r]
              
              #----------------------------------------------------------------------------------------------------            
              
              if ("root" %in% names(r2)){
                
                for (r3 in r2){# For each 3-order root
                  
                  if ("geometry" %in% names(r3)){
                    
                    r<-r+1
                    ns <- r3$geometry$polyline
                    length3<-length(ns)
                    lie.lines<-nrow(lie)
                    
                    age=NULL
                    if (is.character(final.date)==TRUE & "functions" %in% names(r3)){
                      age<-r3$functions
                      for (i in 1:length(age)){if (final.date %in% age[[i]]$.attrs) {time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                    
                    #c(Num, Date, Bran, Apic, Prec, Suiv)
                    if (length(age)==0) {lie[(lie.lines+1):(lie.lines+length3), 1:6]<-c((lie.lines+1):(lie.lines+length3), rep(1, length3), rep(0, length3), rep(0, length3), lie.lines:(lie.lines+length3-1), (lie.lines+2):(lie.lines+length3+1))}
                    if (length(age)!=0) {lie[(lie.lines+1):(lie.lines+length3), 1:6]<-c((lie.lines+1):(lie.lines+length3), time, rep(0, length3), rep(0, length3), lie.lines:(lie.lines+length3-1), (lie.lines+2):(lie.lines+length3+1))}
                    lie[(lie.lines+1):(lie.lines+length3),11]<-rep(r, length3)
                    
                    #c(X,Y,Z)
                    lie[(lie.lines+1):(lie.lines+length3),7]<-sapply(ns, xnodes)
                    lie[(lie.lines+1):(lie.lines+length3),8]<-sapply(ns, ynodes)
                    if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length3),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length3),9]<-0}
                    
                    #c(dist)
                    lie[(lie.lines+1):(lie.lines+length3), 10]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length3)]))^2)))
                    
                    # Search the closest point on the mother root (calculate DBase)
                    scalx<-diff(lie$X[start2:stop2])*(lie$X[start2:(stop2-1)]-lie$X[lie.lines+1])
                    scaly<-diff(lie$Y[start2:stop2])*(lie$Y[start2:(stop2-1)]-lie$Y[lie.lines+1])
                    scalz<-diff(lie$Z[start2:stop2])*(lie$Z[start2:(stop2-1)]-lie$Z[lie.lines+1])
                    d2<-diff(lie$X[start2:stop2])^2+diff(lie$Y[start2:stop2])^2+diff(lie$Z[start2:stop2])^2
                    t<-(-(scalx+scaly+scalz)/d2)
                    
                    if (sum(t>=0 & t<=1)==0){
                      
                      index<-min(which(t<0))
                      dbase<-lie$dist[start2+index-1]
                      if (connect==TRUE){
                        lie[lie.lines+1,5]<-lie$Num[start2+index-1]
                        dist1<-sqrt((lie$X[start2+index-1]-lie$X[lie.lines+1])^2+(lie$Y[start2+index-1]-lie$Y[lie.lines+1])^2+(lie$Z[start2+index-1]-lie$Z[lie.lines+1])^2)}}
                    
                    else {
                      
                      t[t<0]<-NA
                      t[t>1]<-NA
                      xn<-diff(lie$X[start2:stop2])*t+lie$X[start2:(stop2-1)]
                      yn<-diff(lie$Y[start2:stop2])*t+lie$Y[start2:(stop2-1)]
                      zn<-diff(lie$Z[start2:stop2])*t+lie$Z[start2:(stop2-1)]
                      dist1<-sqrt((xn-lie$X[lie.lines+1])^2+(yn-lie$Y[lie.lines+1])^2+(zn-lie$Z[lie.lines+1])^2)
                      
                      if (sum(is.na(dist1)==T)>0) {
                        index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                        dist1<-min(dist1, na.rm=T)}
                      else {
                        index<-as.numeric(match(min(dist1), dist1))
                        dist1<-min(dist1)}
                      
                      dbase<-lie$dist[start2+index-1]+distance3D(x1=lie$X[start2+index-1], x2=xn[index], y1=lie$Y[start2+index-1], y2=yn[index], z1=lie$Z[start2+index-1], z2=zn[index])
                      
                      if (connect==TRUE){
                        lie[nrow(lie)+1,1:11]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie$root[start2+index-1])
                        lie<-lie[order(lie$root, lie$dist),]
                        length1<-length1+1
                        lie.lines<-lie.lines+1
                        stop2<-stop2+1
                        pos<-match(NA, lie$Num)
                        lie$Prec<-match(lie$Prec, lie$Num)
                        lie$Suiv<-match(lie$Suiv, lie$Num)
                        lie$Prec[pos]<-pos-1
                        lie$Suiv[pos]<-pos+1
                        lie$Date[pos]<-lie$Date[lie$Suiv[pos]]
                        lie$Prec[pos+1]<-pos
                        lie$Suiv[pos-1]<-pos
                        lie$Prec[is.na(lie$Prec)==TRUE]<-0
                        lie$Suiv[is.na(lie$Suiv)==TRUE]<-0
                        lie$Num<-match(lie$Num, lie$Num)
                        lie[lie.lines+1,5]<-pos}}
                    
                    if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                    
                    lie[lie.lines+1,3]<-1
                    start3<-as.numeric(lie$Num[lie.lines+1])
                    
                    # Change Suiv and Apic values for the last point of a lateral root
                    lie[lie.lines+length3,c(4,6)]<-c(1, 0)
                    stop3<-as.numeric(lie$Num[lie.lines+length3])
                    
                    # Fill RAC file for the 3-order root
                    
                    if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length3)]))^2))}
                    else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length3)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length3)]))^2))}
                    rac[r, 1:6]<-c(max(rac$Root)+1, currentMother2, 3, dbase, 1, cumulDist)
                    
                    currentMother3<-rac$Root[r]
                    
                    #----------------------------------------------------------------------------------------------------
                    
                    if ("root" %in% names(r3)){
                      
                      for (r4 in r3){# For each 4-order root
                        
                        if ("geometry" %in% names(r4)){
                          
                          r<-r+1
                          ns <- r4$geometry$polyline
                          length4<-length(ns)
                          lie.lines<-nrow(lie)
                          
                          age=NULL
                          if (is.character(final.date)==TRUE & "functions" %in% names(r4)){
                            age<-r4$functions
                            for (i in 1:length(age)){if (final.date %in% age[[i]]$.attrs) {time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                          
                          #c(Num, Date, Bran, Apic, Prec, Suiv)
                          if (length(age)==0) {lie[(lie.lines+1):(lie.lines+length4), 1:6]<-c((lie.lines+1):(lie.lines+length4), rep(1, length4), rep(0, length4), rep(0, length4), lie.lines:(lie.lines+length4-1), (lie.lines+2):(lie.lines+length4+1))}
                          if (length(age)!=0) {lie[(lie.lines+1):(lie.lines+length4), 1:6]<-c((lie.lines+1):(lie.lines+length4), time, rep(0, length4), rep(0, length4), lie.lines:(lie.lines+length4-1), (lie.lines+2):(lie.lines+length4+1))}
                          lie[(lie.lines+1):(lie.lines+length4),11]<-rep(r, length4)

                          #c(X,Y,Z)
                          lie[(lie.lines+1):(lie.lines+length4),7]<-sapply(ns, xnodes)
                          lie[(lie.lines+1):(lie.lines+length4),8]<-sapply(ns, ynodes)
                          if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length4),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length4),9]<-0}
                          
                          #c(dist)
                          lie[(lie.lines+1):(lie.lines+length4), 10]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length4)]))^2)))
                          
                          # Search the closest point on the mother root (calculate DBase)
                          scalx<-diff(lie$X[start3:stop3])*(lie$X[start3:(stop3-1)]-lie$X[lie.lines+1])
                          scaly<-diff(lie$Y[start3:stop3])*(lie$Y[start3:(stop3-1)]-lie$Y[lie.lines+1])
                          scalz<-diff(lie$Z[start3:stop3])*(lie$Z[start3:(stop3-1)]-lie$Z[lie.lines+1])
                          d2<-diff(lie$X[start3:stop3])^2+diff(lie$Y[start3:stop3])^2+diff(lie$Z[start3:stop3])^2
                          t<-(-(scalx+scaly+scalz)/d2)
                          
                          if (sum(t>=0 & t<=1)==0){
                            
                            index<-min(which(t<0))
                            dbase<-lie$dist[start3+index-1]
                            if (connect==TRUE){
                              lie[lie.lines+1,5]<-lie$Num[start3+index-1]
                              dist1<-sqrt((lie$X[start3+index-1]-lie$X[lie.lines+1])^2+(lie$Y[start3+index-1]-lie$Y[lie.lines+1])^2+(lie$Z[start3+index-1]-lie$Z[lie.lines+1])^2)}}
                          
                          else {
                            
                            t[t<0]<-NA
                            t[t>1]<-NA
                            xn<-diff(lie$X[start3:stop3])*t+lie$X[start3:(stop3-1)]
                            yn<-diff(lie$Y[start3:stop3])*t+lie$Y[start3:(stop3-1)]
                            zn<-diff(lie$Z[start3:stop3])*t+lie$Z[start3:(stop3-1)]
                            dist1<-sqrt((xn-lie$X[lie.lines+1])^2+(yn-lie$Y[lie.lines+1])^2+(zn-lie$Z[lie.lines+1])^2)
                            
                            if (sum(is.na(dist1)==T)>0) {
                              index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                              dist1<-min(dist1, na.rm=T)}
                            else {
                              index<-as.numeric(match(min(dist1), dist1))
                              dist1<-min(dist1)}
                            
                            dbase<-lie$dist[start3+index-1]+distance3D(x1=lie$X[start3+index-1], x2=xn[index], y1=lie$Y[start3+index-1], y2=yn[index], z1=lie$Z[start3+index-1], z2=zn[index])
                            
                            if (connect==TRUE){
                              lie[nrow(lie)+1,1:11]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie$root[start3+index-1])
                              lie<-lie[order(lie$root, lie$dist),]
                              length1<-length1+1
                              lie.lines<-lie.lines+1
                              stop3<-stop3+1
                              pos<-match(NA, lie$Num)
                              lie$Prec<-match(lie$Prec, lie$Num)
                              lie$Suiv<-match(lie$Suiv, lie$Num)
                              lie$Prec[pos]<-pos-1
                              lie$Suiv[pos]<-pos+1
                              lie$Date[pos]<-lie$Date[lie$Suiv[pos]]
                              lie$Prec[pos+1]<-pos
                              lie$Suiv[pos-1]<-pos
                              lie$Prec[is.na(lie$Prec)==TRUE]<-0
                              lie$Suiv[is.na(lie$Suiv)==TRUE]<-0
                              lie$Num<-match(lie$Num, lie$Num)
                              lie[lie.lines+1,5]<-pos}}
                          
                          if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                          
                          lie[lie.lines+1,3]<-1
                          start4<-as.numeric(lie$Num[lie.lines+1])
                          
                          # Change Suiv and Apic values for the last point of a lateral root
                          lie[lie.lines+length4,c(4,6)]<-c(1, 0)
                          stop4<-as.numeric(lie$Num[lie.lines+length4])
                          
                          # Fill RAC file for the 4-order root
                          
                          if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length4)]))^2))}
                          else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length4)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length4)]))^2))}
                          rac[r, 1:6]<-c(max(rac$Root)+1, currentMother3, 4, dbase, 1, cumulDist)
                          
                          currentMother4<-rac$Root[r]
                          
                          #----------------------------------------------------------------------------------------------------
                          
                          if ("root" %in% names(r4)){
                            
                            for (r5 in r4){# For each 5-order root
                              
                              if ("geometry" %in% names(r5)){
                                
                                r<-r+1
                                ns <- r5$geometry$polyline
                                length5<-length(ns)
                                lie.lines<-nrow(lie)
                                
                                age=NULL
                                if (is.character(final.date)==TRUE & "functions" %in% names(r5)){
                                  age<-r5$functions
                                  for (i in 1:length(age)){if (final.date %in% age[[i]]$.attrs) {time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                
                                #c(Num, Date, Bran, Apic, Prec, Suiv)
                                if (length(age)==0) {lie[(lie.lines+1):(lie.lines+length5), 1:6]<-c((lie.lines+1):(lie.lines+length5), rep(1, length5), rep(0, length5), rep(0, length5), lie.lines:(lie.lines+length5-1), (lie.lines+2):(lie.lines+length5+1))}
                                if (length(age)!=0) {lie[(lie.lines+1):(lie.lines+length5), 1:6]<-c((lie.lines+1):(lie.lines+length5), rep(1, length5), rep(0, length5), rep(0, length5), lie.lines:(lie.lines+length5-1), (lie.lines+2):(lie.lines+length5+1))}
                                lie[(lie.lines+1):(lie.lines+length5),11]<-rep(r, length5)

                                #c(X,Y,Z)
                                lie[(lie.lines+1):(lie.lines+length5),7]<-sapply(ns, xnodes)
                                lie[(lie.lines+1):(lie.lines+length5),8]<-sapply(ns, ynodes)
                                if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length5),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length5),9]<-0}
                                
                                #c(dist)
                                lie[(lie.lines+1):(lie.lines+length5), 10]<-c(0, cumsum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length5)]))^2)))
                                
                                # Search the closest point on the mother root (calculate DBase)
                                scalx<-diff(lie$X[start4:stop4])*(lie$X[start4:(stop4-1)]-lie$X[lie.lines+1])
                                scaly<-diff(lie$Y[start4:stop4])*(lie$Y[start4:(stop4-1)]-lie$Y[lie.lines+1])
                                scalz<-diff(lie$Z[start4:stop4])*(lie$Z[start4:(stop4-1)]-lie$Z[lie.lines+1])
                                d2<-diff(lie$X[start4:stop4])^2+diff(lie$Y[start4:stop4])^2+diff(lie$Z[start4:stop4])^2
                                t<-(-(scalx+scaly+scalz)/d2)
                                
                                if (sum(t>=0 & t<=1)==0){
                                  
                                  index<-min(which(t<0))
                                  dbase<-lie$dist[start4+index-1]
                                  if (connect==TRUE){
                                    lie[lie.lines+1,5]<-lie$Num[start4+index-1]
                                    dist1<-sqrt((lie$X[start4+index-1]-lie$X[lie.lines+1])^2+(lie$Y[start4+index-1]-lie$Y[lie.lines+1])^2+(lie$Z[start4+index-1]-lie$Z[lie.lines+1])^2)}}
                                
                                else {
                                  
                                  t[t<0]<-NA
                                  t[t>1]<-NA
                                  xn<-diff(lie$X[start4:stop4])*t+lie$X[start4:(stop4-1)]
                                  yn<-diff(lie$Y[start4:stop4])*t+lie$Y[start4:(stop4-1)]
                                  zn<-diff(lie$Z[start4:stop4])*t+lie$Z[start4:(stop4-1)]
                                  dist1<-sqrt((xn-lie$X[lie.lines+1])^2+(yn-lie$Y[lie.lines+1])^2+(zn-lie$Z[lie.lines+1])^2)
                                  
                                  if (sum(is.na(dist1)==T)>0) {
                                    index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                    dist1<-min(dist1, na.rm=T)}
                                  else {
                                    index<-as.numeric(match(min(dist1), dist1))
                                    dist1<-min(dist1)}
                                  
                                  dbase<-lie$dist[start4+index-1]+distance3D(x1=lie$X[start4+index-1], x2=xn[index], y1=lie$Y[start4+index-1], y2=yn[index], z1=lie$Z[start4+index-1], z2=zn[index])
                                  
                                  if (connect==TRUE){
                                    lie[nrow(lie)+1,1:11]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie$root[start4+index-1])
                                    lie<-lie[order(lie$root, lie$dist),]
                                    length1<-length1+1
                                    lie.lines<-lie.lines+1
                                    stop4<-stop4+1
                                    pos<-match(NA, lie$Num)
                                    lie$Prec<-match(lie$Prec, lie$Num)
                                    lie$Suiv<-match(lie$Suiv, lie$Num)
                                    lie$Prec[pos]<-pos-1
                                    lie$Suiv[pos]<-pos+1
                                    lie$Date[pos]<-lie$Date[lie$Suiv[pos]]
                                    lie$Prec[pos+1]<-pos
                                    lie$Suiv[pos-1]<-pos
                                    lie$Prec[is.na(lie$Prec)==TRUE]<-0
                                    lie$Suiv[is.na(lie$Suiv)==TRUE]<-0
                                    lie$Num<-match(lie$Num, lie$Num)
                                    lie[lie.lines+1,5]<-pos}}
                                
                                if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                
                                lie[lie.lines+1,3]<-1
                                
                                # Change Suiv and Apic values for the last point of a lateral root
                                lie[lie.lines+length5,c(4,6)]<-c(1, 0)
                                
                                # Fill RAC file for the 5-order root
                                if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length5)]))^2))}
                                else {cumulDist<-sum(sqrt((diff(lie$X[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Y[(lie.lines+1):(lie.lines+length5)]))^2+(diff(lie$Z[(lie.lines+1):(lie.lines+length5)]))^2))}
                                rac[r, 1:6]<-c(max(rac$Root)+1, currentMother4, 4, dbase, 1, cumulDist)
                                
                              } 
                            }
                          } 
                        } 
                      }
                    } 
                  } 
                }
              } 
            } 
          }
        }
      }
    }
    
    lie$Bran[lie$Bran==0]<-"false"
    lie$Bran[lie$Bran==1]<-"true"
    lie$Apic[lie$Apic==0]<-"false"
    lie$Apic[lie$Apic==1]<-"true"
    if (sum(lie$Z)==0){lie<-lie[,-9]}
    
    #Create TPS file
    dates<-as.numeric(levels(as.factor(lie$Date)))
    
    if (is.null(unittime1)==TRUE){unittime1<-"unittime"}
    
    if (timeserie==FALSE) {
      if (is.null(final.date)==TRUE) {tps<-data.frame(Num=1, Date=1)}
      if (is.null(final.date)==FALSE) {tps<-data.frame(Num=1, Date=final.date)}}
    
    if (timeserie==TRUE) {tps<-data.frame(Num=c(1:length(dates)), Date=dates)}
    
    #Replace Date by Num in LIE file
    date.lie<-lie$Date
    
    for (i in 1:nrow(tps)){
      pos<-which(date.lie %in% tps$Date[i])
      num<-tps$Num[i]
      lie$Date[pos]<-num}
    
    #Make RAC file for time series
    
    rac$Root<-rac$Root-1
    rac$Mother[rac$Mother!=-1]<-rac$Mother[rac$Mother!=-1]-1
    
    if (timeserie==TRUE){
      
      rac<-rac[,-6] #Delete Lengths column because must be replaced by a matrix
      cols<-nrow(tps)
      rows<-nrow(rac)
      length1<-matrix(0, nrow=rows, ncol=cols)
      
      for (i in 1:cols){
        
        for (j in 1:rows){
        
          if ((tps$Num[i] %in% lie$Date[lie$root==j])==TRUE) {length1[j,i]<-max(lie$dist[lie$root==j & lie$Date==tps$Num[i]])}
          if ((tps$Num[i] %in% lie$Date[lie$root==j])==FALSE) {if (i>1) {length1[j,i]<-length1[j,i-1]}}}}
      
      colnames(length1)<-paste(rep("Lengths", cols), c(1:cols), sep="")
      
      rac<-cbind(rac, as.data.frame(length1))}
    
    lie<-lie[,c(-10,-11)] #Remove dist and root columns
    
    #Export RAC, LIE, TPS
    lie.all[[n]]<-lie
    rac.all[[n]]<-rac
    tps.all[[n]]<-tps
    
  }
  result<-list(resolution=resolution, time=unittime1, lie=lie.all, rac=rac.all, tps=tps.all)
  return(result)
}
