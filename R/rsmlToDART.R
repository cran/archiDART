rsmlToDART <- function(rsml.path, final.date, connect){
  
  # Create the plant object
  rsml <- xmlToList(xmlParse(rsml.path))
  plants <- rsml$scene
  
  # Create LIE and RAC files for each root system 
  
  n <- 0 # Initialise number of LIE/RAC files (as 1 RSML file can contain more than 1 first-order root)
  lie.all<-list()
  rac.all<-list()
  
  for (r0 in plants){# For each plant
    
    for (r1 in r0){ # For each first-order root
      
      r <- 0 # Initialise the number of root consisting a root system
      
      if (class(r1)=="list"){# For each first order root, create LIE and RAC files
        
        n<-n+1
        currentMother<-0
        lie<-data.frame(Num=0,Date=0, Bran=FALSE, Apic=FALSE, Prec=0, Suiv=0, X=0, Y=0, dist=0)
        rac<-data.frame(Root=0, Mother=0, Ord=0, DBase=0, DApp=0, Lengths=0)
        ns<-r1$geometry$polyline
        length1<-length(ns)
        r<-r+1
        cumulDist <- 0
        
        for (i in 1:length1){#Fill LIE file for the first order root
          
          lie[i,1]<-i #Num
          lie[i,2]<-1 #Date
          lie[i,3]<-FALSE #Bran
          lie[i,4]<-FALSE #Apic
          lie[i,7] <-  as.numeric(ns[[i]][[1]]) #X
          lie[i,8] <- as.numeric(ns[[i]][[2]]) #Y
          if (i==1){
            dist<-0
            lie[i,2]<-0 #Date
            lie[i,5]<-0 #Prec
            lie[i,6]<-i+1 #Suiv
            lie[i,3]<-TRUE #Bran
            start1<-lie$Num[i]} 
          else {
            if (i==length1){
              lie[i,5]<-i-1 #Prec
              lie[i,6]<-0 #Suiv
              lie[i,4]<-TRUE #Apic
              stop1<-lie$Num[i]}
            else{
              lie[i,5]<-i-1 #Prec
              lie[i,6]<-i+1} #Suiv
            dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]])
            dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]])
            dist <- sqrt(dx^2 + dy^2)}
          cumulDist<-cumulDist+dist
          lie[i,9]<-cumulDist} #dist
        
        # Fill RAC file for the first order root
        
        rac[r,1]<-0 #Root
        rac[r,2]<--1 #Mother
        rac[r,3]<-1 #Ord
        rac[r,4]<-0 #DBase
        rac[r,5]<-0 #DApp
        rac[r,6]<-cumulDist #Lengths
        
        #----------------------------------------------------------------------------------------------------
        
        # If there are lateral roots
        
        if ("root" %in% names(r1)){
          
          for (r2 in r1){# For each 2-order root
            
            if ("geometry" %in% names(r2)){
              
              r<-r+1
              ns <- r2$geometry$polyline
              length2<-length(ns)
              cumulDist <- 0
              lie.lines<-nrow(lie)
              
              for (i in 1:length2){# Fill LIE file for the 2-order root
                lie[lie.lines+i,1]<-lie.lines+i #Num
                lie[lie.lines+i,2]<-1 #Date
                lie[lie.lines+i,3]<-FALSE #Bran
                lie[lie.lines+i,4]<-FALSE #Apic
                lie[lie.lines+i,7]<-as.numeric(ns[[i]][[1]]) #X
                lie[lie.lines+i,8]<-as.numeric(ns[[i]][[2]]) #Y
                
                if (i==1){
                  dist<-0
                  # Search the closest point on the mother root (calculate DBase)
                  min <- 1e9
                  prec <- 0
                  dbase<-0
                  for(j in start1:stop1){
                    dx <- lie$X[lie.lines+i] - lie$X[j] # x
                    dy <- lie$Y[lie.lines+i] - lie$Y[j] # y
                    dist1 <- sqrt(dx^2 + dy^2)
                    if(dist1 < min ){
                      min <- dist1
                      dbase <- lie$dist[j]
                      prec<-lie$Num[j]}}
                  if (connect==FALSE) {lie[lie.lines+i,5]<-lie.lines+i} else {lie[lie.lines+i,5]<-prec} #Prec
                  lie[lie.lines+i,6]<-lie.lines+i+1 #Suiv
                  lie[lie.lines+i,3]<-TRUE #Bran
                  if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+i])^2+(lie$Y[prec]-lie$Y[lie.lines+i])^2)}
                  start2<-lie$Num[lie.lines+i]}
                else { 
                  if (i==length2){
                    lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                    lie[lie.lines+i,6]<-0 #Suiv
                    lie[lie.lines+i,4]<-TRUE #Apic
                    stop2<-lie$Num[lie.lines+i]}
                  else{
                    lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                    lie[lie.lines+i,6]<-lie.lines+i+1} #Suiv
                  dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]])
                  dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]])
                  dist <- sqrt(dx^2 + dy^2)}
                
                cumulDist<-cumulDist+dist
                lie[lie.lines+i,9]<-cumulDist} #dist
              
              # Fill RAC file for the 2-order root
              
              rac[r,1]<-max(rac$Root)+1 #Root
              rac[r,2]<-currentMother #Mother
              rac[r,3]<-2 #Ord
              rac[r,4]<-dbase #DBase
              rac[r,5]<-final.date/2 #DApp
              rac[r,6]<-cumulDist #Lengths
              
              currentMother2<-rac$Root[r]
              
              #----------------------------------------------------------------------------------------------------            
              
              if ("root" %in% names(r2)){
                
                for (r3 in r2){# For each 3-order root
                  
                  if ("geometry" %in% names(r3)){
                    
                    r<-r+1
                    ns <- r3$geometry$polyline
                    length3<-length(ns)
                    cumulDist <- 0
                    lie.lines<-nrow(lie)
                    
                    for (i in 1:length3){# Fill LIE file for the 3-order root
                      lie[lie.lines+i,1]<-lie.lines+i #Num
                      lie[lie.lines+i,2]<-1 #Date
                      lie[lie.lines+i,3]<-FALSE #Bran
                      lie[lie.lines+i,4]<-FALSE #Apic
                      lie[lie.lines+i,7]<-as.numeric(ns[[i]][[1]]) #X
                      lie[lie.lines+i,8]<-as.numeric(ns[[i]][[2]]) #Y
                      
                      if (i==1){
                        dist<-0
                        # Search the closest point on the mother root (calculate DBase)
                        min <- 1e9
                        prec <- 0
                        dbase<-0
                        for(j in start2:stop2){
                          dx <- lie$X[lie.lines+i] - lie$X[j] # x
                          dy <- lie$Y[lie.lines+i] - lie$Y[j] # y
                          dist1 <- sqrt(dx^2 + dy^2)
                          if(dist1 < min ){
                            min <- dist1
                            dbase <- lie$dist[j]
                            prec<-lie$Num[j]}}
                        if (connect==FALSE) {lie[lie.lines+i,5]<-lie.lines+i} else {lie[lie.lines+i,5]<-prec} #Prec
                        lie[lie.lines+i,6]<-lie.lines+i+1 #Suiv
                        lie[lie.lines+i,3]<-TRUE #Bran
                        if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+i])^2+(lie$Y[prec]-lie$Y[lie.lines+i])^2)}
                        start3<-lie$Num[lie.lines+i]}
                      else { 
                        if (i==length3){
                          lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                          lie[lie.lines+i,6]<-0 #Suiv
                          lie[lie.lines+i,4]<-TRUE #Apic
                          stop3<-lie$Num[lie.lines+i]}
                        else{
                          lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                          lie[lie.lines+i,6]<-lie.lines+i+1} #Suiv
                        dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]])
                        dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]])
                        dist <- sqrt(dx^2 + dy^2)}
                      
                      cumulDist<-cumulDist+dist
                      lie[lie.lines+i,9]<-cumulDist} #dist
                    
                    # Fill RAC file for the 3-order root
                    
                    rac[r,1]<-max(rac$Root)+1 #Root
                    rac[r,2]<-currentMother2 #Mother
                    rac[r,3]<-3 #Ord
                    rac[r,4]<-dbase #DBase
                    rac[r,5]<-final.date/2 #DAPP
                    rac[r,6]<-cumulDist #Lengths
                    
                    currentMother3<-rac$Root[r]
                    
                    #----------------------------------------------------------------------------------------------------
                    
                    if ("root" %in% names(r3)){
                      
                      for (r4 in r3){# For each 4-order root
                        
                        if ("geometry" %in% names(r4)){
                          
                          r<-r+1
                          ns <- r4$geometry$polyline
                          length4<-length(ns)
                          cumulDist <- 0
                          lie.lines<-nrow(lie)
                          
                          for (i in 1:length4){# Fill LIE file for the 4-order root
                            lie[lie.lines+i,1]<-lie.lines+i #Num
                            lie[lie.lines+i,2]<-1 #Date
                            lie[lie.lines+i,3]<-FALSE #Bran
                            lie[lie.lines+i,4]<-FALSE #Apic
                            lie[lie.lines+i,7]<-as.numeric(ns[[i]][[1]]) #X
                            lie[lie.lines+i,8]<-as.numeric(ns[[i]][[2]]) #Y
                            
                            if (i==1){
                              dist<-0
                              # Search the closest point on the mother root (calculate DBase)
                              min <- 1e9
                              prec <- 0
                              dbase<-0
                              for(j in start3:stop3){
                                dx <- lie$X[lie.lines+i] - lie$X[j] # x
                                dy <- lie$Y[lie.lines+i] - lie$Y[j] # y
                                dist1 <- sqrt(dx^2 + dy^2)
                                if(dist1 < min ){
                                  min <- dist1
                                  dbase <- lie$dist[j]
                                  prec<-lie$Num[j]}}
                              if (connect==FALSE) {lie[lie.lines+i,5]<-lie.lines+i} else {lie[lie.lines+i,5]<-prec} #Prec
                              lie[lie.lines+i,6]<-lie.lines+i+1 #Suiv
                              lie[lie.lines+i,3]<-TRUE #Bran
                              if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+i])^2+(lie$Y[prec]-lie$Y[lie.lines+i])^2)}
                              start4<-lie$Num[lie.lines+i]}
                            else { 
                              if (i==length4){
                                lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                                lie[lie.lines+i,6]<-0 #Suiv
                                lie[lie.lines+i,4]<-TRUE #Apic
                                stop4<-lie$Num[lie.lines+i]}
                              else{
                                lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                                lie[lie.lines+i,6]<-lie.lines+i+1} #Suiv
                              dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]])
                              dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]])
                              dist <- sqrt(dx^2 + dy^2)}
                            
                            cumulDist<-cumulDist+dist
                            lie[lie.lines+i,9]<-cumulDist} #dist
                          
                          # Fill RAC file for the 4-order root
                          
                          rac[r,1]<-max(rac$Root)+1 #Root
                          rac[r,2]<-currentMother3 #Mother
                          rac[r,3]<-4 #Ord
                          rac[r,4]<-dbase #DBase
                          rac[r,5]<-final.date/2 #DAPP
                          rac[r,6]<-cumulDist #Lengths
                          
                          currentMother4<-rac$Root[r]
                          
                          #----------------------------------------------------------------------------------------------------
                          
                          if ("root" %in% names(r4)){
                            
                            for (r5 in r4){# For each 5-order root
                              
                              if ("geometry" %in% names(r5)){
                                
                                r<-r+1
                                ns <- r5$geometry$polyline
                                length5<-length(ns)
                                cumulDist <- 0
                                lie.lines<-nrow(lie)
                                
                                for (i in 1:length5){# Fill LIE file for the 5-order root
                                  lie[lie.lines+i,1]<-lie.lines+i #Num
                                  lie[lie.lines+i,2]<-1 #Date
                                  lie[lie.lines+i,3]<-FALSE #Bran
                                  lie[lie.lines+i,4]<-FALSE #Apic
                                  lie[lie.lines+i,7]<-as.numeric(ns[[i]][[1]]) #X
                                  lie[lie.lines+i,8]<-as.numeric(ns[[i]][[2]]) #Y
                                  
                                  if (i==1){
                                    dist<-0
                                    # Search the closest point on the mother root (calculate DBase)
                                    min <- 1e9
                                    prec <- 0
                                    dbase<-0
                                    for(j in start4:stop4){
                                      dx <- lie$X[lie.lines+i] - lie$X[j] # x
                                      dy <- lie$Y[lie.lines+i] - lie$Y[j] # y
                                      dist1 <- sqrt(dx^2 + dy^2)
                                      if(dist1 < min ){
                                        min <- dist1
                                        dbase <- lie$dist[j]
                                        prec<-lie$Num[j]}}
                                    if (connect==FALSE) {lie[lie.lines+i,5]<-lie.lines+i} else {lie[lie.lines+i,5]<-prec} #Prec
                                    lie[lie.lines+i,6]<-lie.lines+i+1 #Suiv
                                    lie[lie.lines+i,3]<-TRUE #Bran
                                    if (connect==TRUE){cumulDist<-sqrt((lie$X[prec]-lie$X[lie.lines+i])^2+(lie$Y[prec]-lie$Y[lie.lines+i])^2)}}
                                  else { 
                                    if (i==length5){
                                      lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                                      lie[lie.lines+i,6]<-0 #Suiv
                                      lie[lie.lines+i,4]<-TRUE} #Apic
                                    else{
                                      lie[lie.lines+i,5]<-lie.lines+i-1 #Prec
                                      lie[lie.lines+i,6]<-lie.lines+i+1} #Suiv
                                    dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]])
                                    dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]])
                                    dist <- sqrt(dx^2 + dy^2)}
                                  
                                  cumulDist<-cumulDist+dist
                                  lie[lie.lines+i,9]<-cumulDist} #dist
                                
                                # Fill RAC file for the 5-order root
                                
                                rac[r,1]<-max(rac$Root)+1 #Root
                                rac[r,2]<-currentMother4 #Mother
                                rac[r,3]<-5 #Ord
                                rac[r,4]<-dbase #DBase
                                rac[r,5]<-final.date/2 #DAPP
                                rac[r,6]<-cumulDist #Lengths
                                
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
      lie$Bran[lie$Bran==TRUE]<-"true"
      lie$Bran[lie$Bran==FALSE]<-"false"
      lie$Apic[lie$Apic==TRUE]<-"true"
      lie$Apic[lie$Apic==FALSE]<-"false"
      lie.all[[n]]<-lie
      rac.all[[n]]<-rac
    }
  }
  result<-list(lie=lie.all, rac=rac.all)
  return(result)
}