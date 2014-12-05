latdist<-function(inputrac, output=c("lrd","dtp"), res=NULL, unitlength="px", int.length=NULL, interpol=NULL){
  
  # Errors interception
  
  if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}
  
  if (output=="lrd"|output=="dtp"){} else {stop("The character string in output must be lrd or dtp")}
  
  if (is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0){stop("res must be a positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(interpol)==FALSE) {
    if (mode(interpol)!="numeric"){stop("mode(interpol) must be numeric")}
    if (interpol<=0){stop("interpol must be a positive value")}}
  
  if (is.null(int.length)==FALSE) {
    if (mode(int.length)!="numeric"){stop("mode(int.length) must be numeric")}
    if (int.length<=0){stop("int.length must be a positive value")}}
  
  # Reading of DART output files
  
  filenames.rac<-list.files(path=inputrac, pattern="\\.rac$")
  path.rac<-rep(inputrac, length.out=length(filenames.rac))
  filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
  DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
  print(paste("Number of DART rac files in inputrac:", length(DATA), sep=" "))
  for (i in 1:length(DATA)) {
    colnames(DATA[[i]])<-c()
    colnames(DATA[[i]])[1]<-"Root"
    colnames(DATA[[i]])[2]<-"Mother"
    colnames(DATA[[i]])[3]<-"Ord"
    colnames(DATA[[i]])[4]<-"DBase"
    colnames(DATA[[i]])[5]<-"DApp"
    for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
 
  # Unit conversion
  
  if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
  if (unitlength=="cm") {cunit<-(cm(1)/res)}
  if (unitlength=="px") {cunit<-1}
  for (i in 1:length(DATA)){
    DATA[[i]]$DBase<-DATA[[i]]$DBase*cunit
    for (j in 6:ncol(DATA[[i]])) {DATA[[i]][,j]<-DATA[[i]][,j]*cunit}}
    
  # Calculating lateral root distributions

  rac<-list()
  res<-list()
  
  if (output=="lrd"){
  
    if (is.null(interpol)==TRUE) {
      
      # Calculating lateral root length and density distributions alons each mother root without interpolation
    
        for (i in 1:length(DATA)){
          
          gen<-list()
          LRdensity<-list()
          latrootnum<-c()
          
          for (j in 1:nrow(DATA[[i]])){
            Root<-NULL
            DBase<-NULL
            sub<-subset(DATA[[i]], Mother==j-1, select=c(Root, DBase))
            latrootnum[j]<-nrow(sub)
            if (nrow(sub)==0) {gen[[j]]<-NULL} else {gen[[j]]<-sub}}
          
          finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
          global.latrootdensity<-latrootnum/finalrootlength
          
          rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
          
          for (j in 1:length(gen)){
            
            if (is.null(gen[[j]])==FALSE){
              
              lrdensity<-c()
              dbase<-c()
              lrlength<-c()
              t<-0
            
              for (k in 1:nrow(gen[[j]])){
                
                if (gen[[j]]$DBase[k]>=(int.length/2) & gen[[j]]$DBase[k]<=finalrootlength[j]-(int.length/2)){
                  
                  t<-t+1
                
                  min<-gen[[j]]$DBase[k]-(int.length/2)
                  max<-gen[[j]]$DBase[k]+(int.length/2)
                  
                  dbase[t]<-gen[[j]]$DBase[k]
                  lrdensity[t]<-sum(gen[[j]]$DBase>=min & gen[[j]]$DBase<=max)/int.length
                  lrlength[t]<-sum(finalrootlength[1+gen[[j]]$Root[gen[[j]]$DBase>=min & gen[[j]]$DBase<=max]])/int.length}}
            
              if (length(dbase)!=0) {LRdensity[[j]]<-data.frame(DBase=dbase, LRD=lrdensity, LRL=lrlength)} else {LRdensity[[j]]<-NULL}}
          
            else {
              
              LRdensity[[j]]<-NULL}}
          
          res[[i]]<-LRdensity}
        
          names(rac)<-filenamesrac
          names(res)<-filenamesrac
          outputresults<-list(root=rac, res=res)}
  
    else {
      
      # Calculating lateral root length and density distributions alons each mother root with interpolation
    
      for (i in 1:length(DATA)){
        
        gen<-list()
        LRdensity<-list()
        latrootnum<-c()
        
        for (j in 1:nrow(DATA[[i]])){
          Mother<-NULL
          DBase<-NULL
          Root<-NULL
          sub<-subset(DATA[[i]], Mother==j-1, select=c(Root, DBase))
          latrootnum[j]<-nrow(sub)
          if (nrow(sub)==0) {gen[[j]]<-NULL} else {gen[[j]]<-sub}}
        
        finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
        global.latrootdensity<-latrootnum/finalrootlength
        
        rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
        
        for (j in 1:length(gen)){
          
          if (is.null(gen[[j]])==FALSE){
            
            lrdensity<-c()
            dbase<-c()
            lrlength<-c()
            seqDBase<-seq(from=0, to=finalrootlength[j], by=finalrootlength[j]/(interpol-1))
            t<-0
            
            for (k in 1:length(seqDBase)){
              
              if (seqDBase[k]>=(int.length/2) & seqDBase[k]<=finalrootlength[j]-(int.length/2)){
                
                t<-t+1
                
                min<-seqDBase[k]-(int.length/2)
                max<-seqDBase[k]+(int.length/2)
                
                dbase[t]<-seqDBase[k]
                lrdensity[t]<-sum(gen[[j]]$DBase>=min & gen[[j]]$DBase<=max)/int.length
                lrlength[t]<-sum(finalrootlength[1+gen[[j]]$Root[gen[[j]]$DBase>=min & gen[[j]]$DBase<=max]])/int.length}}
            
            if (length(dbase)!=0) {LRdensity[[j]]<-data.frame(DBase=dbase, LRD=lrdensity, LRL=lrlength)} else {LRdensity[[j]]<-NULL}}
          
          else {
            
            LRdensity[[j]]<-NULL}}
        
        res[[i]]<-LRdensity}
      
      names(rac)<-filenamesrac
      names(res)<-filenamesrac
      outputresults<-list(root=rac, res=res)}}

  
  if (output=="dtp"){
    
    # Calculating the distance between neighbouring roots along each mother root
    
    for (i in 1:length(DATA)){
      
      gen<-list()
      dtp<-list()
      latrootnum<-c()
      
      for (j in 1:nrow(DATA[[i]])){
        DBase<-NULL
        sub<-subset(DATA[[i]], Mother==j-1, select=DBase)
        latrootnum[j]<-nrow(sub)
        if (nrow(sub)<2) {gen[[j]]<-NULL} else {
          gen[[j]]<-sub[order(sub$DBase),]}}
      
      finalrootlength<-DATA[[i]][,ncol(DATA[[i]])]
      global.latrootdensity<-latrootnum/finalrootlength
      
      rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase, LatRootNum=latrootnum, FinalRootLength=finalrootlength, LatRootDensity=global.latrootdensity)
      
      for (j in 1:length(gen)){
        
        if (is.null(gen[[j]])==FALSE) {
        
          difftoprec<-diff(gen[[j]])
          dbase<-gen[[j]][2:length(gen[[j]])]
          
          dtp[[j]]<-data.frame(DBase=dbase, DTP=difftoprec)}
        
        else {
          
          dtp[[j]]<-NULL}}
      
      res[[i]]<-dtp}
    
    names(rac)<-filenamesrac
    names(res)<-filenamesrac
    outputresults<-list(root=rac, res=res)}

return(outputresults)}