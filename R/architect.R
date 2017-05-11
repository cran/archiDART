architect<-function(inputrac=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px", rsml.date=NULL, rsml.connect=FALSE){
  
  # Errors interception
  
  if (is.null(inputrac)==TRUE & is.null(inputtps)==TRUE & is.null(inputrsml)==TRUE){stop("inputrac/inputps and/or inputrsml must be provided")}
  
  if (is.null(inputrac)==FALSE) {if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}}
  
  if (is.null(inputtps)==FALSE) {if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}}
  
  if (is.null(inputrac)==FALSE|is.null(inputtps)==FALSE){
    if (is.null(inputrac)==TRUE|is.null(inputtps)==TRUE){stop("If inputrac/inputtps is not NULL, inputtps/inputrac must be provided")}}
  
  if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
  
  if (is.null(inputrac)==FALSE & is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0|length(res)>1){stop("res must be a single positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(rsml.date)==FALSE){
    if (is.character(rsml.date)==TRUE|is.numeric(rsml.date)==TRUE){} else {stop("If rsml.date is not NULL, rsml.date must be a character string or a positive numeric value")}
    if (is.numeric(rsml.date)==TRUE){if (rsml.date<=0|length(rsml.date)>1){stop("If mode(rsml.date) is numeric, rsml.date must be a single positive value")}}}
  
  if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
  
  # Reading of DART and rsml files
  
  if (is.null(inputtps)==FALSE){
    filenames.tps<-list.files(path=inputtps, pattern="\\.tps$")
    path.tps<-rep(inputtps, length.out=length(filenames.tps))
    filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
    message(paste("Number of DART tps files in inputtps:", length(filenames.tps), sep=" "))}
  
  if (is.null(inputrac)==FALSE){
    filenames.rac<-list.files(path=inputrac, pattern="\\.rac$")
    path.rac<-rep(inputrac, length.out=length(filenames.rac))
    filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
    message(paste("Number of DART rac files in inputrac:", length(filenames.rac), sep=" "))}
  
  if (is.null(inputrsml)==FALSE) {
    filenames.rsml<-list.files(path=inputrsml, pattern="\\.rsml$")
    path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
    filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
    message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
  
  if (is.null(inputrsml)==TRUE){
    if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
    if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}}
  else {
    if (is.null(inputrac)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
    else{
      if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
      if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
      if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}

  if (is.null(inputrsml)==TRUE){ # Only DART files
  
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      if (length(TIME)==1) {
        age<-list()
        obstot<-0
        for (i in 1:length(DATA)) {
          age[[i]]<-TIME[[1]]$Date
          obstot<-obstot+length(age[[i]])}
        for (i in 1:length(DATA)) {if(length(age[[i]])!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
      
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac files in inputrac and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac files and their corresponding tps files must have the same name")}}
          age<-list()
          obstot<-0
          for (i in 1:length(TIME)) {
            age[[i]]<-TIME[[i]]$Date
            obstot<-obstot+length(age[[i]])}
          for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}}
  
  else {
    if (is.null(inputrac)==TRUE){ # Only RSML files
      
      DATA<-list()
      age<-list()
      res1<-c()
      filenamesrac<-c()
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect) # Read RSML files
      obstot<-0
      for (i in 1:length(RSML)){
        res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution), length(RSML[[i]]$lie)))
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        
        l<-length(age)
        for (j in 1:length1) {
          age[[l+j]]<-RSML[[i]]$tps[[j]]$Date
          obstot<-obstot+length(age[[l+j]])}

        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}}
    
    else { # DART and RSML files
      
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE) # Read TPS files
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1) # Read RAC files
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      res1<-rep(res, length(filenames.rac))
      
      if (length(TIME)==1) {
        age<-list()
        obstot<-0
        for (i in 1:length(DATA)) {
          age[[i]]<-TIME[[1]]$Date
          obstot<-obstot+length(age[[i]])}
        for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac files in inputrac and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac files and their corresponding tps files must have the same name")}}
          age<-list()
          obstot<-0
          for (i in 1:length(DATA)) {
            age[[i]]<-TIME[[i]]$Date
            obstot<-obstot+length(age[[i]])}
          for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}
      
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect) # Read RSML files
      
      for (i in 1:length(RSML)){
        res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution), length(RSML[[i]]$lie)))
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        
        l<-length(age)
        for (j in 1:length1) {
          age[[l+j]]<-RSML[[i]]$tps[[j]]$Date
          obstot<-obstot+length(age[[l+j]])}
        
        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}}}
  
  # Creating vectors and matrices for root architecture parameters calculation
  
  FileNames<-c()
  Time<-c()
  FirstOrderRootLength<-c()
  FirstOrderRootNumber<-c()
  TotalRootLength<-c()
  TotalLateralRootNumber<-c()
  TotalLateralRootLength<-c()
  LateralRootDensity<-c()
  GrowthRateFirstOrderRoot<-c()
  GrowthRateTotalRoot<-c()
  
  # Calculation of root architecture parameters
  
  k<-0
  maxord<-max(sapply(DATA, function(x){max(x[[3]])}))
  if (maxord>1) {latroot<-matrix(ncol=4*(maxord-1), nrow=obstot)}
  
  for(i in 1:length(DATA)){
    
    LastFirstOrderRootLength<-DATA[[i]][[ncol(DATA[[i]])]][DATA[[i]]$Ord==1]
    
    for (t in 1:(ncol(DATA[[i]])-5)){
      
      k<-k+1
      
      # Time
      Time[k]<-age[[i]][t]
      
      #Split data per root order
      data.split<-split(DATA[[i]][[paste("Lengths",t,sep="")]],as.factor(DATA[[i]]$Ord))
      
      # File names
      
      FileNames[k]<-filenamesrac[i]
      
      # Unit conversion
      
      if (is.null(inputrsml)==FALSE){res<-res1[i]}
      
      if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
      if (unitlength=="cm") {cunit<-(cm(1)/res)}
      if (unitlength=="px") {cunit<-1}	
      
      # Total root length
      
      TotalRootLength[k]<-sum(DATA[[i]][[paste("Lengths",t,sep="")]])*cunit
      
      # Growth rate of the root system
      
      if (t==1) {GrowthRateTotalRoot[k]<-TotalRootLength[k]/age[[i]][t]}
      if (t>1) {GrowthRateTotalRoot[k]<-(TotalRootLength[k]-TotalRootLength[k-1])/(age[[i]][t]-age[[i]][t-1])}
      
      # First-order root length
      
      FirstOrderRootLength[k]<-sum(data.split$'1')*cunit
      
      # Growth rate of the first-order root
      
      if (t==1) {GrowthRateFirstOrderRoot[k]<-FirstOrderRootLength[k]/age[[i]][t]}
      if (t>1) {GrowthRateFirstOrderRoot[k]<-(FirstOrderRootLength[k]-FirstOrderRootLength[k-1])/(age[[i]][t]-age[[i]][t-1])}
      
      # Total number of first-order roots
      
      FirstOrderRootNumber[k]<-sum(data.split$'1'!=0)
      
      # Total number of lateral roots
      
      if (FirstOrderRootLength[k]==0) {TotalLateralRootNumber[k]<-0}  else {TotalLateralRootNumber[k]<-sum(DATA[[i]][[paste("Lengths",t,sep="")]]!=0)-FirstOrderRootNumber[k]}
      
      # Total length of lateral roots
      
      TotalLateralRootLength[k]<-TotalRootLength[k]-FirstOrderRootLength[k]
      
      # Number, length and growth rate of lateral roots (by branching order)		
      
      if (maxord>1){
      
      for (l in 1:(maxord-1)){
        if (l<=max(DATA[[i]]$Ord)-1) {latroot[k,l]<-sum(data.split[[l+1]]!=0)} else {latroot[k,l]<-0}
        if (l<=max(DATA[[i]]$Ord)-1) {latroot[k,(l+(maxord-1))]<-sum(data.split[[l+1]])*cunit} else {latroot[k,(l+(maxord-1))]<-0}
        if (latroot[k,l]==0){latroot[k,(l+2*(maxord-1))]<-0} else {latroot[k,(l+2*(maxord-1))]<-latroot[k,(l+(maxord-1))]/latroot[k,l]}
        if (t==1) {latroot[k,(l+3*(maxord-1))]<-latroot[k,(l+(maxord-1))]/age[[i]][t]} else {latroot[k,(l+3*(maxord-1))]<-(latroot[k,(l+(maxord-1))]-latroot[k-1,(l+(maxord-1))])/(age[[i]][t]-age[[i]][t-1])}}}
      
      # Density of secondary roots on the first-order root
      
      if (FirstOrderRootLength[k]==0|maxord==1) {LateralRootDensity[k]<-0} else {LateralRootDensity[k]<-latroot[k,1]/FirstOrderRootLength[k]}}}	
  
  # Summarying results in a data frame
  

  if (maxord>1) {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, FirstOrderRootNumber, TotalLateralRootNumber, TotalLateralRootLength, latroot, LateralRootDensity)}
  else {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, FirstOrderRootNumber, TotalLateralRootNumber, TotalLateralRootLength, LateralRootDensity)}
  
  if (maxord>1){
  LRnumberheading<-c()
  LRlengthheading<-c()
  LRmeanlengthheading<-c()
  LRgrowthrateheading<-c()
  
  for (h in 2:maxord){
    LRnumberheading[h-1]<-paste("N", h, "LR", sep="")
    LRlengthheading[h-1]<-paste("L", h, "LR", sep="")
    LRmeanlengthheading[h-1]<-paste("ML", h, "LR", sep="")
    LRgrowthrateheading[h-1]<-paste("GR", h, "L", sep="")}}
  
  if (maxord>1) {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TN1R", "TNLR", "TLRL", t(LRnumberheading), t(LRlengthheading), t(LRmeanlengthheading), t(LRgrowthrateheading), "D2LR")}
  else {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TN1R", "TNLR", "TLRL", "D2LR")}
  
  return(outputresults)}