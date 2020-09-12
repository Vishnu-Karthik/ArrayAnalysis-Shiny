package_list <- c("affy",
                  "affycomp",
                  "affyPLM",
                  "bioDist",
                  "simpleaffy",
                  "affyQCReport",
                  "plier",
                  "yaqcaffy",
                  "gdata",
                  "gplots",
                  "shiny",
                  "shinythemes",
                  "shinyFiles",
                  "biomaRt",
                  "limma",
                  "clusterProfiler",
                  "org.Hs.eg.db",
                  "ArrayTools",
                  "makecdfenv")
lapply(package_list, require, character.only = TRUE)

options(shiny.maxRequestSize = 3000*1024^2)

server <- function(input, output, session) {
  ##### RAW DATA #####
  # Parse the desc file path and read the table
  desc <- reactive({
    req(input$DESCFILE)
    desc <- read.table(input$DESCFILE$datapath, quote = input$QUOTE, sep = input$SEP, header = input$HEADER, fill = FALSE, as.is=TRUE)
  })
  
  # Unzip and read the data file
  readData <- reactive({
    req(input$DATAFILE)
    unzz <- unzip(input$DATAFILE$datapath)
    print(unzz)
    Data <- ReadAffy(filenames = unzz)
    return(Data)
  }) 
  
  # Table of the description file
  output$desc <- renderTable({
    desc()
  })
  
  # Sample of the raw data
  output$rawdata <- renderTable({
    req(readData())
    head(exprs(readData()),10)
  })
  
  # Function to add standard CDF
  addStandardCDFenv <- function(Data, overwrite=FALSE) {
    #if overwrite is FALSE a cdf environment will be kept if already loaded
    #if overwrite is TRUE it will always be overwritten (unless none is found)
    #the first option would be the regular one, used to add cdf environments where
    #   automatic loading failed
    #the second could be used to set back updated cdf files to the standard ones
    
    #start with r0 just in case this could exist
    rev <- 0
    
    #initial value
    CDFenv <- 0
    
    # recall which cdfName was added, in case no other one is found (set back even
    #   if it does not exist)
    presetCDF <- Data@cdfName
    
    #check whether environment is already correct
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    
    #try the annotation plus cdf as a cdf name
    if ((class(CDFenv)!="environment") | overwrite) {
      CDFenv <- 0 #needed for cases where overwrite is TRUE, but CDFenv already
      # was an environment
      Data@cdfName<-paste(Data@annotation,".cdf",sep="")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
      #if no succes also try without a dot
      if (class(CDFenv)!="environment") {
        Data@cdfName<-paste(Data@annotation,"cdf",sep="")
        suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
      }
    }
    
    # don't run the loop if CDF is already known up front, or correct one has been
    # found
    while((class(CDFenv)!="environment") & (rev < 10)) {
      Data@cdfName<-paste(Data@annotation,".r",rev,"cdf",sep="")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
      rev <- rev + 1
    }
    
    if ((class(CDFenv)!="environment")) {
      Data@cdfName <- presetCDF
      warning("could not automatically retrieve CDF environment for this chip type - object kept as is")
    }
    
    cat("current cdf environment loaded:",Data@cdfName,"\n")
    return(Data)
  }
  
  # Match the descriptions to the rawdata sample names and add standard env
  rawData <- reactive({
    rawData <- readData()
    description <- desc()
    
    file_order <- match(description[,1],sampleNames(rawData))
    if(sum(is.na(file_order)) > 0) stop("file names in data directory and file names in description file do not match")
    if(length(unique(file_order)) < length(file_order)) stop("file names in description file are not unique")
    rawData <- rawData[,file_order]
    
    sampleNames(rawData)<- as.character(description[,2])
    rawData <- addStandardCDFenv(rawData)
    return(rawData)
  })
  
  ##### PREREQUISITES FOR QC PLOTS #####
  
  # Get the array type
  aType <- reactive({
    Data <- rawData()
    aType <- "PMMM"  
    
    # Test whether the dataset is of aType "PM-only"
    mismatches <- mm(Data[,1])  
    
    if(is.null(mismatches)) {
      # mm does not exist
      aType <- "PMonly"
    } else {
      if(sum(is.na(mismatches))>=length(mismatches)/2){ 
        # mm is always NA or there are more NA values in the mm probes than 
        # defined values (assuming these would just be controls)
        aType <- "PMonly"
      } else {
        matches <- pm(Data[,1])
        notNA <- !is.na(mismatches) & !is.na(matches)
        if(sum(mismatches[notNA]!=matches[notNA])==0){
          # MM contains a copy of PM, which indicates a PMonly array
          aType <- "PMonly"
        }
      }
    }
    return(aType)
  })
  
  # Get experiment factors
  experimentFactor <- reactive({
    description <- desc()
    experimentFactor <- factor(description[,3])
  })
  
  # Create colors for plots and legend
  colList <- reactive({
    experimentFactor <- experimentFactor()
    
    #check whether a factor has been provided
    if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")
    
    if(length(levels(experimentFactor))==1) {
      #if there is only one group (or no groups are provided) take equally spread colors over the rainbow palette
      plotColors <- rainbow(length(experimentFactor),s=.8,v=.7)
      #set group legend color to white, as there is not a specific group color
      legendColors <- "white"
    } else {
      #compute the number of colors needed for each class
      tab.tmp <- table(experimentFactor)
      
      #set the two extreme colors for each class
      colors.light <- rainbow(length(levels(experimentFactor)),s=1-sapply(tab.tmp,min,5)*.1)
      colors.dark <- rainbow(length(levels(experimentFactor)),v=1-sapply(tab.tmp,min,5)*.14)
      
      #create the colors to plot, and colors for the legend (average one per experimental group)
      plotColors <- NULL
      legendColors <- NULL
      for(l in 1:length(levels(experimentFactor))) {
        colorFun <- colorRampPalette(c(colors.light[l],colors.dark[l]))
        tmpColors <- colorFun(tab.tmp[l])
        plotColors[experimentFactor==levels(experimentFactor)[l]] <- tmpColors
        legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
      }
    }
    return(list(plotColors=plotColors,legendColors=legendColors))
  })
  
  # Plot colors
  plotColors <- reactive({colList()$plotColors})
  
  # Legend colros
  legendColors <- reactive({colList()$legendColors})
  
  # Plot symbols
  plotSymbols <- reactive({18-as.numeric(experimentFactor())})
  
  # Legend symbols
  legendSymbols <- reactive({sort(plotSymbols(), decreasing=TRUE)})
  
  # Width of all the plots
  WIDTH <- 1000
  
  # Height of all the plots
  HEIGHT <- 1414
  
  # Point sizes for the plots
  POINTSIZE <- 24
  POINTSIZE1 <- 12
  
  # Maximum number of arrays that can be computed
  MAXARRAY <- 41 
  
  # Fit a probe level model if a plot requires it
  rawData.pset <- reactive({
    rawData.pset <- NULL
    req(input$Nuse || input$RLE)
    rawData.pset <- fitPLM(rawData())                     
  })

  # Check for lys
  lys <- reactive({
    lys<-NULL
    if(is.null(lys)) {
      try({calls<-detection.p.val(rawData())$call
      lys<-calls[rownames(calls)[grep("lys.*3",rownames(calls), ignore.case = TRUE)],]
      rm(calls)},TRUE)
      if(!is.null(lys)) {
        if(length(lys) > length(sampleNames(rawData()))) 
        { lys<-lys[1,] }
      }
    }
    return(lys)
  })
  
  # Check for sprep
  sprep <- reactive({
    sprep<-NULL
    if(!is.null(lys) && is.null(sprep)) {
      try(yack <- yaqc(rawData()),TRUE) 
      if(exists("yack")) {
        spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' !
                                               rownames(yack@morespikes), ignore.case = TRUE),])
        sprep<-t(yack@morespikes[spnames,])
      }
    }
    return(sprep)
  })
  
  # Check for quality
  quality<- reactive({
    quality<-NULL
    if(is.null(quality)) try(quality <- qc(rawData()),TRUE)
    return(quality)
  })
  
  # Notify user that prerequisites have been loaded
  output$CDFnotify <- renderPrint({
    req(rawData())
    
    print("Standard CDF is now loaded")
    
    if(is.null(lys())){
      print("Values based on the detection.p.val function cannot be computed for this chip type")
    }
    if(is.null(sprep())){
      print("Values based on the yaqc function cannot be computed for this chip type")
    }
    if(is.null(quality())) {
      print("Plots based on the simpleaffy qc function cannot be created for this chip type")
    }
  })
  
  ##### Sample quality #####
  # Inputs for sample quality plots
  output$sqInputs <- renderUI({
    tagList(
      if(!is.null(sprep()) && !is.null(lys())) {
        checkboxInput(inputId = "sampleprep",label =  "Sample prep controls",value = T)
      },
      if(!is.null(quality())) {
        checkboxInput(inputId = "ratioplot",label = "3'/5' ratio",value = T)
      },
      checkboxInput(inputId = "rnadegplot",label = "RNA degradation",value = T)
    )
  })

  # Sample prep plot
  output$sampleprep<-renderImage({
    req(input$sampleprep)
    lys <- lys()
    sprep <- sprep()
    
    #par(parStart)
    lmin<-min(sprep)
    lmax<-max(sprep)+50
    
    # create the plot           
    RawDataSamplePrepControl<- tempfile(fileext='.png')
    
    png(file = RawDataSamplePrepControl, width=WIDTH, height=HEIGHT, pointsize=POINTSIZE)
    if(length(sampleNames(rawData())) < MAXARRAY){    
      par(mfrow=c(1,1),oma=c(13,0,0,0))
    }else{
      par(mfrow=c(1,1),oma=c(0.1,0,0,0))
    }    
    plot(c(lmin,lmax), type = 'n', ann = FALSE, axes = FALSE, 
         frame.plot = TRUE, ylim = c(lmin, lmax), xlim = c(1,5))
    
    dummy<-apply(cbind(sprep,plotColors(),1:length(sampleNames(rawData()))),1,
                 function(x) {
                   par(new=TRUE);
                   plot(as.numeric(x[1:(length(x)-2)]),
                        type = 'b', ann=FALSE, axes=FALSE, pch='.', cex=4, lwd=2, 
                        col=x[length(x)-1], lty = as.numeric(x[length(x)]),
                        ylim=c(lmin,lmax), xlim=c(1,5))})
    
    title(main="Spike-in Sample Prep controls intensities and calls",
          ylab="Intensity",cex.lab=0.8)
    
    label3<-""
    title1<-"Intensity: OK (Lys < Phe < Thr < Dap for all arrays)"
    title2<-"\nLys Present calls: OK (all Lys are called present)"
    
    # test for Lys < Phe < Thr < Dap
    bad1 <-  colnames(t(sprep))[!(sprep[,1]<sprep[,2] 
                                  & sprep[,2]<sprep[,3] 
                                  & sprep[,3]<sprep[,4])]
    if(length(bad1)>0){    
      title1<-paste("Intensity: not OK (some array",ifelse(length(bad1)>1,"s",""),
                    " do",ifelse(length(bad1)>1,"","es")," not follow the rule Lys < Phe < Thr < Dap;",
                    "\nmaybe no Sample Prep controls were spiked on these arrays.)",sep="")
    }
    
    # test for Lys presence
    bad2 <- names(lys[lys != 'P'])
    bad2 <- gsub(".present","",bad2)
    if(length(bad2)>0){
      title2<-paste("Lys Present calls: ",length(bad2),
                    " Lys not called present",sep="")
    }
    if(length(lys[lys=='P'])==0){
      label3 <- 
        "\nApparently no Sample Prep controls were spiked on these arrays."
    } 
    
    title(xlab = c(title1,title2), cex.lab = 0.8)
    legend("topleft", paste("Lys present calls = ", 
                            length(lys[lys == "P"]),"/", length(lys),label3), bty = "n", cex=0.7)
    if(length(sampleNames(rawData())) < (MAXARRAY+20)){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }         
    legend("topright", substr(sampleNames(rawData()),1,20), lwd=2,    
           col = plotColors(), cex =cexval, bty = "n", lt = 1:length(sampleNames(rawData())))
    par(cex.axis = 0.8)
    axis(2)
    axis(1, at=1:4, labels = c("Lys", "Phe", "Thr", "Dap"))
    
    if(length(bad2)>0){
      text(c(rep(0.87,length(lys[lys!='P']))),sprep[lys!='P',1],
           lys[lys!='P'], pos=4,offset=0.2, cex=0.7, col="red")
    }
    dev.off()
    
    list(src = RawDataSamplePrepControl,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
    
  # 3'5' Ratio plot beta actin
  output$ratioplot<-renderImage({
    req(input$ratioplot)
    quality<- quality()
    #par(parStart)
    
    plotFun <- function(i,j,cutoff,Cname){
      ratio35 <- quality@qc.probes[,i]/quality@qc.probes[,j]
      ratioM <- quality@qc.probes[,i]/quality@qc.probes[,i+1]
      
      cMin <- 0
      cMax <- max(max(cbind(ratio35, ratioM)) + 2, 5)
      if(Cname=="beta-actin"){
        symbol=c(17,2) # triangles for beta-actin
      }else{
        symbol=c(19,1) # circles for GAPDH (as in the simpleaffy QC report)
      }
      
      par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
      plot(ratio35, type='n', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, pch='.', cex=10,  ylim=c(cMin,cMax))
      rect(0, 0, length(ratio35)+1, cutoff, col=gray(0.9), border=FALSE)
      
      for (k in 0:cMax){
        abline(h=k,lty=2,col=gray(0.8))
      }
      
      par(new=TRUE)
      plot(ratio35, type='h', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, lwd=3, col=plotColors(), ylim=c(cMin,cMax))
      par(new=TRUE)
      plot(ratio35,type='p', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, pch=symbol[1], col=plotColors(), ylim=c(cMin,cMax))
      par(new=TRUE)
      plot(ratioM,type='p', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, pch=symbol[2], col="black", ylim=c(cMin,cMax))           
      
      title(main= paste("RNA degradation of", Cname),ylab="3'/5' and 3'/M ratios")    
      axis(2)
      par(cex.axis=0.65)
      if(length(sampleNames(rawData()))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
        axis(1,at=1:length(ratio35),las=2,labels=sampleNames(rawData()))
      }
      if(length(levels(experimentFactor()))>1){
        legend("topright", levels(experimentFactor()),
               col = legendColors(), fill = legendColors(), bty = "n",cex = 0.55)
      }
      legend("topleft", c(paste("3\'/5\' ratio (max=", round(max(ratio35),2),")",
                                sep=""), paste("3'/M ratio (max=", round(max(ratioM),2),")",sep="")),
             col = c(gray(0.3),"black"),pch=symbol, cex=0.55, bty = "o")
      t1 <- ratio35 <= cutoff	
      if(length(t1[t1==FALSE])>0 && length(sampleNames(rawData()))>=(MAXARRAY/2)){
        if(length(t1[t1==FALSE])<25){
          textlab<-sampleNames(rawData())[t1==FALSE]
          legend("topleft", c(rep("",4),"Outliers, from left to right",textlab), 
                 text.col = c(rep(1,5), rep("red",length(textlab))), bty = "n",cex=0.38)  
        }else{
          mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
        }
      }
      mtext(paste("Ratios should stand within the grey rectangle [0,",cutoff,"]"), 
            side=4, font=1, cex=0.7)
      
      par(cex.axis=0.8)
      boxplot(cbind(ratio35, ratioM), axes=FALSE, frame.plot=TRUE)
      title(main=paste("Boxplot of", Cname,"ratios"))    
      par(cex.lab=0.8)
      if(max(ratio35 ) < cutoff){
        title(xlab=paste(Cname," QC: OK (all 3'/5' ratios < ",cutoff,")",sep=""))
      }  else{
        title(xlab=paste(Cname," QC: not OK (some 3'/5' ratios >",cutoff,
                         ")\nNote that the threshold of ",cutoff,
                         " \nwas determined for Homo Sapiens.",sep=""), line=4)
      }       
      
      axis(1,at=1:2,labels=c("3'/5' ratio","3'/M ratio"))
      axis(2)
    }
    
    # Graph creation: beta-actin
    RawData53ratioPlot_betaActin<-tempfile(fileext = '.png')
    png(file = RawData53ratioPlot_betaActin,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
    plotFun(1,3,3,"beta-actin")
    dev.off()
    
    list(src = RawData53ratioPlot_betaActin,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # 3'5' Ratio plot GADPH
  output$ratioplot2<-renderImage({
    req(input$ratioplot)
    quality<-quality()
    #par(parStart)
    
    plotFun <- function(i,j,cutoff,Cname){
      ratio35 <- quality@qc.probes[,i]/quality@qc.probes[,j]
      ratioM <- quality@qc.probes[,i]/quality@qc.probes[,i+1]
      
      cMin <- 0
      cMax <- max(max(cbind(ratio35, ratioM)) + 2, 5)
      if(Cname=="beta-actin"){
        symbol=c(17,2) # triangles for beta-actin
      }else{
        symbol=c(19,1) # circles for GAPDH (as in the simpleaffy QC report)
      }
      
      par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
      plot(ratio35, type='n', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, pch='.', cex=10,  ylim=c(cMin,cMax))
      rect(0, 0, length(ratio35)+1, cutoff, col=gray(0.9), border=FALSE)
      
      for (k in 0:cMax){
        abline(h=k,lty=2,col=gray(0.8))
      }
      
      par(new=TRUE)
      plot(ratio35, type='h', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, lwd=3, col=plotColors(), ylim=c(cMin,cMax))
      par(new=TRUE)
      plot(ratio35,type='p', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, pch=symbol[1], col=plotColors(), ylim=c(cMin,cMax))
      par(new=TRUE)
      plot(ratioM,type='p', ann=FALSE, axes=FALSE, 
           frame.plot=TRUE, pch=symbol[2], col="black", ylim=c(cMin,cMax))           
      
      title(main= paste("RNA degradation of", Cname),ylab="3'/5' and 3'/M ratios")    
      axis(2)
      par(cex.axis=0.65)
      if(length(sampleNames(rawData()))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
        axis(1,at=1:length(ratio35),las=2,labels=sampleNames(rawData()))
      }
      if(length(levels(experimentFactor()))>1){
        legend("topright", levels(experimentFactor()),
               col = legendColors(), fill = legendColors(), bty = "n",cex = 0.55)
      }
      legend("topleft", c(paste("3\'/5\' ratio (max=", round(max(ratio35),2),")",
                                sep=""), paste("3'/M ratio (max=", round(max(ratioM),2),")",sep="")),
             col = c(gray(0.3),"black"),pch=symbol, cex=0.55, bty = "o")
      t1 <- ratio35 <= cutoff	
      if(length(t1[t1==FALSE])>0 && length(sampleNames(rawData()))>=(MAXARRAY/2)){
        if(length(t1[t1==FALSE])<25){
          textlab<-sampleNames(rawData())[t1==FALSE]
          legend("topleft", c(rep("",4),"Outliers, from left to right",textlab), 
                 text.col = c(rep(1,5), rep("red",length(textlab))), bty = "n",cex=0.38)  
        }else{
          mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
        }
      }
      mtext(paste("Ratios should stand within the grey rectangle [0,",cutoff,"]"), 
            side=4, font=1, cex=0.7)
      
      par(cex.axis=0.8)
      boxplot(cbind(ratio35, ratioM), axes=FALSE, frame.plot=TRUE)
      title(main=paste("Boxplot of", Cname,"ratios"))    
      par(cex.lab=0.8)
      if(max(ratio35 ) < cutoff){
        title(xlab=paste(Cname," QC: OK (all 3'/5' ratios < ",cutoff,")",sep=""))
      }  else{
        title(xlab=paste(Cname," QC: not OK (some 3'/5' ratios >",cutoff,
                         ")\nNote that the threshold of ",cutoff,
                         " \nwas determined for Homo Sapiens.",sep=""), line=4)
      }       
      
      axis(1,at=1:2,labels=c("3'/5' ratio","3'/M ratio"))
      axis(2)
    }
    
    # Graph creation: two separated graphs
    RawData53ratioPlot_GADPH<-tempfile(fileext = '.png')
    png(file = RawData53ratioPlot_GADPH,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)  
    plotFun(4,6,1.25,"GAPDH")
    dev.off()  
    list(src = RawData53ratioPlot_GADPH,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # RNA degradation plot
  output$rnadegplot<- renderImage({
    req(input$rnadegplot)
    Data.rnadeg=NULL
    
    RawDataRNAdegradation <- tempfile(fileext='.png')
    
    if(is.null(plotColors())) stop("the 'plotColors' parameter is required") 
    if(is.null(Data.rnadeg)) Data.rnadeg <- AffyRNAdeg(rawData())
    
    png(file = RawDataRNAdegradation, width=WIDTH, height=HEIGHT, pointsize=POINTSIZE)
    
    if(length(sampleNames(rawData()))<MAXARRAY){
      cexval <- 0.65
      par(lwd = 2, oma=c(13,0,0,0))	
    } else{
      cexval <- 0.45 
      par(lwd = 2, oma=c(0.1,0,0,0))	
    }    
    par(mar=c(4,4,4,0), cex.axis=0.6, cex.lab=0.75)
    layout(matrix(c(1,2),1,2,byrow=TRUE), c(2,1), 1, FALSE)
    plotAffyRNAdeg(Data.rnadeg, col = plotColors(), lty = 1:length(sampleNames(rawData()))) 
    par(mar=c(4,0,4,0), cex.axis=0.01, cex.lab=0.01)
    plot(1,type = 'n', ann=FALSE, axes=FALSE, frame.plot=FALSE)
    legend("topleft",sampleNames(rawData()), lwd=2, col = plotColors(), cex = cexval, bty = "n", lty = 1)
    dev.off() 
    
    # Return a list containing the file attributes
    list(src = RawDataRNAdegradation,
         #  width = WIDTH,
         #   height = HEIGHT,
         alt = "This is alternate text")
  })
 
   
  ##### Hybridization and signal quality #####
  # Inputs for hybridization and signal quality
  output$hybridInput <- renderUI({
    tagList(
      if(!is.null(quality())) {
        checkboxInput(inputId = "hybrid",label = "Spike-In controls",value = T)
        checkboxInput(inputId = "percentpresent",label = "Percent present",value = T)
        checkboxInput(inputId = "backgroundintensity",label = "Background intensity",value = T)
      },
      checkboxInput(inputId = "posnegdistro",label = "Positive/Negative controls",value = T)
    )
  })
  
  # Spike in hybridization
  output$hybrid<-renderImage({
    req(input$hybrid)
    RawDataSpikeinHybridControl<- tempfile(fileext='.png')
    
    png(file = RawDataSpikeinHybridControl,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
    if(length(sampleNames(rawData())) < MAXARRAY){    
      par(mfrow=c(1,1),oma=c(13,0,0,0))
    }else{
      par(mfrow=c(1,1),oma=c(0.1,0,0,0))
    }
    plot(c(min(quality()@spikes),max(quality()@spikes)),
         type = 'n', ann = FALSE, axes = FALSE, frame.plot = TRUE,
         ylim = c(min(quality()@spikes), max(quality()@spikes)),
         xlim = c(1,5))
    
    dummy<-apply(cbind(quality()@spikes,plotColors(),1:length(sampleNames(rawData()))),1,
                 function(x) {
                   par(new=TRUE);
                   plot(as.numeric(x[1:(length(x)-2)]),
                        type = 'b', ann=FALSE, axes=FALSE, pch='.', cex=4, lwd=2, 
                        col=x[length(x)-1], lty = as.numeric(x[length(x)]),
                        ylim=c(min(quality()@spikes),max(quality()@spikes)), xlim = c(1,5))})
    
    title(main = "Spike-in Hybridization controls intensities and calls",
          ylab="Intensity",cex.lab=0.8) 
    label3<-""
    title1<-"Intensities: OK (bioB < bioC < bioD < creX for all arrays)"
    title2<-"BioB Present calls: OK (indeed all bioB are called present)"
    
    # test for bioB<bioC<bioD<creX
    bad1 <- colnames(t(quality()@spikes))[!(quality()@spikes[,1]<quality()@spikes[,2] 
                                            & quality()@spikes[,2]<quality()@spikes[,3] 
                                            & quality()@spikes[,3]<quality()@spikes[,4])]
    if(length(bad1)>0){    
      title1<-paste("Intensity: not OK (some array",ifelse(length(bad1)>1,"s",""),
                    " do",ifelse(length(bad1)>1,"","es")," not follow the rule bioB < bioC < bioD < creX)",sep="")
    }
    
    # test for bioB presence
    bad2<-gsub(".present","",names(quality()@bioBCalls[quality()@bioBCalls != 'P']))
    if(length(bad2)>0){
      title2<-paste("BioB present calls: not OK (",length(bad2),"/",
                    length(quality()@bioBCalls)," bioB not called present)")
    }
    if(length(quality()@bioBCalls[quality()@bioBCalls =='P'])==0){
      label3 <- 
        "\nApparently no Hybridization controls were spiked on these arrays."
    }
    
    title(xlab = c(title1,title2), cex.lab = 0.8)
    legend("topleft", paste("bioB present calls = ", 
                            length(quality()@bioBCalls[quality()@bioBCalls == "P"]),"/",
                            length(quality()@bioBCalls),label3), bty = "n", cex=0.7)
    
    if(length(sampleNames(rawData())) < (MAXARRAY+20)){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }              
    
    legend("bottomright", substr(sampleNames(rawData()),1,20), lwd=2,
           col = plotColors(), cex = cexval, bty = "n", lt = 1:length(sampleNames(rawData())))
    par(cex.axis = 0.8) 
    axis(2)
    axis(1, at=1:4, labels = names(colnames(quality()@spikes)))
    
    if(length(bad2)>0){
      text(c(rep(0.89,length(quality()@bioBCalls[quality()@bioBCalls!='P']))),
           quality()@spikes[quality()@bioBCalls!='P',1], 
           quality()@bioBCalls[quality()@bioBCalls!='P'], pos=4,offset=0.2, 
           cex=0.8, col="red")
    }  
    dev.off()
    
    list(src = RawDataSpikeinHybridControl,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Percent Present
  output$percentpresent<-renderImage({
    req(input$percentpresent)
    # Define the reference for the grey rectangle (it should minimize the number of outliers)
    ppMin <- min(quality@percent.present)
    ppMax <- max(quality@percent.present)
    ppMean<- mean(quality@percent.present)
    
    ref<-cbind(c(ppMean-5,ppMin,ppMax-10),c(ppMean+5,ppMin+10,ppMax))
    reftext<-c("[mean-5% ; mean+5%]", "[min ; min+10%]", "[max-10% ; max]")
    t1 <- (quality@percent.present >= (ppMean-5) & quality@percent.present <= (ppMean+5))
    t2 <- (quality@percent.present >= ppMin & quality@percent.present <= (ppMin+10))	
    t3 <- (quality@percent.present >= (ppMax-10) & quality@percent.present <= ppMax)	
    outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
    if(outliers[1] == 0){
      ref<- ref[1,]
      reftext<- reftext[1]
    }else{
      ref<- ref[outliers==min(outliers),] # less outliers
      if(length(ref) > 3) { ref<-ref[1,] }
      reftext<- reftext[outliers==min(outliers)] # less outliers
      if(length(reftext) > 1) { reftext<-reftext[1] }		
    }
    testMin <- ref[1]
    testMax <- ref[2]
    
    RawDataPercentPresent<- tempfile(fileext='.png')
    png(file = RawDataPercentPresent,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
    plot(quality@percent.present, type='n', ann=FALSE, axes=FALSE, 
         frame.plot=TRUE, pch='.', cex=10,  ylim=c(0,100))
    rect(0, testMin, length(quality@scale.factors)+1,
         testMax, col=gray(0.9), border=FALSE)
    
    for(i in seq(from=0, to=100, by=20)){
      abline(h=i,lty=2,col=gray(0.8))
    }
    abline(h=50,lty=1,col=gray(0.5))
    
    par(new=TRUE)
    plot(quality@percent.present, type='h', ann=FALSE, axes=FALSE, 
         frame.plot=TRUE, pch='.',lwd=3, col=plotColors(), ylim=c(0,100))
    par(new=TRUE)
    plot(quality@percent.present,type='p', ann=FALSE, axes=FALSE, 
         frame.plot=TRUE, pch='.',cex=10, col=plotColors(), ylim=c(0,100))
    
    t1 <- (quality@percent.present >= testMin & quality@percent.present <= testMax)
    if(length(t1[t1==FALSE])>0 && length(sampleNames(rawData()))>=(MAXARRAY/2)){
      if(length(t1[t1==FALSE])<25){
        textlab<-sampleNames(rawData())[t1==FALSE]
        legend("topleft", c("Outliers, from left to right",textlab), 
               text.col = c(1, rep("red",length(textlab))), bty = "n",cex=0.38)  	
      }else{
        mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
      }			  
    }	
    
    title(main="Plot of percent present",xlab="",ylab="percentage")
    axis(2)	
    par(cex.axis=0.65)
    if(length(sampleNames(rawData()))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
      axis(1,at=1:length(quality@percent.present),las=2,labels=sampleNames(rawData()))
    }	
    
    if(length(levels(experimentFactor()))>1){
      legend("topright", levels(experimentFactor()),
             col = legendColors(), fill = legendColors(), cex = 0.55, bty = "n")
    }
    
    legend("bottomright", c(paste("min = ", round(ppMin,2), sep="" ), 
                            paste("max = ", round(ppMax,2), sep="" ), 
                            paste("max-min = ", round(ppMax-ppMin,2), sep="")), bty = "o",cex = 0.55)
    
    mtext(paste(
      "Data should be in the grey rectangle representing a spread of 10%"
      #"Data should stand within the grey rectangle",reftext
    ), 
    side=4, font=1, cex=0.7)
    par(cex.axis=0.8, cex.lab=0.8)       
    boxplot(quality@percent.present)
    title(main="Boxplot of percent present")
    if(ppMax-ppMin <= 10){
      title(xlab="Percent present QC: OK (spread <= 10%)")
    } else{
      title(xlab="Percent present QC: not OK (spread > 10%)")
    } 
    
    dev.off()
    
    list(src = RawDataPercentPresent,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Positive/Negative Distribution
  output$posnegdistro<-renderImage({
    req(input$posnegdistro)
    RawDataPosNegDistribution<- tempfile(fileext='.png')
    png(file = RawDataPosNegDistribution,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(oma=c(17,0,0,0),srt=90)  
    borderQC1(rawData())
    mtext(paste("All distributions should be similar and","extreme values should not be reached\n"), side=3, cex=0.7)
    dev.off()
    list(src = RawDataPosNegDistribution,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Background intensities
  output$backgroundintensity<-renderImage({
    req(input$backgroundintensity)
    
    # Define the reference for the grey rectangle (it should minimize the number of outliers)
    bgMean<-mean(quality()@average.background)
    bgMin <- min(quality()@average.background)
    bgMax <- max(quality()@average.background)
    ref<-cbind(c(bgMean-10,bgMin,bgMax-20),c(bgMean+10,bgMin+20,bgMax))
    reftext<-c("[mean-10 ; mean+10]", "[min ; min+20]", "[max-20 ; max]")
    t1 <- (quality()@average.background >= (bgMean-10) & quality()@average.background <= (bgMean+10))
    t2 <- (quality()@average.background >= bgMin & quality()@average.background <= (bgMin+20))	
    t3 <- (quality()@average.background >= (bgMax-20) & quality()@average.background <= bgMax)	
    outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
    if(outliers[1] == 0){
      ref<- ref[1,]
      reftext<-reftext[1]
    }else{
      ref<- ref[outliers==min(outliers),] # less outliers
      if(length(ref) > 3) { ref<-ref[1,] }
      reftext<- reftext[outliers==min(outliers)] # less outliers
      if(length(reftext) > 1) { reftext<-reftext[1] }				
    }
    testMin <- ref[1]
    testMax <- ref[2]
    
    ylimit = c(0, max(max(quality()@maximum.background), 
                      testMax)+5)
    
    RawDataBackground<- tempfile(fileext='.png')
    png(file = RawDataBackground,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
    plot(c(quality()@minimum.background, quality()@maximum.background),
         type = 'n', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', cex=10, 
         xlim = c(1,length(quality()@maximum.background+1)),
         ylim = ylimit)
    rect(0,testMin, 
         length(quality()@average.background)+1, 
         testMax, 
         col=gray(0.8), border=TRUE, density=10, lty=6)           
    par(new=TRUE)
    plot(quality()@maximum.background,
         type = 'h', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', lwd=3, 
         col=plotColors(), ylim = ylimit)
    par(new=TRUE)
    plot(quality()@minimum.background,
         type = 'h', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', lwd=4, 
         col='white', ylim = ylimit)        
    par(new=TRUE)
    plot(quality()@maximum.background,
         type = 'p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=6, cex=0.5,
         ylim = ylimit)
    par(new=TRUE)
    plot(quality()@minimum.background,
         type = 'p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=2, cex=0.5,
         ylim=ylimit)
    par(new=TRUE)
    plot(quality()@average.background, 
         type='p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=3, cex=1,
         ylim=ylimit)
    abline(h=testMin, col=gray(0.8), ylim=ylimit)
    abline(h=testMax, col=gray(0.8), ylim=ylimit)	
    title(main="Plot of background intensity",xlab="",ylab="background intensity")
    axis(2)
    par(cex.axis=0.65)
    if(length(sampleNames(rawData()))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
      axis(1,at=1:length(quality()@average.background),las=2,labels=sampleNames(rawData()))
    }		
    if(length(levels(experimentFactor()))>1){    
      legend("bottomright", c(levels(experimentFactor()),"max bg","average bg", 
                              "min bg"), col = c(legendColors(), 1, 1, 1), 
             pch = c(rep(15, length(legendColors())),6,3,2), bty = "n", cex=0.55)
    }else{
      legend("bottomright", c("max bg","average bg", 
                              "min bg"), col = c(1, 1, 1), pch = c(6,3,2), bty = "n", cex=0.55)    
    }
    t1 <- (quality()@average.background >= testMin & quality()@average.background <= testMax)
    if(length(t1[t1==FALSE])>0 && length(sampleNames(rawData()))>=(MAXARRAY/2)){
      if(length(t1[t1==FALSE])<25){
        textlab<-sampleNames(rawData())[t1==FALSE]
        legend("bottomleft", c("Outliers, from left to right",textlab), 
               text.col = c(1, rep("red",length(textlab))), bty = "n",cex=0.38)   
      }else{
        mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
      }	 		  
    }	
    
    mtext(paste(
      "Data should be in the grey rectangle representing a spread of 20"
      #"Data should stand within the dashed grey rectangle",reftext
    ), side=4, font=1, cex=0.7)
    par(cex.axis=0.8, cex.lab=0.8)            
    boxplot(quality()@average.background)
    title(main="Average background intensity")
    
    if(bgMax - bgMin <= 20){
      title(xlab="Background QC: OK (spread <= 20)")
    } else{
      title(xlab="Background QC: not OK (spread > 20)")
    }     
    
    legend("bottomleft", c(paste("min = ", round(bgMin,2), sep="" ), 
                           paste("max = ", round(bgMax,2), sep="" ), 
                           paste("max-min = ", round(bgMax-bgMin,2), sep="")), cex=0.7)        
    
    dev.off()
    
    list(src = RawDataBackground,width = WIDTH,height = HEIGHT,alt = "This is alternate text") 
  })
  
  ##### Signal comparability and bias diagnostic #####
  # Inputs
  output$signalInput <- renderUI({
    tagList(
      if(!is.null(quality())) {
        checkboxInput(inputId = "scalefactors",label = "Scale Factors",value =  T)
      },
      checkboxInput(inputId = "boxplotRaw",label = "Boxplot",value = T),
      checkboxInput(inputId = "densityRaw",label = "Density histogram",value =  T),
      checkboxInput(inputId = "controlplot",label = "Control profiles and affx boxplots",value = T)
    )
  })
  
  # Scale factors
  output$scalefactors<-renderImage({
    req(input$scalefactors)
    
    # Define the reference for the grey rectangle (it should minimize the number of outliers)
    sfMin <- min(log2(quality()@scale.factors))
    sfMax <- max(log2(quality()@scale.factors))
    sfMean<- mean(log2(quality()@scale.factors))
    ref<-cbind(c(sfMean-1.5,sfMin,sfMax-3),c(sfMean+1.5,sfMin+3,sfMax))
    t1 <- (log2(quality()@scale.factors) >= (sfMean-1.5) & log2(quality()@scale.factors) <= (sfMean+1.5))
    t2 <- (log2(quality()@scale.factors) >= sfMin & log2(quality()@scale.factors) <= (sfMin+3))	
    t3 <- (log2(quality()@scale.factors) >= (sfMax-3) & log2(quality()@scale.factors) <= sfMax)	
    outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
    if(outliers[1] == 0){
      ref<- ref[1,]
    }else{
      ref<- ref[outliers==min(outliers),] # less outliers
      if(length(ref) > 3) { ref<-ref[1,] }		
    }
    testMin <- ref[1]
    testMax <- ref[2]
    
    ymin <- min(-3.5, sfMin-1.5, testMin)
    ymax <- max(+3.5, sfMax+1.5, testMax) 
    ylimit <-c(ymin, ymax) 
    
    RawDataScaleFactors<- tempfile(fileext='.png')
    png(file = RawDataScaleFactors,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    
    par(mfrow=c(1,2), oma=c(17,0.5,0,0.5), cex.axis=0.8)
    plot(log2(quality()@scale.factors), type='n', axes=FALSE, frame.plot=TRUE,
         ann=FALSE, pch='.', cex=10, ylim = ylimit)
    rect(0, testMin,
         length(quality()@scale.factors)+1, 
         testMax, col=gray(0.9), border=FALSE)
    
    for(i in seq(from=ceiling(ymin), to=floor(ymax), by=1)){
      abline(h=i,lty=2,col=gray(0.8))
    }
    abline(h=0,lty=1,col=gray(0.5))
    
    par(new=TRUE)
    plot(log2(quality()@scale.factors), 
         type = 'h', axes=FALSE, frame.plot=TRUE, ann=FALSE,
         pch='.', lwd=3, col=plotColors(), ylim = ylimit)
    par(new=TRUE)
    plot(log2(quality()@scale.factors), 
         type = 'p', axes=FALSE, frame.plot=TRUE, ann=FALSE,pch='.',
         cex=10, col=plotColors(), ylim = ylimit)
    t1 <- (log2(quality()@scale.factors) >= testMin & log2(quality()@scale.factors) <= testMax)
    if(length(t1[t1==FALSE])>0 && length(sampleNames(rawData()))>=(MAXARRAY/2)){
      if(length(t1[t1==FALSE])<16){
        textlab<-sampleNames(rawData())[t1==FALSE]
        legend("bottomleft", c("Outliers, from left to right",textlab), 
               text.col = c(1, rep("red",length(textlab))), bty = "n",cex=0.38)  	
      }else{
        mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
      }	 		  
    }
    
    title(main="Plot of Log scale factors", xlab="",ylab="Log2(scale factors)")
    axis(2)
    par(cex.axis=0.65)
    if(length(sampleNames(rawData()))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
      axis(1,at=1:length(quality()@scale.factors),las=2,labels=sampleNames(rawData()))
    }		
    if(length(levels(experimentFactor()))>1){ 
      legend("topright", levels(experimentFactor()),
             col = legendColors(), fill = legendColors(), cex = 0.55, bty = "n")
    } 
    legend("topleft", c(paste("min = ", round(sfMin,2), sep="" ), 
                        paste("max = ", round(sfMax,2), sep="" ), 
                        paste("max - min = ",round((sfMax-sfMin),2),sep="")), bty = "n",
           cex = 0.55)
    
    mtext(paste(
      "Data should be in the grey rectangle representing 3-fold on a log scale"
    ), side=4,  font=1, cex=0.7)
    par(cex.axis=0.8)       
    boxplot(quality()@scale.factors)
    title(main="Boxplot of scale factors")
    mtext("(natural scale)\n", side=3, font=1, cex=0.7)
    par(cex=0.8)
    if((sfMax - sfMin) < 3){
      
      title(xlab="Scale factors QC: OK (spread < 3-fold)")
    } else{
      title(xlab="Scale factors QC: not OK (spread > 3-fold)")
    }         
    dev.off()
    
    list(src = RawDataScaleFactors,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Raw Box plot
  output$boxplotRaw<-renderImage({
    req(input$boxplotRaw)
    
    Type <- "Raw"
    tmain <- "Boxplot of raw intensities"
    tmtext2 <- "Raw log intensity\n\n\n"
    
    DataBoxplot<- tempfile(fileext='.png')
    png(file = DataBoxplot,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)  
    par(oma=c(17,0,0,0), cex.axis=1) 
    suppressWarnings(boxplot(rawData(), col=plotColors() ,main=tmain, axes=FALSE, pch = 20, cex=0.7))
    if(length(levels(experimentFactor()))>1){ 
      legend("topright", levels(experimentFactor()),
             col=legendColors(),fill=legendColors(), cex = 0.7, bg = "white", bty = "o")
    }
    if(length(sampleNames(rawData()))<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }  
    axis(1,at=1:length(sampleNames(rawData())),las=2,labels=sampleNames(rawData()), cex.axis=cexval)        
    axis(2, cex.axis=0.7)  
    mtext(tmtext2, side=2, cex=0.8)  	   
    mtext("Distributions should be comparable between arrays\n", side=3, font=1, 
          cex=0.7)
    dev.off()
    
    list(src = DataBoxplot,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Density Histogram
  output$densityRaw<-renderImage({
    req(input$densityRaw)
     
    Type <- ifelse(class(rawData()) == "AffyBatch","Raw","Norm")
    
    DensityHistogram<- tempfile(fileext='.png')
    png(file = DensityHistogram,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    if(length(sampleNames(rawData()))<MAXARRAY){
      cexval <- 0.65
      par(oma=c(12,0,0,0) )
    }else{
      cexval <- 0.45
      par(oma=c(0.1,0,0,0) )
    }    
    if(Type == "Raw"){
      hist(rawData(), lwd=3, lt = 1:length(sampleNames(rawData())), col = plotColors(), 
           which = "both", main="Density histogram of raw intensities", cex.axis = 0.7, cex.lab=0.8)
    }else{
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      hist(rawData(), lwd=3, lt = 1:length(sampleNames(rawData())), col = plotColors(), 
           main=paste("Density histogram after ", normMeth,"\n", sep=""), cex.axis = 0.7, cex.lab=0.8)
    } 
    
    legend("topright", substr(sampleNames(rawData()),1,20), lwd=3, lt = 1:length(sampleNames(rawData())),
           col = plotColors(), cex = cexval, bty = "n")           
    mtext( "Curves should be comparable between arrays\n", side=3, font=1, 
           cex=0.7)
    dev.off()
    list(src = DensityHistogram,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Control profiles affx
  output$controlplotaffx <- renderImage({
    req(input$controlplot)
    
    affx1 <- unique(probeNames(rawData())[grep("AFFX",probeNames(rawData()),fixed=TRUE)])
    
    dataTable <- paste(substr(rawData()@annotation,1,nchar(rawData()@annotation)-2),"CONTROL",sep="")
    suppressWarnings(eval(parse("",-1,paste("data(",dataTable,")",sep="")))) #ArrayTools
    cntrl <- NULL
    try(cntrl <- get(dataTable),TRUE)
    
    cntrlVars <- NULL
    if(!is.null(cntrl)) {
      for(name in names(table(cntrl[,2])[table(cntrl[,2])>0])) {
        #often names have a format of xxx->yyy, we only want the last part
        name2 <- strsplit(name,"->")[[1]][length(strsplit(name,"->")[[1]])]
        #make sure that there are only letters in name2, to be used as variable name
        name2 <- gsub("[[:punct:]|[:digit:]]","",name2)
        assign(name2,cntrl[grep(name,cntrl[,2]),1])
        assign(name2,get(name2)[get(name2) %in% geneNames(rawData())])
        if(length(get(name2))>0) cntrlVars <- c(cntrlVars,name2)
      }
    }
      
      if(length(affx1)>0 | length(grep("affx",cntrlVars))>0) {
        
        if(length(grep("affx",cntrlVars))>0) {
          #affx2 cannot exists as automatically generated variable names are restricted to letters
          affx2 <- NULL
          for(a in grep("affx",cntrlVars)) {
            affx2 <- c(affx2, get(cntrlVars[a]))
          }
          if(length(affx1)>0) { #both exist
            if(sort(affx2)!=sort(affx1)) {
              #merge both
              affx2 <- unique(c(affx1,affx2))
            }
          }
        } else { #only length(affx1)>0
          affx2 <- affx1
        }		
        #		dimnCoef <- 1 # squared layout
        dimnCoef <- 2
        dimn<-ceiling(sqrt(length(affx2)/dimnCoef))
        number_empty <- dimnCoef*dimn*dimn - length(affx2) 	
        
        RawDataAFFXControlsProfiles <- tempfile(fileext='.png')
        png(file=RawDataAFFXControlsProfiles, width=250*dimn+250, height=dimnCoef*250*dimn, pointsize=POINTSIZE)
        if(number_empty >= dimn){
          lineout<-dimn	
          number_empty <- number_empty - dimn
        }else{
          lineout<-0
        }
        legcol<-dimnCoef*dimn*dimn-lineout+1
        layout(cbind(matrix(1:(dimnCoef*dimn*dimn-lineout),ncol=dimn,byrow=TRUE),legcol))
        par(oma=c(5,3,4,3))
        
        dummy<-sapply(
          affx2,
          function(y) {
            pm<-log(pm(rawData(),as.character(y)),2)
            max_pm<-max(pm)
            min_pm<-min(pm)
            plot(0,type='n',ylim=c(min_pm,max_pm),ylab="expression",xlim=c(1,dim(pm)[1]),
                 main=y,xlab=paste("\naverage: ",round(mean(pm))),cex.lab=1.3,cex.main=1.1)
            assign("col", 0, 1)
            apply(
              pm,
              2,
              function(z) {
                assign("col",get("col",1)+1,1)
                points(z,type='o',col=plotColors()[get("col",1)])
              }
            )
            points(apply(pm,1,mean),type='l',lwd=2,col="black")
          }
        )
        
        
        for(p in 1:number_empty) plot(0,type='n',xaxt='n',yaxt='n',bty='n',xlab="",ylab="",main="")
        par(mar=c(0,0,0,0))
        plot(0,type='n',xaxt='n',yaxt='n',bty='n',xlab="",ylab="",main="")
        legend("topleft",c(sampleNames(rawData()),"","average"),lwd=c(rep(1,length(sampleNames(rawData()))),0,2),
               pch=c(rep(1,length(sampleNames(rawData()))),-1,-1),col=c(plotColors(),"white","black"),cex=1.5,bty="n") 
        
        mtext("affx control profiles", side = 3, outer = TRUE, font = 2, cex = 2)
        
        dev.off()
        list(src = RawDataAFFXControlsProfiles,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
      }
  }, deleteFile = T)
  
  # Control profiles boxplot
  output$controlplotbox <- renderImage({
    req(input$controlplot)
    
    affx1 <- unique(probeNames(rawData())[grep("AFFX",probeNames(rawData()),fixed=TRUE)])
    
    dataTable <- paste(substr(rawData()@annotation,1,nchar(rawData()@annotation)-2),"CONTROL",sep="")
    suppressWarnings(eval(parse("",-1,paste("data(",dataTable,")",sep="")))) #ArrayTools
    cntrl <- NULL
    try(cntrl <- get(dataTable),TRUE)
    
    cntrlVars <- NULL
    if(!is.null(cntrl)) {
      for(name in names(table(cntrl[,2])[table(cntrl[,2])>0])) {
        #often names have a format of xxx->yyy, we only want the last part
        name2 <- strsplit(name,"->")[[1]][length(strsplit(name,"->")[[1]])]
        #make sure that there are only letters in name2, to be used as variable name
        name2 <- gsub("[[:punct:]|[:digit:]]","",name2)
        assign(name2,cntrl[grep(name,cntrl[,2]),1])
        assign(name2,get(name2)[get(name2) %in% geneNames(rawData())])
        if(length(get(name2))>0) cntrlVars <- c(cntrlVars,name2)
      }
    }
    
      if(length(affx1)>0) {
        #affx1 will not already be in controlVars as it contains a number
        cntrlVars <- c("affx1", cntrlVars)
      }
      n <- 0
      for(var in cntrlVars) {
        varName <- paste(var,"Data",sep="")
        assign(varName,NULL)
        try(assign(varName,log(pm(rawData(),as.character(get(var))),2)),TRUE)
        if(!is.null(get(varName))) n<-n+1
      }
      #make boxplots of controls
      if(n > 0) {
        for(i in 1:n) {
          if(!is.null(get(varName))) {
            RawDataAFFXControlsBoxplot <- tempfile(fileext = '.png')
            png(RawDataAFFXControlsBoxplot, width=WIDTH, height=HEIGHT,pointsize=POINTSIZE)
            par(oma=c(17,0,0,0))
            boxplot(get(paste(cntrlVars,"Data",sep="")[i]), main=paste(gsub("1","",
                                                                            cntrlVars[i],fixed=TRUE),"controls"), ylim=c(0,16), col=plotColors(), cex=1, axes=FALSE)
            if(length(sampleNames(rawData()))<MAXARRAY){
              cexval <- 0.65
            }else{
              cexval <- 0.45
            }
            if(length(levels(experimentFactor))>1){ 
              legend("topright", levels(experimentFactor), col = legendColors, 
                     fill = legendColors, cex = 0.7, bg = "white", bty = "o")
            }     					
            axis(1,at=1:length(sampleNames(rawData())),las=2,labels=sampleNames(rawData()), cex.axis=cexval)        
            axis(2, cex.axis=0.7)  
            dev.off()
          }
        }
      }
    list(src = RawDataAFFXControlsBoxplot,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  }, deleteFile = T)
  
  # Layout Plot
  output$layoutplot<-renderImage({
    req(input$layoutplot)
    
    #copy first array as a model
    data.tmp <- rawData()[,1]
    
    #set every probe to intensity 1
    exprs(data.tmp)[] <- 1
    
    #determine which are the perfect match rows in exprs
    rn_pm <- sort(as.numeric(rownames(pm(data.tmp))))
    rn_mm <- NULL
    if(aType()=="PMMM") {
      #determine which are the mismatch rows in exprs if present
      rn_mm <- sort(as.numeric(rownames(mm(data.tmp))))
      #set intensity of all mismatch probes if present to 2
      exprs(data.tmp)[rn_mm,] <- 2
    }
    
    #determine control probes
    affx <- grep("AFFX",probeNames(data.tmp),fixed=TRUE)
    
    library <- paste(substr(rawData()@annotation,1,nchar(rawData()@annotation)-2),"transcriptcluster.db",sep="")
    suppressWarnings(eval(parse("",-1,paste("require(",library,", quietly = TRUE)",sep=""))))
    #get probe annotations
    all<-NULL
    eval(parse("",-1,paste("try(all2<-",substr(library,1,nchar(library)-3),"ACCNUM,TRUE)",sep="")))
    if(exists("all2")) {
      eval(parse("",-1,paste("try(all<-toTable(all2),TRUE)",sep="")))
    }
    
    dataTable <- paste(substr(rawData()@annotation,1,nchar(rawData()@annotation)-2),"CONTROL",sep="")
    suppressWarnings(eval(parse("",-1,paste("data(",dataTable,")",sep="")))) #ArrayTools
    cntrl <- NULL
    try(cntrl <- get(dataTable),TRUE)
    
    #controls<-all[all[,2]=="---",1]
    #sum(controls!=sort(cntrl[,1])) #should be 0, so all IDs are equal
    
    lab_add <- ""
    if(length(affx)>0) {
      #set intensity of control probes to 5 for matches and to 8 for mismatches if present
      pm(data.tmp)[affx,] <- 5
      if(aType()=="PMMM") mm(data.tmp)[affx,] <- 8 
    }
    
    #if no cntrl object given: warn and skip plot
    if(!is.null(cntrl)) {
      #set control probes to intensity 2 (affx, intron, exon)
      #these normally are the affx, intron, and exon probes as ag controls are not present in pm
      pm(data.tmp)[probeNames(data.tmp) %in% cntrl[,1],] <- 5
      if(aType()=="PMMM") mm(data.tmp)[probeNames(data.tmp) %in% cntrl[,1],] <- 8
    }
    
    #if given, set probes missing from list of annotated probes to NA
    if(!is.null(all)) {
      missingAnnotation <- !(probeNames(data.tmp) %in% all[,1])
      if(sum(missingAnnotation)>0) {
        pm(data.tmp)[missingAnnotation,] <- NA
        if(aType()=="PMMM") mm(data.tmp)[missingAnnotation,] <- NA
        lab_add <- " - white: other unannotated probe"
      }
    } 
    
    #set intensity of all probes in exprs, but not in pm (or if present mm) to 15 (all numbers chosen for the colors to work)
    exprs(data.tmp)[setdiff(1:dim(exprs(data.tmp))[1],c(rn_pm,rn_mm)),] <- 15
    
    #determine colors needed
    colors <- c("black","#222222","blue","lightblue","red")[c(TRUE,aType()=="PMMM",TRUE,aType()=="PMMM",TRUE)]
    #set up label
    if(aType()=="PMMM") {
      xlab <- paste("black/gray: regular match/mismatch probe \n blue/light blue: ", 
                    "control match/mismatch probe \n red: unannotated probe (control region)\n",lab_add,sep="")
    } else {
      xlab <- paste("black: regular probe \n blue: control probe \n red: unannotated probe (control region)\n",lab_add,sep="")
    }
    
    #plot the image
    RawDataReferenceArrayLayout <- tempfile(fileext = '.png')
    png(RawDataReferenceArrayLayout,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
    par(oma=c(17,2,0,3))
    image(data.tmp,col=colors,xlab=xlab,main="Array reference layout",cex.lab=0.8)
    dev.off()
    
    list(src = RawDataReferenceArrayLayout,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  }, deleteFile = T)
  
  # Positive/Negative controls
  output$posnegCOI<- renderImage({
    req(input$posnegCOI)
    
    RawDataPosNegPositions <- tempfile(fileext = '.png')
    png(file = RawDataPosNegPositions,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(oma=c(17,0,0,0),srt=90) 
    borderQC2(rawData())
    mtext("COI absolute values should be < 0.5\n", side=3, cex=0.7)
    dev.off()
    
    list(src = RawDataPosNegPositions,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })

  # Spatial images
  
  # PLM images
  
  # Nuse plot
  output$rawNuse <- renderImage({
    req(input$Nuse)
      
    rawdataNUSE <- tempfile(fileext = ".png")
    png(file = rawdataNUSE,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)   
    par(oma=c(17,0,0,0), cex.axis=1)   
    NUSE(rawData.pset(), col = plotColors(),main="Normalized Unscaled Standard Errors (NUSE)", axes = FALSE)
    if(length(sampleNames(rawData()))<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }  
    axis(1,at=1:length(sampleNames(rawData())),las=2,labels=sampleNames(rawData()), cex.axis=cexval)        
    axis(2, cex.axis=0.7)      
    if(length(levels(experimentFactor()))>1){ 
      legend("topright", levels(experimentFactor()), col = legendColors(),
             fill = legendColors(), cex = 0.7, bg = "white", bty = "o")
    }
    
    mtext("NUSE median value should be < 1.1\n", side=3, font=1, cex=0.7)
    dev.off()
    
    list(src = rawdataNUSE,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
    
  })
  
  # Rle plot
  output$rawRle <- renderImage({
    req(input$RLE)
    
    RawDataRLE <- tempfile(fileext = ".png")
    png(file = RawDataRLE,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)   
    par(oma=c(17,0,0,0), cex.axis=1)  
    Mbox(rawData.pset(), col = plotColors(), main = "Relative Log Expression (RLE)", axes = FALSE)
    if(length(sampleNames(rawData()))<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }  
    axis(1,at=1:length(sampleNames(rawData())),las=2,labels=sampleNames(rawData()), cex.axis=cexval)        
    axis(2, cex.axis=0.7)       	   
    if(length(levels(experimentFactor()))>1){ 
      legend("topright", levels(experimentFactor()), col = legendColors(), 
             fill = legendColors(), cex = 0.7, bg = "white", bty = "o")
    }              
    mtext("RLE distributions should be centered around 0\n", side=3, font=1, cex=0.7)
    dev.off()
    
    list(src = RawDataRLE,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
  })
  
  # Corelation plots
  output$correlRaw <- renderImage({
    req(input$correlRaw)
    clusterOption1 <- input$clusteroption1
    clusterOption2 <- input$clusteroption2 
    normMeth <- ""
    
    if(class(rawData()) == "AffyBatch") {
      Type <- "Raw"
      text1 <- "Raw data correlation plot"
    } else {
      if(normMeth == "") stop("When providing a normalized data object, the normMeth parameter is required")
      Type <- "Norm"
      text1 <- paste("Array correlation plot\nafter",normMeth,"normalization")
    } 
    if(length(sampleNames(rawData()))<2) {
      warning("Only one array in dataset, no correlation plot made")
    } else {
      rawdataCorrelation <- tempfile(fileext = ".png")
      png(file = rawdataCorrelation,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      if(length(sampleNames(rawData()))<MAXARRAY) {
        par(oma=c(17,0,0,0),cex.axis=0.7,cex.main=0.8)
        #subval <- 10
      } else {
        par(oma=c(17,0,0,0),srt=90,las=2,cex.axis=0.5,cex.main=0.8)
        #subval <- 16
      }        
      
      #note: for computing array correlation, euclidean would not make sense
      #only use euclidean distance to compute the similarity of the correlation vectors for the arrays
      COpt1 <- "pearson"
      if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
      crp <- cor(exprs(rawData()), use="complete.obs", method=COpt1)
      
      text1 <- paste(text1,"\ncorrelation method:",COpt1,"\ncluster method:",clusterOption2)
      
      switch(tolower(clusterOption1), 
             "pearson" = {
               my.dist <- function(x) cor.dist(x, abs=FALSE)
             },
             "spearman" = {
               my.dist <- function(x) spearman.dist(x, abs=FALSE)
             },
             "euclidean" = {
               my.dist <- function(x) euc(x)
             }
      )
      
      my.hclust <- function(d) hclust(d, method=clusterOption2)
      
      #in order to create some space to put colored symbols as well
      #sampleNames(rawData()) <- paste(sampleNames(rawData())," ")
      
      sideColors <- legendColors()[as.numeric(experimentFactor())]
      
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", symm=TRUE, density.info="density",
                main=text1, dendrogram="row", ColSideColors=sideColors)
      
      #correlationPlot(rawData())    
      #axis(1,side=3,at=seq(from=0.5, to=(length(sampleNames(rawData())))-0.5,by=1),
      #    labels=substr(as.character(sampleNames(rawData())),1,subval),las=2)
      #par(srt=0) 
      #plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE, 
      #    frame.plot = FALSE, xlim = c(0, 2), ylim = c(0,2))
      #text(1,1,text1,cex=1)  
      dev.off()
      list(src = rawdataCorrelation,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
    }
  })

  # PCA raw
  output$PCAraw <- renderImage({
    req(input$PCAraw)
    normMeth=""
    scaled_pca <- TRUE
    namesInPlot <- FALSE 
    
    # Scaled PCA by default
    if(length(sampleNames(rawData()))<3) {
      warning("Only",length(sampleNames(rawData())),"sample(s) in dataset, no PCA plot made")
    } else { 
      
      if(class(rawData()) == "AffyBatch"){
        #raw data
        tmain <- "PCA analysis of Raw data"
        Type="Raw"	
      } else{
        if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
        tmain <- paste("PCA analysis after", normMeth, "normalization")
        Type <- "Norm"
      }
      
      pca1 <- NULL  
      try(pca1 <- prcomp(t(exprs(rawData())[apply(exprs(rawData()),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
      if(is.null(pca1) & scaled_pca) {
        try(pca1 <- prcomp(t(exprs(rawData())[apply(exprs(rawData()),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=FALSE),TRUE)
        if(!is.null(pca1)) warning("pca with scaling unsuccessful, successfully retried without scaling")
      }
      if(!is.null(pca1)) {
        perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)
        
        cex.circle <- 1.5
        cex.text <- 0.7
        cex.legend <- 0.75
        tcol <- "#444444"
        
        rawdataPCAanalysis <- tempfile(fileext = ".png")
        png(file = rawdataPCAanalysis, width=WIDTH+200*(!namesInPlot), height=HEIGHT+283*(!namesInPlot),
            pointsize=POINTSIZE)
        
        if(!namesInPlot) {
          layout(rbind(c(1,1,2,2,5),c(3,3,4,4,5)))
        } else {
          layout(rbind(c(1,1,2,2),c(1,1,2,2),c(3,3,4,4),c(3,3,4,4)))
        }
        par(oma=c(20,0,5,0))
        plot(pca1$x[,1],pca1$x[,2],cex=cex.circle,pch=plotSymbols(),
             col=plotColors(),xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
             ylab=paste("PC2 (",perc_expl1[2],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,1],pca1$x[,2], sampleNames(rawData()),pos=4,cex=cex.text,col=tcol) 
        plot(pca1$x[,1],pca1$x[,3],cex=cex.circle,pch=plotSymbols(),
             col=plotColors(),xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
             ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,1],pca1$x[,3], sampleNames(rawData()),pos=4,cex=cex.text,col=tcol)
        plot(pca1$x[,2],pca1$x[,3],cex=cex.circle,pch=plotSymbols(),
             col=plotColors(),xlab=paste("PC2 (",perc_expl1[2],"%)",sep=""),
             ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,2],pca1$x[,3], sampleNames(rawData()),pos=4,cex=cex.text,col=tcol)
        barplot((100*pca1$sdev^2)/sum(pca1$sdev^2),xlab="components",ylab="% of total variance explained")
        
        if(namesInPlot) {
          if(length(levels(experimentFactor()))>1){ 
            legend("topright",levels(experimentFactor()),
                   pch=legendSymbols(),col=legendColors(),cex=cex.legend)
          }
        } else {
          par(mar=c(0,0,0,0))	
          plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
          if(length(levels(experimentFactor()))>1) {
            legend("topleft",c(levels(experimentFactor()),"",sampleNames(rawData())),
                   #             pch=c(rep(20,length(unique(experimentFactor()))+1),plotSymbols(),
                   pch=c(legendSymbols(),20,plotSymbols()),
                   col=c(legendColors(),"white",plotColors()),cex=(cex.legend+0.1)
                   #             ,fill=c(legendColors(),rep("white",length(experimentFactor())+1)),
                   #             border=c(legendColors(),rep("white",length(experimentFactor())+1))
            )
          } else {
            legend("topleft",sampleNames(rawData()),pch=plotSymbols(),
                   col=plotColors(),cex=0.7, bty = "n")
          }
        }
        
        mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
        dev.off()
        list(src = rawdataPCAanalysis,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
      }
    }
  })

  # Hierarchical clustering
  output$clusterRaw <- renderImage({
    req(input$clusterRaw)
    clusterOption1 <- input$clusteroption1
    clusterOption2 <- input$clusteroption2
    
    if(class(rawData())=="AffyBatch") {
      Type <- "Raw"
      main <- "Cluster dendrogram of raw data"
    } else {
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      Type <- "Norm"
      main <- paste("Cluster dendrogram of",normMeth,"normalized data")
    }  
    if(length(sampleNames(rawData()))<3) {
      warning("Only ",length(sampleNames(rawData()))," sample(s) in dataset, no clustering plot made")
    } else {  
      switch(tolower(clusterOption1), 
             "pearson" = {
               correl <- cor.dist(t(exprs(rawData())),abs=FALSE)
             },
             "spearman" = {
               correl <- spearman.dist(t(exprs(rawData())),abs=FALSE)
             },
             "euclidean" = {
               correl <- euc(t(exprs(rawData())))
             }
      )
      if(tolower(clusterOption2)!="ward.d2" & tolower(clusterOption2)!="ward.d") {
        clust <- hclust(correl, method = tolower(clusterOption2))
      } else {
        if(tolower(clusterOption2)=="ward.d2") {
          clust <- hclust(correl, method = "ward.D2")
        } else {
          clust <- hclust(correl, method = "ward.D")
        }
      }
      rawClustering <- tempfile(fileext = ".png")
      png(file = rawClustering,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      if(length(sampleNames(rawData()))<MAXARRAY) {
        cexval1 <- 0.75
        cexval2 <- 1.23
        cexval3 <- 0.55
      } else {
        cexval1 <- 0.55
        cexval2 <- 1.6
        cexval3 <- 0.41
      }   
      par(cex=cexval1,oma=c(14,1,0,0))	
      par(cex.axis=cexval2,cex.lab=cexval2,cex.main=cexval2)	       
      plot(clust, hang=-1, main=main, xlab=paste("distance:",clusterOption1), sub=paste(" cluster method:",clusterOption2))
      points(1:length(clust$order),rep(0,length(clust$order)),pch=15,col="white",cex=1.5)
      points(1:length(clust$order),rep(0,length(clust$order)),pch=plotSymbols()[clust$order],col=plotColors()[clust$order])
      if(length(levels(experimentFactor()))>1) { 
        legend("topright",levels(experimentFactor()),
               pch=legendSymbols(),col=legendColors())
      }
      par(cex=cexval3)    
      dev.off()
      
      list(src = rawClustering,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
    }
  })
  
  ##### PREPROCESSING #####
  # Upload CDF
  uploadcdfenv <- function(Data,cdf_filename,cdf_path){
    #spp <- c("Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu","Mm","Os","Rn",
    #         "Sc","Sp","Ss")
    #names(spp) <- c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
    #                "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
    #                "Drosophila melanogaster","Gallus gallus","Homo sapiens",
    #                "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
    #                "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa")
    
    #species <- spp[tolower(names(spp))==tolower(species)]
    
    #initial value
    CDFenv <- 0
    
    # recall which cdfName was added, in case no updated one is found (set back
    # even if it does not exist)
    presetCDF <- Data@cdfName
    
    #try without a version number
    suppressWarnings(try(CDFenv <- make.cdf.env(filename = cdf_filename,cdf.path=cdf_path,compress=FALSE,verbose = f),TRUE)) 
    
    if ((class(CDFenv)!="environment")) {
      Data@cdfName <- presetCDF
      warning("Could not load custom CDF environment for this chip type - object kept as is")
    }
    
    cat("current cdf environment loaded:",Data@cdfName,"\n")
    return(Data)
  }
  
  # Custom CD env
  addUpdatedCDFenv <- function(Data, species=NULL, type="ENSG") {
    # note: this function will add an updated cdf environment to the data object
    # and will overwrite a possible already loaded environment, unless no updated
    # cdf environment is found
    
    # developer's note: it may be of interest to find out whether available
    # species and types can be retrieved automatically from the brainarray website
    
    # Match the species to two letter symbols  
    spp <- c("Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu","Mm","Os","Rn",
             "Sc","Sp","Ss")
    names(spp) <- c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
                    "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
                    "Drosophila melanogaster","Gallus gallus","Homo sapiens",
                    "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
                    "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa")
    
    species <- spp[tolower(names(spp))==tolower(species)]

    
    #initial value
    CDFenv <- 0
    
    # recall which cdfName was added, in case no updated one is found (set back
    # even if it does not exist)
    presetCDF <- Data@cdfName
    
    #try to find updated cdf file of choice
    print(Data@cdfName<-paste(Data@annotation,species,type,sep="_"))
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    #try without a version number
    print(Data@cdfName<-paste(gsub("v[0-9]$","",Data@annotation),species,type,sep="_"))
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE)) 
    
    #if it hasn't loaded, try to download
    if ((class(CDFenv)!="environment")) {
      install.packages(tolower(paste(Data@annotation,species,type,"cdf",sep="")),
                       repos="http://brainarray.mbni.med.umich.edu/bioc")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    }
    
    #if it hasn't loaded, try to download without version number
    if ((class(CDFenv)!="environment")) {
      install.packages(tolower(paste(gsub("v[0-9]$","",Data@annotation),species,type,"cdf",sep="")),
                       repos="http://brainarray.mbni.med.umich.edu/bioc")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    }
    
    if ((class(CDFenv)!="environment")) {
      Data@cdfName <- presetCDF
      warning("Could not automatically retrieve CDF environment for this chip type - object kept as is")
    }
    
    cat("current cdf environment loaded:",Data@cdfName,"\n")
    return(Data)
  }
  
  # Normalize data
  normData <- eventReactive(eventExpr = input$preprocessing, {
    Data <- rawData()
    normMeth <- input$normMeth
    CDFtype <- input$CDFtype
    species <- input$species
    customCDF <- FALSE
    
    if(input$annotations=="Custom annotations"){
      customCDF <- TRUE
    } 
    
    if(input$perGroup == "group"){
      perGroup <- TRUE
    } else if (input$perGroup == "dataset"){
      perGroup <- FALSE
    }
    
    normMeth <- toupper(normMeth)
    #if customCDF option is chosen, apply to copy of Data, in order not to change the original data object
    Data.copy <- Data
    if(customCDF){
      print ("Change CDF before pre-processing")
      Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
    }
    print ("Pre-processing is running")
    
    nGroups <- 1
    if(perGroup) {
      nGroups <- max(1,length(levels(experimentFactor())))
      if(nGroups==1) warning("normalization per group requested, but no groups indicated in data set")
    }
    
    #if per group normalization required, or a method selected that does not return an ExpressionSet object,
    #make a model of class ExpressionSet to paste real values in, use the relatively fast RMA method
    #note that binding of ExpressionSet objects is NOT possible
    if((nGroups>1)) { # || (normMeth=="MAS5")) {
      normData <- rma(Data.copy)
      exprs(normData)[] <- NA
    }
    
    for(group in 1:nGroups) {
      if(nGroups==1) {
        Data.tmp <- Data.copy
      } else {
        Data.tmp <- Data.copy[,experimentFactor()==(levels(experimentFactor())[group])]
      }
      switch(normMeth, 
             "MAS5" = {
               #doesn't work
               normData.tmp <- mas5(Data.tmp) 
             },
             "GCRMA" = {
               if(customCDF) {
                 #probe library needed, first try whether this has been intalled, otherwise do so
                 probeLibrary <- tolower(paste(Data@annotation,species,CDFtype,"probe",sep=""))
                 loaded <- suppressWarnings(try(eval(parse("",-1,paste("library(",probeLibrary,")", sep=""))),TRUE))
                 if(class(loaded)=="try-error") {
                   install.packages(probeLibrary, repos="http://brainarray.mbni.med.umich.edu/bioc")
                 }
               }
               if(aType() == "PMMM") ntype = "fullmodel"
               if(aType() == "PMonly") ntype = "affinities" # good results if most of the genes are not expressed
               normData.tmp <- gcrma(Data.tmp, type=ntype, fast = FALSE)
             },
             "RMA" = {
               normData.tmp <- rma(Data.tmp)
             },
             "PLIER" = {
               if(aType() == "PMMM") ntype = "together"
               if(aType() == "PMonly") ntype = "pmonly"
               normData.tmp <- justPlier(Data.tmp, normalize=TRUE, norm.type = ntype)
             }
      )
      if(nGroups==1) {
        normData <- normData.tmp
        if(normMeth=="MAS5") exprs(normData)<-log2(exprs(normData))
      } else {
        try(
          if(normMeth=="MAS5"){
            exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- log2(exprs(normData.tmp))
          }else{
            exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- exprs(normData.tmp)
          },TRUE)
      }
      rm(normData.tmp, Data.tmp)
    }
    
    rm(Data.copy)
    return(normData)
  })
  
  # Normalized data table
  normDataTable <- eventReactive(eventExpr = input$preprocessing, {
    customCDF <- input$customCDF
    species <- input$species
    CDFtype <- input$CDFtype
    normData <- normData()
    
    if(is.null(customCDF)) stop("The customCDF parameter must be provided")
    if(customCDF) {
      if(species=="" || is.null(species)) {
        warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
        species <- deduceSpecies(rawData@annotation)
      }
      if(species=="" || is.null(species)){
        warning("Could not define species; the CDF will not be changed")
        customCDF<-FALSE
      }
    }
    
    #add column of IDs and normalized data to normDataTable
    normDataTable<-cbind(rownames(exprs(normData)),exprs(normData))
    
    #remove "_at" from custom probeset IDs to get to real ID
    if(customCDF) {
      control_rows <- grep("affx",tolower(normDataTable[,1]))
      if(length(control_rows)>0) {
        normDataTable[-control_rows,1]<-substring(normDataTable[-control_rows,1],1,nchar(normDataTable[-control_rows,1])-3)
      } else {
        normDataTable[,1]<-substring(normDataTable[,1],1,nchar(normDataTable[,1])-3)
      }
    }
    
    #add column names to normDataTable
    colnames(normDataTable)<-c(ifelse(customCDF,paste(CDFtype,"_ID",sep=""),"Probeset_ID"),colnames(exprs(normData)))
    
    #add gene name and description in case ensembl IDs have been used (otherwise there is no 1 to 1 mapping)
    if(customCDF && CDFtype=="ENSG") {
      #load gene name and description annotations
      library(biomaRt)
      
      spName <- ""
      if(species=="Ag" || species=="Anopheles gambiae") spName <- "agambiae"
      if(species=="At" || species=="Arabidopsis thaliana") spName <- "athaliana"
      if(species=="Bt" || species=="Bos taurus") spName <- "btaurus"
      if(species=="Ce" || species=="Caenorhabditis elegans") spName <- "celegans"
      if(species=="Cf" || species=="Canis familiaris") spName <- "cfamiliaris"
      if(species=="Dr" || species=="Danio rerio") spName <- "drerio"
      if(species=="Dm" || species=="Drosophila melanogaster") spName <- "dmelanogaster"
      if(species=="Gg" || species=="Gallus gallus") spName <- "ggallus"
      if(species=="Hs" || species=="Homo sapiens") spName <- "hsapiens"
      if(species=="MAmu" || species=="Macaca mulatta") spName <- "mmulatta"
      if(species=="Mm" || species=="Mus musculus") spName <- "mmusculus"
      if(species=="Os" || species=="Oryza sativa") spName <- "osativa"
      if(species=="Rn" || species=="Rattus norvegicus") spName <- "rnorvegicus"
      if(species=="Sc" || species=="Saccharomyces cerevisiae") spName <- "scerevisiae"
      if(species=="Sp" || species=="Schizosaccharomyces pombe") spName <- "spombe"
      if(species=="Ss" || species=="Sus scrofa") spName <- "sscrofa"
      
      try(ensembl <- useMart("ensembl", dataset = paste(spName,"_gene_ensembl",sep="")))
      if(exists("ensembl")) {
        try(annotationTable<-getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), mart=ensembl, uniqueRows=TRUE),TRUE)
      }
      if(exists("annotationTable")) {
        normDataTable <- as.data.frame(normDataTable,stringsAsFactors=FALSE)
        suppressWarnings(normDataTable<-cbind(normDataTable,annotationTable[match(normDataTable[,1],annotationTable[,1]),2:(dim(annotationTable)[2])]))
        normDataTable[,2:(dim(exprs(normData))[2]+1)] <- apply(normDataTable[,2:(dim(exprs(normData))[2]+1),drop=FALSE],2,as.numeric)
      } else {
        warning("No gene names and annotation could be retrieved from BioMart for this species or no connection could be established, gene information not added to normalized data table")
      }
    }
    
    return(normDataTable)
  })
  
  # Normalized data table output
  output$normdatatable <- renderTable({
    req(normDataTable())
    head(normDataTable(), 20)
  })
  
  # Normalized boxplot and density plots
  output$boxplotNorm <- renderImage({
    req(input$boxplotNorm)
    normMeth <- input$normMeth
    
    Type <- "Norm"
    tmain <- paste("Boxplot after ", normMeth, sep="")
    tmtext2 <- "Normalized log intensity\n\n\n"
    
    normBoxplot<- tempfile(fileext='.png')
    png(file = normBoxplot,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)  
    par(oma=c(17,0,0,0), cex.axis=1) 
    suppressWarnings(boxplot(normData(), col=plotColors() ,main=tmain, axes=FALSE, pch = 20, cex=0.7))
    if(length(levels(experimentFactor()))>1){ 
      legend("topright", levels(experimentFactor()),
             col=legendColors(),fill=legendColors(), cex = 0.7, bg = "white", bty = "o")
    }
    if(length(sampleNames(normData()))<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }  
    axis(1,at=1:length(sampleNames(normData())),las=2,labels=sampleNames(normData()), cex.axis=cexval)        
    axis(2, cex.axis=0.7)  
    mtext(tmtext2, side=2, cex=0.8)  	   
    mtext("Distributions should be comparable between arrays\n", side=3, font=1, 
          cex=0.7)
    dev.off()
    
    list(src = normBoxplot,
         width = WIDTH,
         height = HEIGHT,
         alt = "This is alternate text")
  })
  
  output$densityNorm <- renderImage({
    req(input$densityNorm)
    normMeth <- input$normMeth
    
    Type <- "Norm"
    
    normDensityHistogram<- tempfile(fileext='.png')
    png(file = normDensityHistogram,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    if(length(sampleNames(normData()))<MAXARRAY){
      cexval <- 0.65
      par(oma=c(12,0,0,0) )
    }else{
      cexval <- 0.45
      par(oma=c(0.1,0,0,0) )
    }    
    if(Type == "Raw"){
      hist(normData(), lwd=3, lt = 1:length(sampleNames(normData())), col = plotColors(), 
           which = "both", main="Density histogram of raw intensities", cex.axis = 0.7, cex.lab=0.8)
    }else{
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      hist(normData(), lwd=3, lt = 1:length(sampleNames(normData())), col = plotColors(), 
           main=paste("Density histogram after ", normMeth,"\n", sep=""), cex.axis = 0.7, cex.lab=0.8)
    } 
    
    legend("topright", substr(sampleNames(normData()),1,20), lwd=3, lt = 1:length(sampleNames(normData())),
           col = plotColors(), cex = cexval, bty = "n")           
    mtext( "Curves should be comparable between arrays\n", side=3, font=1, 
           cex=0.7)
    dev.off()
    list(src = normDensityHistogram,
         width = WIDTH,
         height = HEIGHT,
         alt = "This is alternate text")
  })
  
  output$correlNorm <- renderImage({
    req(input$correlNorm)

    clusterOption1 <- input$clusteroption1
    clusterOption2 <- input$clusteroption2 
    normMeth <- input$normMeth

    
    Type <- "Norm"
    text1 <- paste("Array correlation plot\nafter",normMeth,"normalization")
    
    if(length(sampleNames(normData()))<2) {
      warning("Only one array in dataset, no correlation plot made")
    } else {
      normdataCorrelation <- tempfile(fileext = ".png")
      png(file = normdataCorrelation,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      if(length(sampleNames(normData()))<MAXARRAY) {
        par(oma=c(17,0,0,0),cex.axis=0.7,cex.main=0.8)
        #subval <- 10
      } else {
        par(oma=c(17,0,0,0),srt=90,las=2,cex.axis=0.5,cex.main=0.8)
        #subval <- 16
      }        
      
      #note: for computing array correlation, euclidean would not make sense
      #only use euclidean distance to compute the similarity of the correlation vectors for the arrays
      COpt1 <- "pearson"
      if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
      crp <- cor(exprs(normData()), use="complete.obs", method=COpt1)
      
      text1 <- paste(text1,"\ncorrelation method:",COpt1,"\ncluster method:",clusterOption2)
      
      switch(tolower(clusterOption1), 
             "pearson" = {
               my.dist <- function(x) cor.dist(x, abs=FALSE)
             },
             "spearman" = {
               my.dist <- function(x) spearman.dist(x, abs=FALSE)
             },
             "euclidean" = {
               my.dist <- function(x) euc(x)
             }
      )
      
      my.hclust <- function(d) hclust(d, method=clusterOption2)
      
      sideColors <- legendColors()[as.numeric(experimentFactor())]
      
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", symm=TRUE, density.info="density",
                main=text1, dendrogram="row", ColSideColors=sideColors)
      
      dev.off()
      list(src = normdataCorrelation,
           width = WIDTH,
           height = HEIGHT,
           alt = "This is alternate text")
    }
  })
  
  output$PCAnorm <- renderImage({
    req(input$PCAnorm)
    normMeth <- input$normMeth
    scaled_pca <- TRUE
    namesInPlot <- FALSE 
    
    # Scaled PCA by default
    if(length(sampleNames(normData()))<3) {
      warning("Only",length(sampleNames(normData())),"sample(s) in dataset, no PCA plot made")
    } else { 
      
      tmain <- paste("PCA analysis after", normMeth, "normalization")
      Type <- "Norm"
      
      pca1 <- NULL  
      try(pca1 <- prcomp(t(exprs(normData())[apply(exprs(normData()),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
      if(is.null(pca1) & scaled_pca) {
        try(pca1 <- prcomp(t(exprs(normData())[apply(exprs(normData()),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=FALSE),TRUE)
        if(!is.null(pca1)) warning("pca with scaling unsuccessful, successfully retried without scaling")
      }
      if(!is.null(pca1)) {
        perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)
        
        cex.circle <- 1.5
        cex.text <- 0.7
        cex.legend <- 0.75
        tcol <- "#444444"
        
        normDataPCAanalysis <- tempfile(fileext = ".png")
        png(file = normDataPCAanalysis, width=WIDTH+200*(!namesInPlot), height=HEIGHT+283*(!namesInPlot),
            pointsize=POINTSIZE)
        
        if(!namesInPlot) {
          layout(rbind(c(1,1,2,2,5),c(3,3,4,4,5)))
        } else {
          layout(rbind(c(1,1,2,2),c(1,1,2,2),c(3,3,4,4),c(3,3,4,4)))
        }
        par(oma=c(20,0,5,0))
        plot(pca1$x[,1],pca1$x[,2],cex=cex.circle,pch=plotSymbols(),
             col=plotColors(),xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
             ylab=paste("PC2 (",perc_expl1[2],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,1],pca1$x[,2], sampleNames(normData()),pos=4,cex=cex.text,col=tcol) 
        plot(pca1$x[,1],pca1$x[,3],cex=cex.circle,pch=plotSymbols(),
             col=plotColors(),xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
             ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,1],pca1$x[,3], sampleNames(normData()),pos=4,cex=cex.text,col=tcol)
        plot(pca1$x[,2],pca1$x[,3],cex=cex.circle,pch=plotSymbols(),
             col=plotColors(),xlab=paste("PC2 (",perc_expl1[2],"%)",sep=""),
             ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,2],pca1$x[,3], sampleNames(normData()),pos=4,cex=cex.text,col=tcol)
        barplot((100*pca1$sdev^2)/sum(pca1$sdev^2),xlab="components",ylab="% of total variance explained")
        
        if(namesInPlot) {
          if(length(levels(experimentFactor()))>1){ 
            legend("topright",levels(experimentFactor()),
                   pch=legendSymbols(),col=legendColors(),cex=cex.legend)
          }
        } else {
          par(mar=c(0,0,0,0))	
          plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
          if(length(levels(experimentFactor()))>1) {
            legend("topleft",c(levels(experimentFactor()),"",sampleNames(normData())),
                   #             pch=c(rep(20,length(unique(experimentFactor()))+1),plotSymbols(),
                   pch=c(legendSymbols(),20,plotSymbols()),
                   col=c(legendColors(),"white",plotColors()),cex=(cex.legend+0.1)
                   #             ,fill=c(legendColors(),rep("white",length(experimentFactor())+1)),
                   #             border=c(legendColors(),rep("white",length(experimentFactor())+1))
            )
          } else {
            legend("topleft",sampleNames(normData()),pch=plotSymbols(),
                   col=plotColors(),cex=0.7, bty = "n")
          }
        }
        
        mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
        dev.off()
        list(src = normDataPCAanalysis,
             width = WIDTH,
             height = HEIGHT,
             alt = "This is alternate text")
      }
    }
  })
  
  output$clusterNorm <- renderImage({
    req(input$clusterNorm)
    
    normMeth <- input$normMeth
    clusterOption1 <- input$clusteroption1
    clusterOption2 <- input$clusteroption2

    Type <- "Norm"
    main <- paste("Cluster dendrogram of",normMeth,"normalized data")
    
    if(length(sampleNames(normData()))<3) {
      warning("Only ",length(sampleNames(normData()))," sample(s) in dataset, no clustering plot made")
    } else {  
      switch(tolower(clusterOption1), 
             "pearson" = {
               correl <- cor.dist(t(exprs(normData())),abs=FALSE)
             },
             "spearman" = {
               correl <- spearman.dist(t(exprs(normData())),abs=FALSE)
             },
             "euclidean" = {
               correl <- euc(t(exprs(normData())))
             }
      )
      if(tolower(clusterOption2)!="ward.d2" & tolower(clusterOption2)!="ward.d") {
        clust <- hclust(correl, method = tolower(clusterOption2))
      } else {
        if(tolower(clusterOption2)=="ward.d2") {
          clust <- hclust(correl, method = "ward.D2")
        } else {
          clust <- hclust(correl, method = "ward.D")
        }
      }
      normClustering <- tempfile(fileext = ".png")
      png(file = normClustering,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      if(length(sampleNames(normData()))<MAXARRAY) {
        cexval1 <- 0.75
        cexval2 <- 1.23
        cexval3 <- 0.55
      } else {
        cexval1 <- 0.55
        cexval2 <- 1.6
        cexval3 <- 0.41
      }   
      par(cex=cexval1,oma=c(14,1,0,0))	
      par(cex.axis=cexval2,cex.lab=cexval2,cex.main=cexval2)	       
      plot(clust, hang=-1, main=main, xlab=paste("distance:",clusterOption1), sub=paste(" cluster method:",clusterOption2))
      points(1:length(clust$order),rep(0,length(clust$order)),pch=15,col="white",cex=1.5)
      points(1:length(clust$order),rep(0,length(clust$order)),pch=plotSymbols()[clust$order],col=plotColors()[clust$order])
      if(length(levels(experimentFactor()))>1) { 
        legend("topright",levels(experimentFactor()),
               pch=legendSymbols(),col=legendColors())
      }
      par(cex=cexval3)    
      dev.off()
      
      list(src = normClustering,
           width = WIDTH,
           height = HEIGHT,
           alt = "This is alternate text")
    }
  })
  
  
  ##### DE ANALYSIS #####
  # Design matrix
  designMatrix <- reactive({
    validate(need(input$DESCFILE, message = "Description file required"))
    descdata <- desc()
    arrayNames <-descdata[,1]
    arrayNames <- make.names(arrayNames)
    
    experimentFactor <- experimentFactor()
    experimentFactor <- make.names(experimentFactor)
    experimentFactor <- as.factor(experimentFactor)
    
    design <- model.matrix(~ 0 + experimentFactor)
    rownames(design) <- arrayNames
    colnames(design) <- gsub("experimentFactor","",colnames(design))
    return(design)
  })
  
  # Initial contrast matrix
  m <- reactive({
    design <- designMatrix()
    matrix(cbind(c(1,-1),c(-1,1)), 
            ncol = length(colnames(design)),
            nrow = 2, 
            dimnames = list(NULL, gsub("experimentFactor","",colnames(design)))
      )
  })
  output$contrasts <- renderUI({
                  matrixInput(inputId = "contMatrix",
                              value = m(),
                              rows = list(extend = TRUE),
                              cols = list(names = TRUE)
                  )
  })

  # List of differential analysis tables
  diffTable <- eventReactive(eventExpr = input$startStat, {
    normDataTable <- normDataTable()
    description <- desc()
    
    genes <- normDataTable[,1]
    
    normDataTable <- apply(as.matrix(description[,2]),1,function(x) {
      as.numeric(normDataTable[,which(colnames(normDataTable)==x)])})
    normDataTable <- as.data.frame(normDataTable)
    colnames(normDataTable) <- description[,2]
    rownames(normDataTable) <- genes
    
    if(is.null(dim(normDataTable))) {
      stop("could not match array names from description file to normalized data
		file")
    }
    
    fit <- lmFit(normDataTable, design = designMatrix())
    contMatrix <- t(input$contMatrix)
    contMatrix <- contMatrix[,1:(ncol(contMatrix)-1)]
    colnames(contMatrix) <- paste("comparison",1:ncol(contMatrix))
    storage.mode(contMatrix) <- "numeric"
    fit <- contrasts.fit(fit, contMatrix)
    fit <- eBayes(fit)
    
    lapply(1:(ncol(t(input$contMatrix))-1), function(i){
      tab <- topTable(fit,number = Inf, coef = i)
    })
  })
 
  # Outputs for the stat page
  output$statOutput <- renderUI({
    l <- length(diffTable())
    tagList(
    # Stat tables
      lapply(diffTable(), function(i){
        output[[l]] <- renderTable({head(i,10)}, digits = 15, include.rownames = TRUE)
      }),
      
    # p-value 
      lapply(diffTable(), function(i){
        renderImage({
          pvalhistogram <- tempfile(fileext = ".png")
          adjPval <- i$adj.P.Val
          png(pvalhistogram,width=WIDTH,height=HEIGHT)
          hist (adjPval,main=paste("Adjusted P.Value Histogram"), xlab = "adj.p-values",col="blue")
          dev.off()
       
          list(src = pvalhistogram,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
        },outputArgs = list(width = "100%", height = "auto"))
      }),
    
    # logFC histograms
      lapply(diffTable(), function(i){
        renderImage({
          logFChistogram <- tempfile(fileext = ".png")
          logFCval <- i$logFC
          png(logFChistogram,width=1000,height=1000)
          hist (logFCval,main=paste("adapted fold change histogram"),xlab = "adapted fold changes",col="green", breaks=60)
          dev.off()
        
          list(src = logFChistogram,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
        },outputArgs = list(width = "100%", height = "auto"))
      })
    )
  })
  
  ##### GO ANALYSIS #####
  geneList <- eventReactive(eventExpr = input$goAnalysis, {
    lapply(diffTable(), function(i){
      statData <- rownames(i)
      statData <- cbind(statData,i)
    })
  })
  sig.genes <- eventReactive(eventExpr = input$goAnalysis, {
      lapply(geneList(), function(i){i[i$adj.P.Value <= input$Pvalcutoff,1]})
  })
    
  bkgd.genes <- eventReactive(eventExpr = input$goAnalysis, {
        lapply(geneList(), function(i){i[,1]})
  })
  

  output$goOutput <- renderUI({
    k <- length(sig.genes())
    tagList(
      mapply(function(i,j){
        ego <- clusterProfiler::enrichGO(
          gene     = i, 
          universe = j, 
          OrgDb    = org.Hs.eg.db,
          ont      = input$ont,
          pAdjustMethod = "fdr",
          pvalueCutoff = input$gopvalCutoff,
          readable = TRUE)
        
        output[[k]] <- renderText(head(ego))
        
        output[[k+1]] <- renderImage({barplot(ego, showCategory = 20)})
    
      }, sig.genes(),bkgd.genes())
    )
  })
  
  ##### PATHWAY ANALYSIS #####
  ##### NETWORK ANALYSIS #####
  ##### NAVIGATION #####
  observeEvent(input$next2, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Quality Control"
    )
  }) 
  observeEvent(input$next3, {
    updateTabsetPanel(session,
                      inputId = "QCsubset", 
                      selected = "Hybridization and overall signal quality"
    )
  })
  observeEvent(input$next4, {
    updateTabsetPanel(session,
                      inputId = "QCsubset", 
                      selected = "Signal comparability and bias diagnostic"
    )
  })
  observeEvent(input$next5, {
    updateTabsetPanel(session,
                      inputId = "QCsubset", 
                      selected = "Array correlation"
    )
  })
  observeEvent(input$next6, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Preprocessing and evaluation"
    )
  })
  observeEvent(input$next7, {
    updateTabsetPanel(session,
                      inputId = "Presubset", 
                      selected = "Signal comparability and bias of normalized intensities"
    )
  })
  observeEvent(input$next8, {
    updateTabsetPanel(session,
                      inputId = "Presubset", 
                      selected = "Normalized array correlation"
    )
  })
  observeEvent(input$next9, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Differential Analysis"
    )
  })
  observeEvent(input$next10, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "GO Enrichment"
    )
  })
  
  observeEvent(input$prev1, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Data Input"
    )
  }) 
  observeEvent(input$prev2, {
    updateTabsetPanel(session,
                      inputId = "QCsubset", 
                      selected = "Sample Quality"
    )
  }) 
  observeEvent(input$prev3, {
    updateTabsetPanel(session,
                      inputId = "QCsubset", 
                      selected = "Hybridization and overall signal quality"
    )
  }) 
  observeEvent(input$prev4, {
    updateTabsetPanel(session,
                      inputId = "QCsubset", 
                      selected = "Signal comparability and bias diagnostic"
    )
  }) 
  observeEvent(input$prev5, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Quality Control"
    )
  }) 
  observeEvent(input$prev6, {
    updateTabsetPanel(session,
                      inputId = "Presubset", 
                      selected = "Normalization method and annotation"
    )
  }) 
  observeEvent(input$prev7, {
    updateTabsetPanel(session,
                      inputId = "Presubset", 
                      selected = "Signal comparability and bias of normalized intensities"
    )
  }) 
  observeEvent(input$prev8, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Preprocessing and evaluation"
    )
  }) 
  observeEvent(input$prev9, {
    updateTabsetPanel(session,
                      inputId = "Mainset", 
                      selected = "Differential Analysis"
    )
  }) 
}
