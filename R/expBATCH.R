expBATCH <-function(D,batchCL,Conf=NA,mindim=1,maxdim=2,method="ppcca")
{
  if (missing(D)) {
    stop("Expression data are required for batch correction.\n")
  }
  if (missing(batchCL)) {
    stop("Provide variable batchCL which indicates batches.\n ")
  }
  if (nrow(D) != length(batchCL)) {
    stop("Expression data and batch variable should have the same number of samples.\n")
  }
  if (missing(mindim)) {
    mindim <- 1
  }
  if (missing(maxdim)) {
    maxdim <- 2
  }
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(mindim)) {
    stop("mindim should be a whole number.\n")
  }
  if (!is.wholenumber(maxdim)) {
    stop("maxdim should be a whole number.\n")
  }
  if (mindim > maxdim) {
    stop("mindim can not be greater than maxdim.\n")
  }
  if (maxdim > ncol(D)) {
    stop("maxdim can not be greater than the number of variables.\n")
  }
  if (maxdim > 10) {
    cat("Warning! Model fitting may become very slow for maxdim > 10.\n\n")
  }
  if(any(is.na(Conf))){Conf=rep(1,nrow(D))}
  type=names(table(Conf))
  nt=length(type)
  conF=as.matrix(Conf)
  ####################################
  ## installing missing R packages
  installed <- installed.packages()[,1]                              
  required <- c("mvtnorm", "mclust", "MetabolAnalyze", "sva", "stats", "ggplot2","RColorBrewer")
  toinstall <- required[!(required %in% installed)]
  if(length(toinstall) != 0){
    #install.packages(toinstall)
    source("https://bioconductor.org/biocLite.R")
    biocLite(toinstall)
  }
  lapply(required, require, character.only = TRUE)
  ## theme for ggplot
  theme <- theme(strip.text.y = element_text(),#rotate strip text to horizontal
                 axis.text = element_text(colour = "black", family="Arial", size=18),
                 axis.title.x = element_text(colour = "black", family="Arial",face="bold", size=20),
                 axis.title.y = element_text(colour = "black", family="Arial",face="bold", size=20),
                 axis.text.x = element_text(colour = "black", family="Arial", size=18),#remove y axis labels
                 legend.background = element_rect(fill = "white"),
                 panel.grid.major = element_line(colour = "white"),
                 legend.title = element_text(colour = "black", family="Arial",face="bold", size=20),
                 legend.text = element_text(colour="black", family="Arial", size = 18, face = "bold"),
                 strip.text.x = element_text(colour = "black", family="Arial", size = 18),
                 plot.title = element_text(size=20, family="Arial", face="bold"),
                 panel.border=element_rect(colour="black",size=2)) #format title
  ## colors
  grps<-names(table(batchCL))
  ngrps<-length(grps)
  jBrewColors <- brewer.pal(n = 8, name = "Dark2")
  if(ngrps<9)
  {
    cGsub<-jBrewColors[1:ngrps]
  }else{
    cGsub<-c(1:ngrps)
  }
  ## creating folder for storing results and setting the working directory
  dir.create("exploBATCHresults", showWarnings = FALSE)
  dir.create("exploBATCHresults/pcaBeforeCorrection", showWarnings = FALSE)
  dir.create("exploBATCHresults/ppccaBeforeCorrection", showWarnings = FALSE)
  dir.create("exploBATCHresults/findBATCH", showWarnings = FALSE)
  dir<-getwd()   
  setwd(paste0(dir,"/exploBATCHresults"))
  
  Ys <- as.matrix(D)
  ############## pca before batch correction
  rs1<-prcomp(Ys,scale=TRUE)  #pca
  pcaBeforeCorrection(rs1,grps,cGsub,batchCL,Conf,type,Ys,theme)
  ####
  
  ############## ppcca before batch correction
  if(nt==1)
  {
    designX<-as.matrix(model.matrix(~factor(batchCL))[,-1])
  }else{designX<-as.matrix(model.matrix(~factor(batchCL)+factor(conF))[,-1])}
  res1<-ppcca.metabol(Ys,designX,mindim,maxdim,scale="unit",plot.BIC=FALSE)      # run ppcca
  ppccaBeforeCorrection(res1,grps,cGsub,batchCL,Conf,type,Ys,theme)
  ####

  ############## detect batch and biological effect using ppcca
  print("Detecting batch effect via findBATCH....")
  res2<-tryCatch(ppcca.metabol.jack(Ys,batchCL,res1$q,res1$q,scale="unit"),error = function(e) NULL)        # run ppcca via Jacknifing
  rerun = FALSE
  if(nt!=1){
    rerun<-tryCatch(ppcca.metabol.jack(Ys,conF,res1$q,res1$q,scale="unit"),error = function(e) NULL)        # run ppcca via Jacknifing
  }
  LCI<-findBATCH(res2,res1,nt,rerun,theme)
  ####
  
  if(sum(LCI)>0){ ## correct batch effect if its significant
    dir.create("correctComBat", showWarnings = FALSE)
    ############## correct batch effect using combat
    comres<-correctComBat(Ys,batchCL,grps,cGsub,Conf,theme,type)

    ## assess combat correction
    dir.create("assessComBat", showWarnings = FALSE)
    print("Assessing CamBat corrected data via findBATCH....")
    rerun1<-tryCatch(ppcca.metabol.jack(t(comres),batchCL,res1$q,res1$q,scale="unit"),error = function(e) NULL)
    rerun12=FALSE
    if(nt!=1){
      rerun12<-tryCatch(ppcca.metabol.jack(t(comres),conF,res1$q,res1$q,scale="unit"),error = function(e) NULL)
    }
    assessComBat(rerun1,res1,nt,rerun12,theme)
    print("FINISHED ASSESSING COMBAT CORRECTED DATA!!")
    ####
    
    if(method=="ppcca")
    {
      dir.create("correctBATCH", showWarnings = FALSE)
      ############## correct batch effect using ppcca
      mindim=maxdim=69
      designX<-as.matrix(model.matrix(~factor(batchCL))[,-1])
      res11<-ppcca.metabol(Ys,designX,mindim,maxdim,scale="unit",plot.BIC=FALSE)     
      ppccaDat<-correctBatch(res11,designX,Ys,batchCL,ngrps,grps,Conf,nt,cGsub,comres,theme,type)
      print("correctBATCH FINISHED CORRECTING FOR BATCH EFFECT!!")
      
      ## assess PPCCA correction
      dir.create("assessCorrectBATCH", showWarnings = FALSE)
      print("Assessing correctBATCH corrected data via findBATCH....")
      rerun2<-tryCatch(ppcca.metabol.jack(ppccaDat,matrix(batchCL),res1$q,res1$q,scale="unit"),error = function(e) NULL)
      rerun22=FALSE
      if(nt!=1){
        rerun22<-tryCatch(ppcca.metabol.jack(ppccaDat,matrix(conF),res1$q,res1$q,scale="unit"),error = function(e) NULL)
      }
      assessCorrectBATCH(rerun2,res1,nt,rerun22,theme)
      print("FINISHED ASSESSING correctBATCH CORRECTED DATA!!")
      ####
    }
  }
  sessionInfo()
  save.image(paste0(Sys.Date(),"_exploBATCH_results.RData"))
  print("exploBATCH completed running!! Output in your working directory.")
}
