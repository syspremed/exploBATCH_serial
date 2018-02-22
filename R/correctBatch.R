correctBatch <-function(res11,designX,Ys,batchCL,ngrps,grps,Conf,nt,cGsub,comres,theme,type)
{
  XD <- standardize(designX)
  Covars <- cbind(rep(1,nrow(Ys)),XD)
  if(ngrps==2)
  {
    resid<-( t(res11$scores)-matrix(res11$coefficients[,2])%*%t(matrix(Covars[,2])) )
  }else{
    resid<-( t(res11$scores)-res11$coefficients[,2:ngrps]%*%t(Covars[,2:ngrps]) )
  }
  mumut <- matrix(colMeans(Ys),nrow=ncol(Ys),ncol=nrow(Ys),byrow=FALSE)
  simData <- array(NA,c(ncol(Ys),nrow(Ys),100))
  for(i in 1:100)
  {
    GuassianSpray <- t(rmvnorm( nrow(Ys), rep(0,ncol(Ys)), res11$sig*diag(ncol(Ys)) ))
    simData[,,i] <- (res11$loadings%*%resid + mumut + GuassianSpray)
  }
  ppccaDat <- apply(simData,c(1,2),median)
  ppccaDat <- t(ppccaDat)  # batch corrected expression data
  rownames(ppccaDat)<-rownames(Ys)
  write.table(ppccaDat,paste0("correctBatch/",Sys.Date(),"ppccaCorrectedData.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  ############## compare correctBATCH data and observed data
  x <- c(ppccaDat)
  y <- c(Ys)
  sdat3<-tdat3<-data.frame(x,y)
  colnames(tdat3)<-c("PPCCA_predictedData","ObservedData")
  write.table(tdat3,paste0("correctBatch/",Sys.Date(),"PPCCA_predicted_vs_ObservedData.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig4<-ggplot(data= sdat3, aes(x, y,colour=c("orange")[rep(1,length(y))])) +
    geom_point(aes(shape=factor(16)),size=3,alpha=.2)+
    scale_colour_manual("",values="orange")+theme(legend.key = element_blank())+
    labs(title=paste0("Correlation: ",round(cor(x,y),2)),x="PPCCA predicted data",y="Observed data") +theme_bw()+theme+
    theme(legend.position = c(10000, 10000))+theme(legend.key = element_blank()) +theme(legend.title=element_blank())+
    guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="correctBatch/PPCCA_predicted_vs_ObservedData.pdf", device=cairo_pdf,plot = fig4)
  ### pca on ppcca batch corrected data
  res4<-prcomp(ppccaDat, scale=TRUE )
  x <- res4$x[,1]
  y <- res4$x[,2]
  sdat5 <- data.frame(x=x,y=y,group=grps[batchCL])
  tdat5<-data.frame(sdat5[,1:2],grps[batchCL],Conf)
  colnames(tdat5)<-c("PC-1","PC-2","Batch","Type")
  rownames(tdat5)<-rownames(Ys)
  write.table(tdat5,paste0("correctBatch/",Sys.Date(),"pca_After_PPCCA_correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig6<-ggplot(data=sdat5, aes(x, y,colour=cGsub[batchCL])) + 
    geom_point(aes(shape=factor(c(16,15)[Conf],labels=type)),size=8,alpha=.7)+
    scale_shape_manual(values=c(16,15),guide=guide_legend(override.aes=aes(size=6)))+
    labs(title="",x="PC-1",y="PC-2")+theme_bw()+theme+
    scale_colour_manual("Batch",labels=grps,values=cGsub)+theme(legend.key = element_blank())+
    theme(legend.title=element_blank())+guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="correctBatch/pca_After_PPCCA_correction.pdf", device=cairo_pdf,plot = fig6)
  ############## ppcca vs combat correction
  x <- c(t(comres))
  y <- c(ppccaDat)
  sdat7 <- tdat7<-data.frame(x=x,y=y)
  colnames(tdat7)<-c("CombatPredictedData","PPCCApredictedData")
  write.table(tdat7,paste0("correctBatch/",Sys.Date(),"Combat_vs_PPCCA_predictedData"),quote=TRUE,sep="\t",row.names=TRUE)
  fig8<-ggplot(data=sdat7, aes(x, y,colour=c("orange")[rep(1,length(y))])) + 
    geom_point(aes(shape=factor(16)),size=3,alpha=.2)+scale_colour_manual("",values="orange")+theme(legend.key = element_blank())+
    labs(title=paste0("Correlation: ",round(cor(x,y),2)),x="Combat-predicted data",y="PPCCA-predicted data") +theme_bw()+theme+
    theme(legend.position = c(10000, 10000))+theme(legend.key = element_blank())+theme(legend.title=element_blank())+
    guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="correctBatch/Combat_vs_PPCCA_predictedData.pdf", device=cairo_pdf,plot = fig8)
  ppccaDat
}
