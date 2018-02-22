correctComBat <-function(Ys,batchCL,grps,cGsub,Conf,theme,type)
{
  comdat<-t(Ys); 
  comres<-ComBat(dat=comdat,batch=batchCL,mod=NULL, par.prior=TRUE, prior.plots=FALSE) # run combat
  write.table(t(comres),paste0("correctComBat/",Sys.Date(),"combatCorrectedData.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  ## compare combat predicted data with observed data
  x <- c(comres)
  y <- c(comdat)
  tdat4<-sdat4<-data.frame(x,y)
  colnames(tdat4)<-c("Combat_predictedData","ObservedData")
  write.table(tdat4,paste0("correctComBat/",Sys.Date(),"Combat_predicted_vs_ObservedData.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig5<-ggplot(data=sdat4, aes(x, y,colour=c("orange")[rep(1,length(y))])) +
    geom_point(aes(shape=factor(16)),size=3,alpha=.2)+scale_colour_manual("",values="orange")+theme(legend.key = element_blank())+
    labs(title=paste0("Correlation: ",round(cor(x,y),2)),x="Combat predicted data",y="Observed data") +theme_bw()+theme+
    theme(legend.position = c(10000, 10000))+theme(legend.key = element_blank())+theme(legend.title=element_blank())+
    guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="correctComBat/Combat_predicted_vs_ObservedData.pdf", device=cairo_pdf,plot = fig5)
  ### pca on combat batch corrected
  res5<-prcomp(t(comres),scale=TRUE)
  x <- res5$x[,1]
  y <- res5$x[,2]
  sdat6 <- data.frame(x=x,y=y,group=grps[batchCL])
  tdat6<-data.frame(sdat6[,1:2],grps[batchCL],Conf)
  colnames(tdat6)<-c("PC-1","PC-2","Batch","Type")
  rownames(tdat6)<-rownames(Ys)
  write.table(tdat6,paste0("correctComBat/",Sys.Date(),"pca_After_Combat_Correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig7<-ggplot(data= sdat6, aes(x, y,colour=cGsub[batchCL])) + 
    geom_point(aes(shape=factor(c(16,15)[Conf],labels=type)),size=8,alpha=.7)+
    scale_shape_manual(values=c(16,15),guide=guide_legend(override.aes=aes(size=6)))+
    labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
    scale_colour_manual("Batch",labels=grps,values=cGsub)+theme(legend.key = element_blank())+
    theme(legend.title=element_blank())+guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="correctComBat/pca_After_Combat_Correction.pdf", device=cairo_pdf,plot = fig7)
  comres
}
