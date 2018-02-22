pcaBeforeCorrection <-function(rs1,grps,cGsub,batchCL,Conf,type,Ys,theme)
{
  #### proportion of variation
  variation <- c(rs1$sdev^2/sum(rs1$sdev^2))
  x <- 1:length(variation)
  y <- variation
  variat <- data.frame(x=x, y=y)
  colnames(variat)<-c("numberPCs","PropVariation")
  rownames(variat)<-paste0("PC-",1:length(variation))
  write.table(variat,paste0("pcaBeforeCorrection/",Sys.Date(),"_ProportionVariationPCAbeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  prop <- ggplot(data=variat, aes(x, y)) + geom_point(size=3)+geom_line()+
    labs(title=" ",x="Number of PCs",y="Proportion of variation")+theme_bw()+theme
  ggsave(filename ="pcaBeforeCorrection/ProportionVariationPCAbeforeCorrection.pdf", device=cairo_pdf,plot = prop)
  #### pca scores
  x <- rs1$x[,1]
  y <- rs1$x[,2]
  sdat1 <- data.frame(x=x,y=y,group=grps[batchCL])
  tdat1<-data.frame(sdat1[,1:2],grps[batchCL],Conf)
  colnames(tdat1)<-c("PC-1","PC-2","Batch","Type")
  rownames(tdat1)<-rownames(Ys)
  write.table(tdat1,paste0("pcaBeforeCorrection/",Sys.Date(),"_pcaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig1<-ggplot(data= sdat1, aes(x, y,colour=cGsub[batchCL])) +
    geom_point(aes(shape=factor(c(16,15)[Conf],labels=type)),size=8,alpha=.7)+
    scale_shape_manual(values=c(16,15),guide=guide_legend(override.aes=aes(size=6)))+
    labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
    scale_colour_manual("Batch",labels=grps,values=cGsub)+theme(legend.key = element_blank())+
    theme(legend.title=element_blank())+guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="pcaBeforeCorrection/pcaBeforeCorrection.pdf", device=cairo_pdf,plot = fig1)
  #### pairs plot for pca
  pd<-rs1$x
  colnames(pd)<-paste0("PC-",1:ncol(rs1$x))
  write.table(pd,paste0("pcaBeforeCorrection/",Sys.Date(),"pairspcaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pdf(file="pcaBeforeCorrection/pairsplotpcaBeforeCorrection.pdf",onefile=TRUE)
  pairs(pd[,1:10],col=cGsub[batchCL],pch=c(15,15)[Conf])
  dev.off()
}
