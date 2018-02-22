ppccaBeforeCorrection <-function(res1,grps,cGsub,batchCL,Conf,type,Ys,theme)
{
  #### number of ppcs
  x <- c(1:length(res1$BIC))
  y <- res1$BIC
  bic <- bicp <- data.frame(x=x,y=y)
  colnames(bicp)<-c("Number-pPCs","BIC")
  rownames(bicp)<-paste0("pPC-",1:length(res1$BIC))
  write.table(bicp,paste0("ppccaBeforeCorrection/",Sys.Date(),"_NumberOfpPCs.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig2 <- ggplot(data=bic, aes(x, y)) + geom_point(size=3)+geom_line()+labs(title=" ",x="Number of pPCs", 
          y="BIC")+theme_bw()+theme+geom_vline(xintercept=res1$q,lty=5,col=2,size=1.5)
  ggsave(filename ="ppccaBeforeCorrection/NumberOfpPCs.pdf", device=cairo_pdf,plot=fig2)
  #### ppcca scoresplot before batch correction 
  x <- res1$scores[,1]
  y <- res1$scores[,2]
  sdat2 <- data.frame(x=x,y=y,group=grps[batchCL])
  tdat2 <- data.frame(sdat2[,1:2],grps[batchCL],Conf)
  colnames(tdat2)<-c("pPC-1","pPC-2","Batch","Type")
  rownames(tdat2)<-rownames(Ys) 
  write.table(tdat2,paste0("ppccaBeforeCorrection/",Sys.Date(),"ppccaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig3<-ggplot(data= sdat2, aes(x, y,colour=cGsub[batchCL])) +
  geom_point(aes(shape=factor(c(16,15)[Conf],labels=type)),size=8,alpha=.7)+
  scale_shape_manual(values=c(16,15),guide=guide_legend(override.aes=aes(size=6)))+
  labs(title="",x="pPC-1",y="pPC-2") +theme_bw()+theme+
  scale_colour_manual("Batch",labels=grps,values=cGsub)+theme(legend.key = element_blank())+
  theme(legend.title=element_blank())+guides(colour = guide_legend(override.aes = list(linetype=0)))
  ggsave(filename ="ppccaBeforeCorrection/ppccaBeforeCorrection.pdf", device=cairo_pdf,plot = fig3)
  #### pairs plot for ppcca
  pdat<-res1$scores
  colnames(pdat)<-paste0("pPC-",1:res1$q)
  write.table(pdat,paste0("ppccaBeforeCorrection/",Sys.Date(),"_pairsppccaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pdf(file="ppccaBeforeCorrection/pairsplotppccaBeforeCorrection.pdf",onefile=TRUE)
  pairs(pdat[,1:res1$q],col= cGsub[batchCL],pch=c(16,15)[Conf])
  dev.off()
}
