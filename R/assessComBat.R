assessComBat <-function(rerun1,res1,nt,rerun12,theme)
{
  if(nt!=1){
    #### Detect biological effect
    blci<-data.frame(matrix(paste0("pPC-",1:res1$q)),matrix(rerun12$coefficients[,2]),rerun12$coeffCI[,3:4])
    colnames(blci)<-c("x","y","ylo","yhi")
    bldat<-blci[,-1]
    rownames(bldat)<-blci[,1]
    colnames(bldat)<-c("Effect","LowerCI","UpperCI")
    write.table(bldat,paste0("assessComBat/",Sys.Date(),"BioEFFECT.txt"),quote=TRUE,sep="\t",row.names=TRUE)
    bpl <- ggplot(blci, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+labs(title="",x="",y="")+
    theme_bw()+theme+coord_flip()+geom_hline(aes(x=0,yintercept=0),lty=5,col=2,size=1.5)
    ggsave(filename="assessComBat/BioEFFECT.pdf", device=cairo_pdf,plot = bpl)
  }
  #### Detect batch effect
  lci<-data.frame(matrix(paste0("pPC-",1:res1$q)),matrix(rerun1$coefficients[,2]),rerun1$coeffCI[,3:4])
  colnames(lci)<-c("x","y","ylo","yhi")
  ldat<-lci[,-1]
  rownames(ldat)<-lci[,1]
  colnames(ldat)<-c("Effect","LowerCI","UpperCI")
  write.table(ldat,paste0("assessComBat/",Sys.Date(),"batchEffect.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pl<-ggplot(lci, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+labs(title="",x="",y="")+
    theme_bw()+theme+coord_flip()+geom_hline(aes(x=0,yintercept=0),lty=5,col=2,size=1.5)
  ggsave(filename="assessComBat/BatchEffect.pdf", device=cairo_pdf,plot = pl)
}
