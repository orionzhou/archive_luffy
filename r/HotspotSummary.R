##file: location of .sum files to input
### if suf=F DONOT INCLUDE .sum prefix
### Can search for multiple files; can include "*" usage for finding files
##:LR: threshold for finding hotspot
##LRm lower threshold for hotspot region

HotspotSummary=function(file,LR=12,LRm=4,plot=F,output=T,ofile="temp.txt",title=file) {
    temp=system(paste("ls ",file,sep=""),intern=T)
    nf=length(temp)
    ##READ IN DATA
    Res=NULL
    if(nf==0){
      cat("No files to read-- Used Command \n\n",paste("ls ",file,sep=""),"\n")
      return()
    }
    for(i in 1:nf){
      cat("Reading information from ",temp[i],"\n")
      res=read.table(temp[i],skip=1,fill=T)
      for(ii in 2:dim(res)[1]){
        for(j in 1:dim(res)[2])
        if(is.na(res[ii,j])) res[ii,j]<-res[ii-1,j]
      }
      Res=rbind(Res,res)
    }
    ##ORDER
    index=order(Res[,1])
    Res=Res[index,1:3]
    LRst=Res[,3]
    Regions=NULL
    Hotspots=NULL
    while(max(LRst)>LR){
      jj=which.max(LRst)
      Hotspots=rbind(Hotspots,Res[jj,])
      j=jj
      FLAG=1
      while(j>1 && FLAG) {
        if(LRst[j-1]>LRm) {
          j=j-1
        } else {
          if(j>2) {
            index=1:(j-2)
            index=index[Res[index,2]>Res[j-1,1]] ##regions to left that overlap
            if(length(index)>0) {
              if(max(LRst[index])>LR) { ##does region overlap with 
                j=j-1
              }else{
                FLAG=0 ##j-1 is not in hotspot region
              }
            }
          }else{
            FLAG=0
          }
        }
      }
      jl=j
      j=jj
      FLAG=1
      while(j<length(LRst) && FLAG){
        if(LRst[j+1]>LRm){
          j=j+1
        }else{
          if(j+1<length(LRst)){
            index=(j+2):length(LRst)
            index=index[Res[index,1]<Res[j+1,2]] ##regions to right that overlap
            if(length(index)>0){
              if(max(LRst[index])>LR){ ##does region overlap with 
                j=j+1
              }else{
                FLAG=0 ##j+1 is not in hotspot region
              }
            }
          }else{
            FLAG=0
          }
        }
      }
      jr=j
      Regions=rbind(Regions,c(Res[jl,1],Res[jr,2]))
      LRst[jl:jr]=0
                       
    }
    ###Output summary
    nh=dim(Hotspots)[1]
    if(output) sink(ofile)
    if(length(nh)==0 || nh==0){
      cat("No hotpots\n")
    }else{
      index=order(Hotspots[,1])
      Hotspots=Hotspots[index,]
      Regions=Regions[index,]
      dim(Regions)=c(nh,2)
      cat("Hotspot Start\t Hotspot Stop\t LR \t Region Start\t Region Stop\n")
      for(i in 1:nh){
        cat(Hotspots[i,1],"\t",Hotspots[i,2],"\t",Hotspots[i,3],"\t",Regions[i,1],"\t",Regions[i,2],"\n")
      }
    }
    if(output) sink()
    if(plot){
      m=dim(Res)[1]
      plot(c(Res[1,1],Res[m,2])/1e3,c(0,max(Res[,3])),type="n",xlab="Position (kb)",ylab="LR",main=title)
      for(i in 1:m){
        if(Res[i,3]>0){
          lines(Res[i,1:2]/1e3,c(Res[i,3],Res[i,3]),lwd=2)
        }
      }
    if(length(nh)!=0 && nh!=0){
      for(i in 1:nh){
        lines(c(Hotspots[i,1],Hotspots[i,1])/1e3,c(Hotspots[i,3],0),lwd=3,col=2)
        lines(c(Hotspots[i,2],Hotspots[i,2])/1e3,c(Hotspots[i,3],0),lwd=3,col=2)
        lines(Regions[i,]/1e3,c(0,0),col=3,lwd=3)
      }
    }
      abline(h=LR,col=2)
      
    }
}
