## function to read the data
readData=function(dirPath, ValAdjust=NULL){
  pnames=read.csv("scripts/paramNames.csv", header=T, stringsAsFactors=F)[ ,1]
  pnames=gsub(" ", "", pnames)
  FileIn=dir(dirPath)
 # FilesHer2=dir(her2dist)
  #tempSort=read.csv(StainFile, header=T, stringsAsFactors=FALSE)
  FilesInS=unlist(regmatches(FileIn, regexec("[0-9]+_[a-z]+[1-5]*_[0-5]", FileIn)))
  FilesInB=substr(FileIn, 1, 4)
  FilesInR=unlist(regmatches(FileIn, regexec("[a-z]+[1-5]*_[0-5]", FileIn)))
  ### Initialise storage data frames or something
  DSum=data.frame(matrix(NA, round(length(FileIn)*4), ncol=19), stringsAsFactors=FALSE)
  DSumName=rep(NA, round(length(FileIn)*4))
  ESum=data.frame(matrix(NA, length(FilesInS), ncol=73))
  VSum=data.frame(matrix(NA, length(FilesInS), ncol=73))
  SampleStorage= vector(mode = "list", length = length(unique(FilesInB)))
  names(SampleStorage)=unique(FilesInB)
  for(i in names(SampleStorage)){
    x1=unique(FilesInR[FilesInB==i])
    SampleStorage[[i]]=vector(mode="list", length=length(x1))
    names(SampleStorage[[i]])=x1  
  }
  ncount=0
   # check whether all files are available for a certain sample, otherwise print out an error message
  # Need to check CellLabels, presence of the Spot List, and that it is the right length
  #:length(FilesInS)
  for (i in 1:length(FilesInS)){
    
    tempA=read.table(sprintf('%s%s', dirPath, FileIn[i]), header=T)
  #  tempB=read.table(sprintf('%s%s', her2dist, FilesHer2[i]), header=T)
    colnames(tempA)=gsub("NuclearAdjustedInt", "NucAdjustedInt", colnames(tempA))
    colnames(tempA)=gsub("memb_temp", "memb_BackAdjustedInt", colnames(tempA))
    ## DOUBLE CHECK WE CAN REMOVE! xb=grep("^NucAdjustedInt", colnames(tempA))
    ##if (length(xb)>0){
    ##  tempA=tempA[ , -xb]
    ##}
 
    ## ERROR MSG: CELL LABELS
    if (length(grep("CellLabels", colnames(tempA)))==0){
      print(sprintf('Missing Cell Labels: sample %s', FilesInS[i]))
      next
    }   
    # order the tempA samples to be in dapi, memb, nucl, cent, Her2
    
    tidx=unique(c(grep("rownames", colnames(tempA)),
            grep("Loc", colnames(tempA)),
           grep("NN", colnames(tempA)),
                 grep("dapi_", colnames(tempA)),
                 grep("memb_", colnames(tempA)),
           grep("nucl_", colnames(tempA)),
                grep("cent_", colnames(tempA)),
                grep("Her2_", colnames(tempA)),
           grep("CellLabels", colnames(tempA))))
    tempA=tempA[ ,tidx]
    
   # browser()
 #   tempA=cbind(tempA, tempB[ ,grep("cent_Dist", colnames(tempB))])
#    c3=cbind(tempA, tempB[ ,grep("cent_", colnames(tempB))])
    
    if (ncol(tempA)!=length(pnames)){
    #  browser()
      dv=setdiff( pnames, colnames(tempA))
      tempA[ , (ncol(tempA)+1):(ncol(tempA)+length(dv))]=NA
      colnames(tempA)[(ncol(tempA)-length(dv)+1):(ncol(tempA))]=dv
     }
    
    x1=tempA[ ,grep("SpotDist", colnames(tempA))]/tempA$dapi_MajorAxisLength
    colnames(x1)=paste(colnames(x1), "Norm", sep="") 
    tempA=cbind(tempA, x1)
     
    tempA$cent_ConvexAreaNorm=tempA$cent_ConvexArea/tempA$dapi_Area
    tempA$Her2_ConvexAreaNorm=tempA$Her2_ConvexArea/tempA$dapi_Area
    tempA$Her2_SpotAreaRange=tempA$Her2_SpotAreaMax-tempA$Her2_SpotAreaMin
    tempA$cent_SpotAreaRange=tempA$cent_SpotAreaMax-tempA$cent_SpotAreaMin
    
    tempA$memb_pos=(tempA$memb_MembraneIntensity-tempA$memb_MeanIntensity)/
      tempA$memb_MembraneIntensity
    tempA$memb_PixVar=(tempA$memb_OverLapIntensity-tempA$memb_nonLapIntensity)/
      tempA$memb_OverLapIntensity
    tempA$memb_PixVar[which(tempA$memb_PixVar=="-Inf")]=0
    ## Adjust the values if specified
    
    if (!is.null(ValAdjust)){
      tempA=DataTidy(tempA, ValAdjust)
    }
   # browser()
    tempA$Her2_centRatio=tempA$Her2_Area/tempA$cent_Area
    #browser()
    ## Add the information to a new list:
    SampleStorage[[FilesInB[i]]][[FilesInR[i]]]=tempA
    ## Add to a summary data.frame
    tidx=grep('Tumour', tempA$CellLabels)
    fidx=c(grep('Fib', tempA$CellLabels), grep('Lymph', tempA$CellLabels))
    didx=grep('DCIS', tempA$CellLabels)
  #  oidx=grep("Other", tempA$CellLabels)
    colidx=c(grep("_area$", colnames(tempA), ignore.case=TRUE),
             grep("BackAdjustedInt", colnames(tempA)),
             grep("BackgroundAdjust", colnames(tempA)),
             grep("NucAdjustedInt", colnames(tempA)),
             grep("PixelCV", colnames(tempA)),
             grep("Ratio", colnames(tempA)))
    tt=colnames(tempA)[colidx]
    if (length(tt) != 15) { 
      stop(sprintf("Infile %s filtered to %i cols instead of 14", FilesInS[i], length(tt)))
    }
    colidx=colidx[order(tt)]
#browser()
    if (length(tidx)>0){
      ncount=ncount+1
      DSum[ncount, 1:15]=colMeans(tempA[tidx, colidx], na.rm=T)
      DSum[ncount, 16]=length(tidx)
      DSum[ncount, 17]="Tumour"
      DSum[ncount, 18]=FilesInS[i]
      DSum[ncount, 19]=substr(FilesInS[i], 1,4)
      DSumName[ncount]=sprintf('%s_Tum', FilesInS[i])
    }
#     if (length(oidx)>0){
#       ncount=ncount+1
#       DSum[ncount, 1:15]=colMeans(tempA[oidx, colidx], na.rm=T)
#       DSum[ncount, 16]=length(oidx)
#       DSum[ncount, 17]="Other"
#       DSum[ncount, 18]=FilesInS[i]
#       DSum[ncount, 19]=substr(FilesInS[i], 1,4)
#       DSumName[ncount]=sprintf('%s_Other', FilesInS[i])
#     }
    if (length(fidx)>0){
      ncount=ncount+1
    DSum[ncount, 1:15]=colMeans(tempA[fidx, colidx], na.rm=T)
      DSum[ncount, 16]=length(fidx)
      DSum[ncount, 17]="Fib"
      DSum[ncount, 18]=FilesInS[i]
      DSum[ncount, 19]=substr(FilesInS[i], 1,4)
      DSumName[ncount]=sprintf('%s_Fib', FilesInS[i])
    }
    if (length(didx)>0){
      ncount=ncount+1
      DSum[ncount, 1:15]=colMeans(tempA[didx, colidx], na.rm=T)
      DSum[ncount, 16]=length(didx)
      DSum[ncount, 17]="DCIS"
      DSum[ncount, 18]=FilesInS[i]
      DSum[ncount, 19]=substr(FilesInS[i], 1,4)
      DSumName[ncount]=sprintf('%s_DCIS', FilesInS[i])
    }
    a1=order(colnames(tempA))
    tempA=tempA[ ,a1]
    ESum[ i, ]=colMeans(tempA[tidx,-c(1,75:77)], na.rm=T)
    VSum[ i, ]=sapply(2:74,function(x) sd(tempA[tidx,x], na.rm=T))
  }
  x1=!is.na(DSumName)
  DSum=DSum[!is.na(DSumName), ]
  rownames(DSum)=DSumName[!is.na(DSumName)]
  colnames(DSum)=c(sort(tt), "Ncells", "CellType", "ID", "Pat")
  colnames(ESum)=colnames(tempA)[c(2:74)]
  rownames(ESum)=FilesInS; rownames(VSum)=FilesInS
  colnames(VSum)=paste("Var", colnames(tempA)[c(2:74)], sep="_")

#browser()
  Output=list(Summary=data.frame(DSum), ImSample=SampleStorage, Ext.Sum=data.frame(ESum),
              Var.Sum=data.frame(VSum)) #, SpotSamples=SpotStorage)
  return(Output)
}

DataTidy=function(data, SArea){
  # TIDY UP DATA BY REMOVING NEGATIVE INTENSITIES
  # Also, shift the spot size distribution so all cells have at least 1 spot.
  # SArea = estimated spot area
  
  ## ADJUST VALUES TO 0
  less0=which(data<0, arr.ind=T)
  data[less0]=0
  
  ## SHUFFLE spot count by 1
  ## Shuffle the centromere count by 1
  data$cent_CN=data$cent_CN+1
  data$cent_Area=data$cent_Area+SArea
  data$Her2_CN=data$Her2_CN+1
  data$Her2_Area=data$Her2_Area+SArea
  data$Her2_ArSpotFrac=data$Her2_Area/data$Her2_CN
  data$cent_ArSpotFrac=data$cent_Area/data$cent_CN
  data$Her2_centRatio=data$Her2_Area/data$cent_Area
  data
}

### Function to read the Manual Count Samples
ManCountRead=function(data){
  Input=read.csv(data, header=T)
  CNames=colnames(Input)
  SampNames=substr(CNames, 2, 5)
  tempNames=unique(SampNames)
  SampleStorage=vector("list", length(tempNames))
  #browser()
  names(SampleStorage)=tempNames
  SampleSumm=matrix(NA, length(tempNames), 10)
  colnames(SampleSumm)=c("Ecent17", "Rcent17", "EHer2", "RHer2", "Her2CentR", 
                         "pcScat", "pcMix", "pcClust", "pcER", "pcHER2")
  rownames(SampleSumm)=tempNames
  for (i in tempNames){
    x1=grep(i, SampNames)
    SampleStorage[[i]]=na.omit(Input[ ,x1])
    tnames=colnames(SampleStorage[[i]])
    t1=regmatches(tnames, regexpr("_[[:alnum:]]+", tnames))
    colnames(SampleStorage[[i]])=regmatches(t1, regexpr("[[:alnum:]]+", t1))
    SampleStorage[[i]]$ER=gsubfn(".", list("0"="neg", "1"="neg", "2"="neg",
                          "3"="pos", "4"="pos", "5"="neg", "6"="pos", "7"="pos"), 
                                 as.character(SampleStorage[[i]]$ER))
    SampleStorage[[i]]$ER=factor(SampleStorage[[i]]$ER)
    SampleStorage[[i]]$HER2prot=gsubfn(".", list("1"="pos", "2"="neg", "0"="weak"), 
                                 as.character(SampleStorage[[i]]$HER2prot)) 
    SampleStorage[[i]]$HER2prot=factor(SampleStorage[[i]]$HER2prot)
  
    SampleStorage[[i]]$arr=gsubfn(".", list("1"="clust", "11"="clust", "2"="scat",
                                            "3"="mix"), as.character(SampleStorage[[i]]$arr))
    SampleStorage[[i]]$arr=factor(SampleStorage[[i]]$arr)
   # browser()
    SampleSumm[i, ]=c(mean(SampleStorage[[i]]$cent17), diff(range(SampleStorage[[i]]$cent17)),
                      mean(SampleStorage[[i]]$Her2), diff(range(SampleStorage[[i]]$Her2)),
                      mean(SampleStorage[[i]]$Her2/SampleStorage[[i]]$cent17),
                      length(which(SampleStorage[[i]]$arr=="scat"))/nrow(SampleStorage[[i]]),
                      length(which(SampleStorage[[i]]$arr=="mix"))/nrow(SampleStorage[[i]]),
                      length(which(SampleStorage[[i]]$arr=="clust"))/nrow(SampleStorage[[i]]),
                      length(which(SampleStorage[[i]]$ER=="pos"))/nrow(SampleStorage[[i]]),
                      length(which(SampleStorage[[i]]$HER2prot=="pos"))/nrow(SampleStorage[[i]]))

    t1=which(SampleStorage[[i]]$Her2>20)
    SampleStorage[[i]]$Her2=as.character(SampleStorage[[i]]$Her2)
    SampleStorage[[i]]$Her2[t1]="amp"
    SampleStorage[[i]]$Her2=factor(SampleStorage[[i]]$Her2)
  }

  Output=list(SampleStorage=SampleStorage, Summary=data.frame(SampleSumm))
  Output
}


## Function to read in the distance files onlY:
TestSampsRead=function(dirPath, truePath, GFPath){
  FileIn=dir(dirPath)
  File2=dir(truePath)
  File3=dir(GFPath)
  Storage=list()
  tNames=substr(FileIn, 1, 4)

  for (i in 1:length(FileIn)){
  #browser()
      tempB=read.csv(sprintf('%s%s', dirPath, FileIn[i])) #,",")
  tempA=read.table(sprintf("%s%s", truePath, File2[i]), sep="\t", header=T)
  tempD=read.table(sprintf("%s%s", GFPath, File3[i]), sep="\t", header=T)
      tempD$memb_PixVar=(tempD$memb_OverLapIntensity-tempD$memb_nonLapIntensity)/
        tempD$memb_OverLapIntensity
      tempD$memb_PixVar[tempD$memb_PixVar=="-Inf"]=0
      tempD$memb_pos=(tempD$memb_MembraneIntensity-tempD$memb_MeanIntensity)/
        tempD$memb_MembraneIntensity
 #   browser()
      tempA$HER2.membrane.score[which(tempA$HER2.membrane.score==0)]="abs"
      tempA$HER2.membrane.score[which(tempA$HER2.membrane.score==1)]="com"
      tempA$HER2.membrane.score[which(tempA$HER2.membrane.score==2)]="bro"
      
    # tempD$HER2.membrane.score=factor(tempD$HER2.membrane.score)
   #   tempD$HER2.membrane.score<- factor(tempD$HER2.membrane.score)
   #   tempD$HER2.membrane.score=factor
  tempD2=merge(tempD, tempA, by.x="rownames", by.y="Manual_GF", all.x=T)
      colnames(tempD2)[which(colnames(tempD2)=="memb_temp")]="memb_BackAdjustedInt"
     # browser()
 # tempD$Mem2_PixelCV=tempD$Her2_range/tempD$Mem2_MembraneIntensity
   ## add the first two columns to the dataframe
 # L1=sapply(tempB, function(x) as.numeric(x[1:2]))  
#  rownames(L1)=c("Her2_ArVar", "Her2_ArVar_norm")
  #L1=data.frame(t(L1))
#   L1$trueARR=NA
#   idxV=as.numeric(as.character(na.omit(tempA$Manual_GF)))
#   idxT=which(is.na(tempA$Manual_GF)==F)
#     NaRm=which(is.na(idxV)==T)
#     if (length(NaRm)>0){
#     idxV=idxV[-NaRm]
#     idxT=idxT[-NaRm]
#     }

#   L1$trueARR[idxV]=tempA$HER2ARR[idxT]
#   L1$HER2CN=NA
#   L1$HER2CN[idxV]=tempA$HER2CN[idxT]
#   L1$true.memb=NA
#   L1$true.memb[idxV]=tempA$HER2.membrane.score[idxT]
#   L1$memb.Pos=NA
#   L1$memb.Pos[idxV]=tempA$HER2prot.score[idxT]
#     
#   tempC=sapply(tempB, function(x) length(x)-3)
#   L3=sapply(1:length(tempB), function(x) as.numeric(tempB[[x]][3:round(2+tempC[x]/2)]))
#   L2=sapply(1:length(tempB), function(x) as.numeric(tempB[[x]][round(3+tempC[x]/2):(2+tempC[x])]))
#   L4=rep(1:length(L3), times=unlist(lapply(L3, length)))  
#   df=data.frame(samp=tNames[i], cell=L4, DistEst=unlist(L3), CellSize=unlist(L2))
  
#   colnames(tempD)=gsub("Spot5", "Her2", colnames(tempD))    
#   colnames(tempD)=gsub("Mem2", "memb", colnames(tempD))  
#   colnames(tempD)[grep("memb_range", colnames(tempD))]="memb_PixelSD"  
#   L1=cbind(L1, tempD[ ,c(grep("Her2", colnames(tempD)), grep("memb", colnames(tempD)))], 
#            cent_Area=tempD[ ,"Spot3_Area"])
#   
#   Storage$ArVar[[i]]=L1 
   Storage[[i]]=tempD2
  }
#   names(Storage$ArVar)=tNames
   names(Storage)=tNames
  Storage
  
  ## create a list for specific cell spot areas and area fracs
  
}


QCFiles=function(dirP){
  dirN=dir(dirP)
  tempB=read.csv(sprintf('%s/%s', dirP, dirN[1]), header=T, stringsAsFactors=F)
  FN=unlist(regmatches(dirN[1], regexec("_[a-zA-Z]+", dirN[1])))
  colnames(tempB)[2:6]=paste(colnames(tempB)[2:6], FN, sep="")
  
  for (i in 2:length(dirN)){
    tempA=read.csv(sprintf('%s/%s', dirP, dirN[i]), header=T, stringsAsFactors=F)
    FN=unlist(regmatches(dirN[i], regexec("_[a-zA-Z]+", dirN[i])))
    colnames(tempA)[2:6]=paste(colnames(tempA)[2:6], FN, sep="")
   # browser()
    tempB=merge(tempB, tempA, by.x="Sample", by.y="Sample")
  }
  tempB$Sample=unlist(regmatches(tempB$Sample, regexec("[0-9]+_[a-z]+[1-5]*_[0-5]", tempB$Sample)))
  tempB$Pat=substr(tempB$Sample, 1, 4)
  tempB$type=substr(tempB$Sample, 6,6)
  tempB$type[tempB$type=="l"]="m"
  tempB
}


ReadManualData=function(data){
  Input=read.table(data, header=T, sep="\t")
  MeasNames=c("ER_IHC", "HER2_IHC", "HER2_IHC_score", "ER_IFISH", "ER_pc",
              "HER2_IFISH","HER2_pc", "Clust", "Scat", "Mix", "CellNo",
              "CNcent17", "CNHer2")
  colnames(Input)=c("ID", rep(MeasNames, 3))
  Input2=rbind(cbind(Input[ ,2:(length(MeasNames)+1)], type="Pre"),
               cbind(Input[ ,(length(MeasNames)+2):(2*length(MeasNames)+1)], type="Post"),
               cbind(Input[ ,(2*length(MeasNames)+2):(3*length(MeasNames)+1)], type="Met"))
  Input[1:5 ,"ID"]=paste("00", Input[1:5 ,"ID"], sep="" )
  Input2$ID=Input$ID
    x3=rowSums(Input2[ ,-sapply(c("ID", "type"), function(x) grep(x, colnames(Input2)))], na.rm=T)
  Input2=Input2[-which(x3==0), ]
  nidx=unlist(sapply(c("_pc", "Clust", "Scat", "Mix", "HER2_IHC_score"), function(x) grep(x, colnames(Input2))))
  x1=which(Input2[, -nidx]==1, arr.ind=T)
  Input2[ ,-nidx][x1]="pos"
  x1=which(Input2[, -nidx]==2, arr.ind=T)
  Input2[ ,-nidx][x1]="neg"
    #Input2=na.omit(Input2)
  Input2$Clust=as.numeric(Input2$Clust)
  Input2$Scat=as.numeric(Input2$Scat)
  Input2$Mix=as.numeric(Input2$Mix)
  Input2$Her2_dist="Het"  
  Input2$Her2_dist[Input2$Clust>50]="Clust"  
  Input2$Her2_dist[Input2$Scat>50]="Scat"
  Input2$Her2_dist[Input2$Mix>50]="Mix"
  Input2
}

ContTable=function(tab, title, chisqtest=F, xlabL="xlabel",ylabL="ylabel"){
  library("RColorBrewer")
  if (chisqtest==T){
    a1=chisq.test(tab)
    tit2=sprintf("Chisq = %g", round(a1$p.value*100)/100)
  } else {
    tit2=" "
  }
  nr=nrow(tab)
  nc=ncol(tab)
  if (chisqtest==T){
    l1=(a1$observed-a1$expected)
    l1[which(l1<0, arr.ind=T)]=0
    image(l1, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab=xlabL, ylab=ylabL,
          main=sprintf("%s %s", title, tit2))
  }else{
    image(scale(tab), col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab=xlabL, ylab=ylabL,
          main=sprintf("%s %s", title, tit2))}
  if (nr==1){
    xval=0
  }else{
  xval=seq(0, 1, 1/(nr-1))
  }
  if (nc==1){
    yval=0
  }else{
    yval=seq(0, 1, 1/(nc-1))
  }
  #browser()
  axis(1, at=xval, rownames(tab))
  axis(2, at=yval, colnames(tab))
  for(i in 1:nr){
    for (j in 1:nc){
      text(xval[i], yval[j], tab[i,j], cex=1.5)
    }
  }
}

EstCNVar=function(dataFile){
  x1=read.csv(dataFile, header=T, sep=",")
  cep17Col=seq(1, ncol(x1), 5)
  Her2Col=seq(2, ncol(x1), 5)
  
  TabSum=matrix(0, nrow=length(cep17Col), ncol=4)
  
  for (i in 1:length(cep17Col)){
    tc=cut(x1[ ,cep17Col[i]], c(0, 3, 1000), labels=c("norm", "amp"))
    td=cut(x1[ ,Her2Col[i]], c(0, 3, 1000), labels=c("norm", "amp"))
    Out=paste(td, tc, sep="_")
    t2=summary(factor(Out))/length(na.omit(tc))
    if ("amp_amp" %in% names(t2)){
    TabSum[i, 1]=t2["amp_amp"];}
    if ("amp_norm" %in% names(t2)){
    TabSum[i, 2]=t2["amp_norm"];}
    if ("norm_amp" %in% names(t2)){
    TabSum[i, 3]=t2["norm_amp"];}
    if ("norm_norm" %in% names(t2)){
    TabSum[i, 4]=t2["norm_norm"];}
  }
  rownames(TabSum)=substr(colnames(x1)[cep17Col], 2, 5)
  colnames(TabSum)=c("aa", "an", "na", "nn")
  TabSum
}