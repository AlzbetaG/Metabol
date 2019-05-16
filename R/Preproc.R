#' Data preprocessing
#'
#' Counting zeros, interpolation by LOESS, checking of coefficient of variation (CV) of quality control samples (QCs), zero imputatuion, summary of data preprocessing.
#' @param data Data table sorted according to batch order (in columns) with variables (metabolites) in rows.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param CVlim Limit of coefficient of variation (CV) of quality control samples (QCs). All metabolites with CV of QCs higher than this limit will be omitted from the analysis. The default value is 30.
#' @details Up to fifteen different groups can be distinguished in data (including QCs).
#' @return LOESS curves.
#' @return Excel file with preprocessed data.
#' @import openxlsx
#' @importFrom stats sd loess loess.smooth
#' @examples data=batchSort
#' name="Preprocessing"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' Preproc(data,name,groupnames)
#' @export
Preproc=function(data,name,groupnames,CVlim=30){
#################################################################################################################################################
# Excel with outputs
#################################################################################################################################################
dirout2=paste(getwd(),"/",sep = "")

wb <- createWorkbook()
sheet1  <- addWorksheet(wb, sheetName = "00 Summary")
sheet2  <- addWorksheet(wb, sheetName = "01 Zeros")
sheet3  <- addWorksheet(wb, sheetName = "02 Big zeros off")
sheet4  <- addWorksheet(wb, sheetName = "06 Interpolated")
sheet5  <- addWorksheet(wb, sheetName = "07 Interpolated cleaned")
sheet6  <- addWorksheet(wb, sheetName = "08 CV")
sheet7  <- addWorksheet(wb, sheetName = "09 CV below 30")
sheet8  <- addWorksheet(wb, sheetName = "10 Final Zeros")
sheet9  <- addWorksheet(wb, sheetName = "11 Final no Zeros")

file00=paste(dirout2, "Preproc_",name,".xlsx",sep="")

################################################################################################################################
# Groups of data
################################################################################################################################
# Done automatically (color, Group etc.).
# Needed structure: string of letters (group) and numbers (sample), e.g. con1, con2, ..., pat1, ...
# Names must be unique, overlaps are permitted,
# Maximal number of gorups is 6!!!

sum1=nrow(data)

blank=c(grep("Blank",colnames(data)),grep("blank",colnames(data)))
blankname=c(grep("Blank",groupnames),grep("blank",groupnames))

if (length(blank)!=0){
  data=data[,-blank]
  if(length(blankname)!=0){
    groupnamesall=groupnames[-blankname]
  }else{
    groupnamesall=groupnames
  }
}else{
  groupnamesall=groupnames
}

countall=length(groupnamesall)                                      # Number of groups

#################################################################################################################################################
# Deleting of zero metabolites
#################################################################################################################################################

zeros=as.matrix(apply(data, 1, function(c)(sum(c==0)/ncol(data))*100))
colnames(zeros)="% of zeros"
dataz=cbind(zeros,data)

dll=NULL
for (i in 1:countall){
    Gr=grep(groupnamesall[[i]],colnames(dataz))
    dl=as.matrix(c(rep(NA,nrow(dataz))))
    colnames(dl)=paste(groupnamesall[[i]],"-",length(Gr))
    for(j in 1:nrow(dataz)){
        if (length(which(dataz[j,Gr]==0))==0) {dl[j]=NA
        } else {
            dl[j]=length(which(dataz[j,Gr]==0))}
    }
    dll=cbind(dll,dl)
}

dataz=cbind(dll,dataz)

writeData(wb,sheet2,dataz,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

datazf=NULL
for (i in 1:nrow(dataz)){
    if (dataz[i,(countall+1)]<30){
        datazf=rbind(datazf,dataz[i,-c(1:(countall+1))])
    }
}

writeData(wb,sheet3,datazf,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

data=datazf
RowNames=rownames(data)

sum2=nrow(datazf)

##########################################################################################################################################

QCi=grep("QC",colnames(data))
QCi
dataQC=as.matrix(data[,QCi])

groupnames=groupnamesall[-grep("QC",groupnamesall)]
count=length(groupnames)

basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown","gray40","darkslateblue","firebrick","darkcyan","darkmagenta",
            "deeppink1","limegreen")           # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

Group=matrix(rep(NA,count*3),ncol=3)
colnames(Group)=c("min","max","length")
rownames(Group)=groupnames
ggroup=NULL
datag=NULL
color=NULL
for (i in 1:count){
    Gr=grep(groupnames[i],colnames(data))
    ggroup=c(ggroup,Gr)
    Group[i,3]=length(Gr)
    if (i==1){
        Group[i,1]=1
        Group[i,2]=Group[i,3]
    }
    else {
        Group[i,1]=(1+Group[(i-1),2])
        Group[i,2]=(Group[(i-1),2]+Group[i,3])
    }
    datag1=as.matrix(data[,Gr])
    datag=cbind(datag,datag1)
    cl=rep(basecolor[i],length(Gr))
    color=c(color,cl)
}

##########################################################################################################################################
# Basic graphs
##########################################################################################################################################

pdf(paste("03_Raw_",name,".pdf",sep=""))           # Original data
for(i in 1:nrow(data)){
  limx=c(0,ncol(data))
  ma=max(data[i,])
  mi=min(data[i,])
  konst=abs(ma-mi)
  m=max(data[i,])+konst/10
  limy=c(0,m)
  plot(t(data[i,]),pch=16,main=RowNames[i],xlab="Index",ylab="Area",xlim=limx,ylim=limy)
  par(new=TRUE)
  plot(dataQC[i,]~QCi,pch=17,xlab="",xlim=limx,ylim=limy,col="red",main="",ylab="",cex=1.5)
  par(new=TRUE)
  plot(datag[i,]~ggroup,pch=16,xlab="",xlim=limx,ylim=limy,col=color,main="",ylab="")
  legend("topleft",legend = c(groupnames,"QC"), pch = c(rep(16,count),17), col = c(unique(color),"red"),cex=0.75)
}
dev.off()

##################   LOESS   ###############################################

pdf(paste("04_Smooth_",name,".pdf",sep=""))        # LOESS curve
for(i in 1:nrow(data)){
  limx=c(0,ncol(data))
  ma=max(data[i,])
  mi=min(data[i,])
  konst=abs(ma-mi)
  m=max(data[i,])+konst/10
  limy=c(0,m)
  plot(t(data[i,]),pch=16,main=RowNames[i],xlab="Index",ylab="Area",xlim=limx,ylim=limy)
  par(new=TRUE)
  plot(datag[i,]~ggroup,pch=16,xlab="",xlim=limx,ylim=limy,col=color,main="",ylab="")
  par(new=TRUE)
  plot(dataQC[i,]~QCi,pch=17,xlab="",xlim=limx,ylim=limy,col="red",main="",ylab="",cex=1.5)
  fitd=loess(dataQC[i,]~QCi,span=0.75,degree=2)
  pred=predict(fitd,se=TRUE)
  lines(pred$fit~QCi,xlim=limx,ylim=limy,col=6)
  smooth=loess.smooth(QCi,dataQC[i,],span=0.75,degree=2, evaluation = ncol(data),family = "gaussian")
  par(new=TRUE)
  plot(smooth,pch=1,xlim=limx,ylim=limy,xlab="", ylab="",col=4)
  #smooth2=loess.smooth(QCi,dataQC[i,],span=0.75,degree=2, evaluation = ncol(data),family = "symmetric")
  #par(new=TRUE)
  #plot(smooth2,pch=1,xlim=limx,ylim=limy,xlab="", ylab="",col=3)
  legend("topleft",legend = c(groupnames,"QC","LOESS curve"), pch = c(rep(16,count),17,1), col = c(unique(color),"red",4),cex=0.75)
  }
dev.off()

##############################################################################################################

data=as.matrix(data)            # Important!!!

#################################################################################################################################################
# Smoothing factors
#################################################################################################################################################

r=matrix(rep(0,ncol(data)),ncol(data),nrow(data))
row.names(r)=colnames(data)
colnames(r)=rownames(data)
for(i in 1:nrow(data)){
    for (j in 1:ncol(data)){
        dat=data[i,]
        smooth=loess.smooth(QCi,dataQC[i,],span=0.75,degree=2, evaluation = ncol(data),family = "gaussian")
        smo=smooth$y
        r[j,i]=dat[j]/smo[j]
    }}

rr=matrix(rep(0,nrow(data)),nrow(data),1)
colnames(rr)="sm.fact"
for(i in 1:nrow(data)){
    dat=data[i,]
    smooth=loess.smooth(QCi,dataQC[i,],span=0.75,degree=2, evaluation = ncol(data),family = "gaussian")
    smomax=max(smooth$y)
    smomin=min(smooth$y)
    rr[i,1]=smomax/smomin
}

newtable=cbind(rr,t(r))

writeData(wb,sheet4,newtable,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

##########################################################################################################################################
# Negative smoothing factors
#################################################################################################################################################

dataSet=NULL
for (i in 1:nrow(newtable)){
    if (newtable[i,1]>0){
        dataSet=rbind(dataSet,newtable[i,-1])
        rownames(dataSet)[i]=rownames(newtable)[i]
    }
    else{
        ddataset=rep(0,(ncol(newtable)-1))
        dataSet=rbind(dataSet,ddataset)
    }
}

idx = which(rowSums(dataSet) == 0)

if (length(idx) != 0) {
    dataSet=dataSet[-idx,]
} else {
    dataSet = dataSet}

sum3=nrow(dataSet)

writeData(wb,sheet5,dataSet,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

##########################################################################################################################

dataSetQC=as.matrix(dataSet[,QCi])
RowNames=rownames(dataSet)

pdf(paste("05_Interpol_",name,".pdf",sep=""))
for(i in 1:nrow(dataSet)){
    dat=dataSet[i,]
    smooth=loess.smooth(QCi,dataSetQC[i,],span=0.75,degree=2, evaluation = ncol(dataSet),family = "gaussian")
    smo=smooth$y
    r=rep(0,length(smo))
    for (j in 1:ncol(dataSet)) {r[j]=dat[j]/smo[j]}
    mm=max(r)
    limyy=c(0,mm)
    rQC=r[QCi]
    rother=r[-QCi]
    other=(1:ncol(dataSet))[-QCi]
    plot(rQC~QCi,xlim=limx,ylim=limyy,xlab="Index",ylab="Original data/smooth parameter",pch=17, col=2,main=RowNames[i],cex=1.5)
    par(new=TRUE)
    plot(rother~other,xlim=limx,ylim=limyy,xlab="",ylab="",pch=16,main="")
    legend("topleft",legend = c("QC","Samples"), pch = c(17,16), col = c("red","black"),cex=0.75)
}
dev.off()

#################################################################################################################################################
# QCs with CV > 30%
#################################################################################################################################################

QC=grep("QC", colnames(dataSet))
QC2=dataSet[,QC]

co.var<-function(x){                                # Funkcion for CV
  100*apply(x,1,sd)/apply(x,1,mean)
}

CV=co.var(QC2)
QC007=cbind(CV,dataSet)

writeData(wb,sheet6,QC007,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

QC30=matrix(rep(0,(nrow(QC007))*ncol(QC007)),ncol=ncol(QC007))
dim(QC30)

colnames(QC30)=colnames(QC007)
rownamesQC30=rownames(dataSet)
row.names(QC30)=rownamesQC30

for(j in 1:nrow(QC007)){
  for(k in 1:ncol(QC007)){
    if ( CVlim>QC007[j,1]) { QC30[j,k]=QC007[j,k]}
  }}

##############################################

idx = which(rowSums(QC30) == 0)

if (length(idx) != 0) {
    QC30b=QC30[-idx,]
} else {
    QC30b = QC30}

colnames(QC30b)=colnames(QC007)
QC30b=as.data.frame(QC30b)
CV=round(QC30b$CV,2)
QC30b=QC30b[,-1]
CV=matrix(CV,nrow=nrow(QC30b))
rownames(CV)=rownames(QC30b)
dataSet=QC30b

#################################################################################################################################################
# New ordering according to groups
#################################################################################################################################################

dataSetg=NULL
for (i in 1:count){
    poradi=ggroup[Group[i,1]:Group[i,2]]
    dataSetg=cbind(dataSetg,as.matrix(dataSet[,poradi]))
}

dataQC=as.matrix(dataSet[,QCi])

dataSet2=cbind(dataSetg,dataQC)

sum4=nrow(dataSet2)

writeData(wb,sheet7,dataSet2,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

#################################################################################################################################################
# Transponation
#################################################################################################################################################
dataSetM=as.matrix(dataSet2)
dataSetM=t(dataSetM)

writeData(wb,sheet8,dataSetM,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

#################################################################################################################################################
# Zero imputation
#################################################################################################################################################

for (k in 1:count){
    vyber=Group[k,1]:Group[k,2]
    dl=c(rep(0,ncol(dataSetM)))
    for(j in 1:ncol(dataSetM)){
        sl=dataSetM[vyber,j]
        dl[j]=min(sl[sl>0])
    }

    for(i in vyber){
        for(j in 1:ncol(dataSetM)){
            if(dataSetM[i,j]==0){
                dataSetM[i,j]=(2/3)*dl[j]
            }
        }
    }
}

# Zeros in QCs
qqc=length(QCi)
groupQC=nrow(dataSetM)-qqc+1:qqc

dl=c(rep(0,ncol(dataSetM)))
for(j in 1:ncol(dataSetM)){
    sl=dataSetM[groupQC,j]
    dl[j]=min(sl[sl>0])
}

for(i in groupQC){
    for(j in 1:ncol(dataSetM)){
        if(dataSetM[i,j]==0){
            dataSetM[i,j]=(2/3)*dl[j]
        }
    }
}

writeData(wb,sheet9,dataSetM,colNames = TRUE, rowNames = TRUE)
#saveWorkbook(wb, file = file00,overwrite = TRUE)

write.table(dataSetM, file = paste("11_final_",name,".csv",sep=""), sep=";",col.names =NA)     # Final table

pos <- 1
envir = as.environment(pos)
assign("dataSetM",get("dataSetM"), envir =envir)

#################################################################################################################################################
# Summary
#################################################################################################################################################

sumtab1=rbind(sum1,sum2,sum3,sum4)
sumtab2=rbind(NA,sum1-sum2,sum2-sum3,sum3-sum4)
sumtab3=rbind(NA,round(100*(sum1-sum2)/sum1,2),round(100*(sum2-sum3)/sum2,2),round(100*(sum3-sum4)/sum3,2))
sumtab4=rbind(NA,sum1-sum2,sum1-sum3,sum1-sum4)
sumtab5=rbind(NA,round(100*(sum1-sum2)/sum1,2),round(100*(sum1-sum3)/sum1,2),round(100*(sum1-sum4)/sum1,2))
sumtab=cbind(sumtab1,sumtab2,sumtab3,sumtab4,sumtab5)
rownames(sumtab)=c("At the beginning","Zero metabolites off","Negative smoothing factors off","CV > 30% off")
colnames(sumtab)=c("Number","Off_absolute","Off_percet","Total_absolute","Total_percent")

writeData(wb,sheet1,sumtab,colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = file00,overwrite = TRUE)

groupnames=c(groupnames,"QC")

pos <- 1
envir = as.environment(pos)
assign("groupnames",get("groupnames"), envir =envir)

save(dataSetM,groupnames,file=paste("Preproc_",name,".RData",sep=""))

#save.image(file=paste("Preproc_",name,".RData",sep=""),.GlobalEnv)

}
