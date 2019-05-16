#' Outlier analysis of PCA and PLS-DA
#'
#' Outlier analysis of PCA and PLS-DA.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param type Definition which type of method do you want to use: "pls" (default) and "pca". If two groups are not nicely separated, the algorithm returns error. In this case set this parameter to "err".
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return Plot of outliers.
#' @import ropls
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @references Thevenot E. A. et al., Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses,Journal of Proteome Research 14, 3322-3335, 2015.
#' @export

GraphsOutlier=function(data,name,groupnames,tsf="clr",type="pls",QCs=FALSE){

    dirout=getwd()

    dirout2 = paste(dirout,"/","Outlier",sep = "")
    dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
    setwd(dirout2)
  ##########################################################################################################################
  if (tsf=="clr"){
    dataM=cenLR(data)$x.clr
  }

  if (tsf=="log"){
    dataM=log(data)
  }

  if (tsf=="log10"){
    dataM=log10(data)
  }

  if (tsf=="pareto"){
    dataM=scale(data, scale=TRUE, center=FALSE)
  }

    PQN1=function (x){
        xref=apply(x,2,median)
        podil=x
        for (i in 1:ncol(x)){
            podil[,i]=x[,i]/xref[i]
        }
        s=apply(podil,1,median)
        PQN=x
        for (j in 1:nrow(x)){
            PQN[j,]=x[j,]/s[j]
        }
        return(PQN)
    }

    if (tsf=="PQN"){
        dataM=PQN1(data)
    }

    if (tsf=="lnPQN"){
        dataM=log(PQN1(data))
    }

    if (tsf=="none"){
        dataM=data
    }
  #################################################################################################
  count=length(groupnames)
  basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown",
              "gray40","darkslateblue","firebrick","darkcyan","darkmagenta", "deeppink1","limegreen","gold2","bisque2",
              "lightcyan3","red","darkolivegreen3")           # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
  basemarks=c(15,17,18,8,11,2,0,16,5,6,4,10,3,7,9,12)

  Group=matrix(rep(NA,count*3),ncol=3)
  colnames(Group)=c("min","max","length")
  rownames(Group)=groupnames
  groups=NULL
  marks=NULL
  color=NULL
  for (i in 1:count){
    Gr=grep(groupnames[i],rownames(dataM))
    Group[i,1]=min(Gr)
    Group[i,2]=max(Gr)
    Group[i,3]=length(Gr)
    gr=rep(i,length(Gr))
    groups=c(groups,gr)
    zn=rep(basemarks[i],length(Gr))
    marks=c(marks,zn)
    cl=rep(basecolor[i],length(Gr))
    color=c(color,cl)
  }

  ##########################################################################################################################
  if (QCs==FALSE){
    QC=grep("QC",rownames(dataM))
    if (length(QC)!=0){
      dataM=dataM[-QC,]
      color=color[-QC]
      marks=marks[-QC]
      groups=groups[-QC]
      groupnames=groupnames[unique(groups)]
      count=length(groupnames)
      Group=matrix(rep(NA,count*3),ncol=3)
      colnames(Group)=c("min","max","length")
      rownames(Group)=groupnames
      for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataM))
        Group[i,1]=min(Gr)
        Group[i,2]=max(Gr)
        Group[i,3]=length(Gr)
      }
    }
  }

  #################################################################################################
  if (type=="pca"){
      randOPLSDA=opls(dataM,predI=NA,orthoI = 0,plotL=FALSE,scaleC="pareto",printL=FALSE)
      plot(randOPLSDA, typeVc = "outlier",file.pdfC = paste("Outlier_PCA_all_",name,".pdf",sep=""))

      if (count>2){
      for (i in 1:(count-1)){
          for (j in (i+1):count){
              OP=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
              dataOP=dataM[c((Group[i,1]:Group[i,2]),(Group[j,1]:Group[j,2])),]
              noteOP=paste(name,"_",groupnames[i],"_",groupnames[j],sep="")
              randOPLSDA=opls(dataOP,predI=NA,orthoI = 0,plotL=FALSE,scaleC="pareto",printL=FALSE)
              plot(randOPLSDA, typeVc = "outlier",file.pdfC = paste("Outlier_PCA_",noteOP,".pdf",sep=""))
              remove(OP)
              remove(dataOP)
              remove(noteOP)
          }
      }
  }
}

  if (type=="pls"){
      randOPLSDA=opls(dataM,groups,predI=NA,orthoI = 0,plotL=FALSE,scaleC="pareto",printL=FALSE)
      plot(randOPLSDA, typeVc = "outlier",file.pdfC = paste("Outlier_PLS_all_",name,".pdf",sep=""))

      if (count>2){
      for (i in 1:(count-1)){
          for (j in (i+1):count){
              OP=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
              dataOP=dataM[c((Group[i,1]:Group[i,2]),(Group[j,1]:Group[j,2])),]
              noteOP=paste(name,"_",groupnames[i],"_",groupnames[j],sep="")
              randOPLSDA=opls(dataOP,OP,predI=NA,orthoI = 0,plotL=FALSE,scaleC="pareto",printL=FALSE)
              plot(randOPLSDA, typeVc = "outlier",file.pdfC = paste("Outlier_PLSDA_",noteOP,".pdf",sep=""))
              remove(OP)
              remove(dataOP)
              remove(noteOP)
          }
      }
      }
  }

  if (type=="err"){
      randOPLSDA=opls(dataM,groups,predI=2,orthoI = 0,plotL=FALSE,scaleC="pareto",printL=FALSE)
      plot(randOPLSDA, typeVc = "outlier",file.pdfC = paste("Outlier_PLS_all_",name,".pdf",sep=""))

      if (count>2){
          for (i in 1:(count-1)){
              for (j in (i+1):count){
                  OP=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
                  dataOP=dataM[c((Group[i,1]:Group[i,2]),(Group[j,1]:Group[j,2])),]
                  noteOP=paste(name,"_",groupnames[i],"_",groupnames[j],sep="")
                  randOPLSDA=opls(dataOP,OP,predI=2,orthoI = 0,plotL=FALSE,scaleC="pareto",printL=FALSE)
                  plot(randOPLSDA, typeVc = "outlier",file.pdfC = paste("Outlier_PLSDA_",noteOP,".pdf",sep=""))
                  remove(OP)
                  remove(dataOP)
                  remove(noteOP)
              }
          }
      }
  }


  setwd(dirout)

}
