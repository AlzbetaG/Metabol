#' Permutation test of PLS-DA and OPLS-DA
#'
#' Makes permutation test of PLS-DA and OPLS-DA, displays score plots (for OPLS-DA only), permutation plots and summary table.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param type Definition which type of method do you want to use. If "opls" is set (default), OPLS-DA is performed and the number of orthogonal components is automatically computed by using the cross-validation (with a maximum of 9 orthogonal components). If two groups are not nicely separated, the algorithm returns error. In this case set this parameter to "err".  When "pls" is set, PLS-DA is done.
#' @param permut Number of random permutations of response labels to estimate R2Y and Q2Y significance by permutation testing. Default is 100.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return Score plots and permutation plots of OPLS-DA.
#' @return Excel file with summary.
#' @import ropls
#' @importFrom openxlsx write.xlsx
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @references Thevenot E. A. et al., Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses,Journal of Proteome Research 14, 3322-3335, 2015.
#' @export

GraphsPermut=function(data,name,groupnames,tsf="clr",type="opls",permut=100,QCs=FALSE){

  dirout=getwd()

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
  dirout2 = paste(dirout,"/","Permutation_test",sep = "")
  dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
  setwd(dirout2)

  #################################################################################################

  if (type=="opls"){
      pre=1
      ort=NA
  }

  if (type=="err"){
      pre=1
      ort=1
  }

  if (type=="pls"){
      pre=NA
      ort=0
  }

  randOPLSDAall=opls(dataM,groups,predI=pre,permI = permut,orthoI = ort,scaleC="pareto",printL=FALSE,plotL=FALSE)
  plot(randOPLSDAall, typeVc = "permutation",file.pdfC = paste("Permutation_",type,"_all_",name,".pdf",sep=""))
  if (type=="opls" | type=="err"){plot(randOPLSDAall, typeVc = "x-score",file.pdfC = paste("Scores_",type,"_all_",name,".pdf",sep=""))}
  summyall=randOPLSDAall@summaryDF
  rownames(summyall)="All"

  write.xlsx(summyall,file = paste("Permutation_",type,"_all_",name,".xlsx",sep=""),
             sheetName="Permutation",col.names=TRUE, row.names=TRUE, append=FALSE,
             showNA=TRUE)


  table=NULL
  if (count>2){
  for (i in 1:(count-1)){
    for (j in (i+1):count){
      OP=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
      dataOP=dataM[c((Group[i,1]:Group[i,2]),(Group[j,1]:Group[j,2])),]
      colors2=color[c(Group[i,1],Group[j,1])]
      marks2=marks[c((Group[i,1]:Group[i,2]),(Group[j,1]:Group[j,2]))]
      groupnames2=groupnames[c(i,j)]
      noteOP=paste(name,"_",groupnames[i],"_",groupnames[j],sep="")
      randOPLSDA=opls(dataOP,OP,predI=pre,permI = permut,orthoI = ort,plotL=FALSE,scaleC="pareto",printL=FALSE)
      plot(randOPLSDA, typeVc = "permutation",file.pdfC = paste("Permutation_",type,"_",noteOP,".pdf",sep=""))
      if (type=="opls" | type=="err"){plot(randOPLSDA, typeVc = "x-score",file.pdfC = paste("Scores_",type,"_",noteOP,".pdf",sep=""))}
      summy=randOPLSDA@summaryDF
      rownames(summy)=paste(groupnames[i],"_",groupnames[j],sep="")
      table=rbind(table,summy)
      remove(OP)
      remove(dataOP)
      remove(noteOP)
     }
  }

   write.xlsx(table,file = paste("Permutation_",type,"_",name,".xlsx",sep=""),
             sheetName="Permutation",col.names=TRUE, row.names=TRUE, append=FALSE,
             showNA=TRUE)

  }
  setwd(dirout)
}
