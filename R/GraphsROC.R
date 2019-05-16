#' ROC curves
#'
#' Display ROC curves with confidence intervals for every variable.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details ROC curves can be used only for comparison of two groups. If there is more groups in data, all possible combinations of pairs are evaluated.
#' @return ROC curves with confidence intervals of area under the curve (AUC).
#' @return Excel sheet with summary of AUCs.
#' @importFrom  pROC ci.se plot.roc
#' @importFrom robCompositions cenLR
#' @import openxlsx
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsROC(data,name,groupnames)
#' @export
GraphsROC=function(data,name,groupnames,tsf="clr",QCs=FALSE){

    ############################################################################################
    dirout=getwd()
    dirout2 = paste(dirout,"/","ROC",sep = "")
    dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
    setwd(dirout2)

    ##########################################################################################################################
    if (tsf=="clr"){
        dataM=cenLR(data)$x.clr
        dataS=scale(dataM, scale=FALSE, center=TRUE)
    }

    if (tsf=="log"){
        dataM=log(data)
        dataS=scale(dataM, scale=FALSE, center=TRUE)
    }

    if (tsf=="log10"){
        dataM=log10(data)
        dataS=scale(dataM, scale=FALSE, center=TRUE)
    }

    if (tsf=="pareto"){
        dataS=scale(data, scale=TRUE, center=TRUE)
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
        dataS=scale(dataM, scale=FALSE, center=TRUE)
    }

    if (tsf=="lnPQN"){
        dataM=log(PQN1(data))
        dataS=scale(dataM, scale=FALSE, center=TRUE)
    }

    if (tsf=="none"){
        dataM=data
        dataS=scale(data, scale=TRUE, center=TRUE)
    }
    ############################################################################################

    count=length(groupnames)
    Group=matrix(rep(NA,count*3),ncol=3)
    colnames(Group)=c("min","max","length")
    rownames(Group)=groupnames
    groups=NULL
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataS))
        Group[i,1]=min(Gr)
        Group[i,2]=max(Gr)
        Group[i,3]=length(Gr)
        gr=rep(i,length(Gr))
        groups=c(groups,gr)
    }

##########################################################################################################################
if (QCs==FALSE){
    QC=grep("QC",rownames(dataS))
    if (length(QC)!=0){
        dataS=dataS[-QC,]
        groups=groups[-QC]
        groupnames=groupnames[unique(groups)]
        count=length(groupnames)
        Group=matrix(rep(NA,count*3),ncol=3)
        colnames(Group)=c("min","max","length")
        rownames(Group)=groupnames
        for (i in 1:count){
            Gr=grep(groupnames[i],rownames(dataS))
            Group[i,1]=min(Gr)
            Group[i,2]=max(Gr)
            Group[i,3]=length(Gr)
        }
    }
}

    ############################################################################################

    for (i in 1:(count-1)){
        for (j in (i+1):count){
    groups2=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
    vyber=c(Group[i,1]:Group[i,2],Group[j,1]:Group[j,2])
    datavyb=dataS[vyber,]

    PDF=paste("ROC_",name,"_",groupnames[i],"_",groupnames[j],".pdf",sep="")

    wb <- createWorkbook()
    sheet1  <- addWorksheet(wb, sheetName = paste("Summary_",groupnames[i],"_",groupnames[j],sep=""))
    sheet2  <- addWorksheet(wb, sheetName = paste("AUC_",groupnames[i],"_",groupnames[j],sep=""))

    file00=paste("ROC_",name,"_",groupnames[i],"_",groupnames[j],".xlsx", sep="")

    nazvy=colnames(datavyb)

    tableAUC=matrix(rep(0,ncol(datavyb)),ncol=1)
    colnames(tableAUC)="AUC"
    met=NULL

    pdf((PDF),width=10,height=10)
    for (k in 1:ncol(datavyb)){
        rocobj <- plot.roc(groups2, datavyb[,k],main=nazvy[k],percent=TRUE,ci=TRUE,print.auc=FALSE)
        sens.ci <- ci.se(rocobj, specificities=seq(0, 100, 1),boot.n=100,conf.level=0.95, stratified=FALSE)
        plot(sens.ci, type="shape", col="lavender")
        text(12,0,paste("AUC: ",round(rocobj$auc,2),"% (",round(rocobj$ci[1],1),"% - ",round(rocobj$ci[3],1),"%)",sep=""),font=2)

        tableAUC[k,1]=rocobj$auc
        met=rbind(met,nazvy[k])
    }
    dev.off()

    ord=matrix(1:ncol(datavyb),ncol=1)
    colnames(ord)="Order"
    colnames(met)="Metabolite"

    tableAUCfin=list(ord,met,tableAUC)

    writeData(wb,sheet2,tableAUCfin,colNames = TRUE, rowNames = TRUE)

    ################################################################################################################################
    #summary:

    ssum1= tableAUC >= 70 & tableAUC < 80
    sum1=length(tableAUC[ssum1==TRUE])

    ssum2=tableAUC >= 80 & tableAUC < 90
    sum2=length(tableAUC[ssum2==TRUE])

    ssum3=tableAUC >= 90
    sum3=length(tableAUC[ssum3==TRUE])

    sumA=rbind(sum3,sum2,sum1,sum1+sum2+sum3)
    sumB=rbind(round((100*sum3/ncol(datavyb)),2),round(100*sum2/ncol(datavyb),2),round(100*sum1/ncol(datavyb),2),
               round(100*(sum1+sum2+sum3)/ncol(datavyb),2))
    sum=cbind(sumA,sumB)
    rownames(sum)=c("AUC in [90,100]","AUC in [80,90) ","AUC in [70,80)","Sum")
    colnames(sum)=c("Number","Percent")

    writeData(wb,sheet1,sum,colNames = TRUE, rowNames = TRUE)
    saveWorkbook(wb, file = file00,overwrite = TRUE)

        }
    }

    setwd(dirout)
}


