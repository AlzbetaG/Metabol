#' ROC curves for combination of variables
#'
#' Display ROC curves with confidence intervals for chosen combination of variables. Variables are combined by logistic regression.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param selm Select sequence number of metabolites, their combination by logistic regression will be plotted in ROC curve. Up to six metabolites can be chosen.
#' @param selg Select sequence number of groups of data for which ROC curve will be plotted. The dafult is c(1,2).
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @return ROC curve of chosen combination of variables with confidence intervals of area under the curve (AUC).
#' @importFrom  pROC ci.se plot.roc
#' @importFrom stats binomial glm
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsROCcomb(data,name,groupnames,selm=c(2,3))
#' @export
GraphsROCcomb=function(data,name,groupnames,selm,selg=c(1,2),tsf="clr",QCs=FALSE){

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

    i=selg[1]
    j=selg[2]

    pocet=length(selm)

    groups2=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
    datav=c(Group[i,1]:Group[i,2],Group[j,1]:Group[j,2])
    datavyb=dataS[datav,]

    nazvy=colnames(datavyb)

    if (pocet==2) {
        regr = glm(groups2 ~ datavyb[,selm[1]]*datavyb[,selm[2]], family=binomial(logit))
        met=paste(nazvy[selm[1]],"&",nazvy[selm[2]],sep=" ")
        met2=paste(nazvy[selm[1]],"_",nazvy[selm[2]],sep="")
    }

    if (pocet==3) {
        regr = glm(groups2 ~ datavyb[,selm[1]]*datavyb[,selm[2]]*datavyb[,selm[3]], family=binomial(logit))
        met=paste(nazvy[selm[1]],"&",nazvy[selm[2]],"&",nazvy[selm[3]],sep=" ")
        met2=paste(nazvy[selm[1]],"_",nazvy[selm[2]],"_",nazvy[selm[3]],sep="")
    }

    if (pocet==4) {
        regr = glm(groups2 ~ datavyb[,selm[1]]*datavyb[,selm[2]]*datavyb[,selm[3]]*datavyb[,selm[4]], family=binomial(logit))
        met=paste(nazvy[selm[1]],"&",nazvy[selm[2]],"&",nazvy[selm[3]],"&",nazvy[selm[4]],sep=" ")
        met2=paste(nazvy[selm[1]],"_",nazvy[selm[2]],"_",nazvy[selm[3]],"_",nazvy[selm[4]],sep="")
    }

    if (pocet==5) {
        regr = glm(groups2 ~ datavyb[,selm[1]]*datavyb[,selm[2]]*datavyb[,selm[3]]*datavyb[,selm[4]]*datavyb[,selm[5]],
                   family=binomial(logit))
        met=paste(nazvy[selm[1]],"&",nazvy[selm[2]],"&",nazvy[selm[3]],"&",nazvy[selm[4]],"&",nazvy[selm[5]],sep=" ")
        met2=paste(nazvy[selm[1]],"_",nazvy[selm[2]],"_",nazvy[selm[3]],"_",nazvy[selm[4]],"_",nazvy[selm[5]],sep="")
    }

    if (pocet==6) {
        regr = glm(groups2 ~ datavyb[,selm[1]]*datavyb[,selm[2]]*datavyb[,selm[3]]*datavyb[,selm[4]]*datavyb[,selm[5]]*datavyb[,selm[6]],
                   family=binomial(logit))
        met=paste(nazvy[selm[1]],"&",nazvy[selm[2]],"&",nazvy[selm[3]],"&",nazvy[selm[4]],"&",nazvy[selm[5]],"&",nazvy[selm[6]],sep=" ")
        met2=paste(nazvy[selm[1]],"_",nazvy[selm[2]],"_",nazvy[selm[3]],"_",nazvy[selm[4]],"_",nazvy[selm[5]],"_",nazvy[selm[6]],sep="")
    }

    PDF=paste("ROC_",name,"_",groupnames[i],"_",groupnames[j],"_",met2,".pdf",sep="")

    pdf((PDF),width=10,height=10)
    rocobj <- plot.roc(groups2, regr$fitted.values,main=met,percent=TRUE,ci=TRUE,print.auc=FALSE)
    sens.ci <- ci.se(rocobj, specificities=seq(0, 100, 5),boot.n=100,conf.level=0.95, stratified=FALSE)
    plot(sens.ci, type="shape", col="lavender")
    text(12,0,paste("AUC: ",round(rocobj$auc,2),"% (",round(rocobj$ci[1],1),"% - ",round(rocobj$ci[3],1),"%)",sep=""),font=2)
    dev.off()

    setwd(dirout)
}


