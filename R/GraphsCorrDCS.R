#' Correlation graphs with class separation
#'
#' Evaluate correlations between all pairs of variables and degree of class separation.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return Correlation graph with evaluated clusters.
#' @return Excel file with degree of class separation (DCS) of all possible pairs of groups in data for all pairs of comparisons of two variables.
#' @importFrom robCompositions cenLR
#' @importFrom stats cor cor.test
#' @importFrom openxlsx write.xlsx
#' @import car
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @references Pierce, K.M. et al. (2005) Classification of gasoline data obtained by gas chromatography using a piecewise alignment algorithm combined with feature selection and principal component analysis, J CHROMATOGR A 1096, p. 101-110.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsCorrDSC(data,name,groupnames)
#' @export
GraphsCorrDSC=function(data,name,groupnames,tsf="clr",QCs=FALSE){

    ####################################################################################################
    dirout = paste(getwd(),"/",sep = "")
    PDF=paste(dirout,"Corr_DCS_",name,".pdf",sep="")

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
    ##########################################################################################################################
    basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown",
                "gray40","darkslateblue","firebrick","darkcyan","darkmagenta", "deeppink1","limegreen","gold2","bisque2",
                "lightcyan3","red","darkolivegreen3")           # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    basemarks=c(15,17,18,8,11,2,0,16,5,6,4,10,3,7,9,12)

    count=length(groupnames)
    groups=NULL
    marks=NULL
    color=NULL
    groupss=matrix(rep(NA,count*2),ncol=2)
    colnames(groupss)=c("min","max")
    rownames(groupss)=groupnames
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataM))
        groupss[i,1]=min(Gr)
        groupss[i,2]=max(Gr)
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
            groupss=matrix(rep(NA,count*2),ncol=2)
            colnames(groupss)=c("min","max")
            rownames(groupss)=groupnames
            for (i in 1:count){
                Gr=grep(groupnames[i],rownames(dataM))
                groupss[i,1]=min(Gr)
                groupss[i,2]=max(Gr)
            }
        }
    }


    ####################################################################################################
    names=rownames(dataM)

    pdf(PDF)
    DCStable=NULL
    DCSmet=NULL
    for (i in 1:(ncol(dataM)-1)){
        for(j in (i+1):ncol(dataM)){
            mmean=matrix(c(rep(0,count*2)),nrow=count)
            dev=matrix(c(rep(0,count)),nrow=count)
            for(k in 1:count){
                x=dataM[(groupss[k,1]:groupss[k,2]),i]
                y=dataM[(groupss[k,1]:groupss[k,2]),j]
                scores=cbind(x,y)
                mmean[k,1]=mean(scores[,1])
                mmean[k,2]=mean(scores[,2])
                dev[k,]=var(dist(rbind(mmean[k,],scores),method = "euclidean")[1:nrow(scores)])
            }

            DCS=NULL
            for(l in 1:(count-1)){
                for(m in (l+1):count){
                    d=dist(rbind(mmean[l,],mmean[m,]),method = "euclidean")[1]
                    DCS2=d/(sqrt(dev[l,]+dev[m,]))
                    DCS=cbind(DCS,DCS2)
                }
            }

            nazvy=NULL
            for(l in 1:(count-1)){
                for(m in (l+1):count){
                    col=paste(groupnames[l],"_",groupnames[m],sep="")
                    nazvy=c(nazvy,col)
                }
            }

            DCSmet=rbind(DCSmet,cbind(colnames(dataM)[i],colnames(dataM)[j]))
            DCStable=rbind(DCStable,DCS)

            gr=as.factor(groups)
            dataEllipse(dataM[,i],dataM[,j],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE, center.cex=0.2, col=unique(color),
                        pch=unique(marks),center.pch=3,main=paste(colnames(dataM)[i],"&",colnames(dataM)[j],sep=" "),xlab=colnames(dataM)[i],ylab=colnames(dataM)[j],grid=FALSE)
            legend("topleft",legend = paste("DCS_",nazvy[1:count], "=",round(DCS,2)),cex=0.75)
            text(dataM[,i],dataM[,j],names,pos=1,cex=0.5)
        }
    }
    dev.off()

    pocet=ncol(dataM)*(ncol(dataM)-1)/2

    DCSord=matrix(1:pocet,ncol=1)
    colnames(DCStable)=c(nazvy)
    colnames(DCSmet)=c("Met_1","Met_2")
    colnames(DCSord)="Order"

    DCStablefin=list(DCSord,DCSmet,DCStable)

    write.xlsx(DCStablefin,file = paste(dirout,"Corr_DCS_",name,".xlsx", sep=""),sheetName="DCS",
               col.names=TRUE, row.names=FALSE, append=FALSE, showNA=TRUE)
}

