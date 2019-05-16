#' Correlation graphs
#'
#' Evaluate correlations (Pearson and Spearman) between all pairs of variables, display correlation plots.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @details Cor_S denotes Spearman correlation coefficient.
#' @return Correlation graph.
#' @return Excel file with Pearson and Spearman correlations in standard tables and also sorted in one column. Significance tests are also evaluated.
#' @import openxlsx
#' @importFrom robCompositions cenLR
#' @importFrom stats cor cor.test
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsCorr(data,name,groupnames)
#' @export
GraphsCorr=function(data,name,groupnames,tsf="clr",QCs=FALSE){

    ####################################################################################################
    dirout = paste(getwd(),"/",sep = "")
    ####################################################################################################
    PDF1 = paste(dirout,"Corr_",name,".pdf",sep="")

    wb <- createWorkbook()
    sheet1  <- addWorksheet(wb, sheetName = "Corr P")
    sheet2  <- addWorksheet(wb, sheetName = "Corr S")
    sheet3  <- addWorksheet(wb, sheetName = "Test P")
    sheet4  <- addWorksheet(wb, sheetName = "Test S")
    sheet5  <- addWorksheet(wb, sheetName = "List corr P")
    sheet6  <- addWorksheet(wb, sheetName = "List corr S")

    file00=paste(dirout, "Corr_",name,".xlsx",sep="")

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
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataM))
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
        }
    }

    ####################################################################################################

    corelace=cor(dataM)

    writeData(wb,sheet1,corelace,colNames = TRUE, rowNames = TRUE)

    corelaces=cor(dataM,method="spearman")

    writeData(wb,sheet2,corelaces,colNames = TRUE, rowNames = TRUE)

    nazvy=colnames(dataM)
    names=rownames(dataM)

    corelace2=NULL
    corelaces2=NULL
    cortestp=NULL
    cortests=NULL
    met=NULL

    pdf(PDF1)
    for (i in 1:(ncol(dataM)-1)){
        for(j in (i+1):ncol(dataM)){
            x=dataM[,i]
            y=dataM[,j]

            co=cor(x,y)
            corelace2=rbind(corelace2,co)
            coo=round(co,3)
            cos=cor(x,y,method="spearman")
            corelaces2=rbind(corelaces2,cos)
            coos=round(cos,3)

            testp=cor.test(x,y,alternative = "two.sided", method = "pearson", conf.level = 0.95)$p.value
            cortestp=rbind(cortestp,testp)
            tests=cor.test(x,y,alternative = "two.sided", method = "spearman", conf.level = 0.95)$p.value
            cortests=rbind(cortests,tests)

            plot(x,y,pch=marks,xlab=nazvy[i],ylab=nazvy[j],main=paste(nazvy[i],"&",ylab=nazvy[j],sep=" "),col=color)
            text(x,y,names,pos=1,cex=0.5)
            legend("topleft",c(paste("cor=",coo),paste("cor_s=",coos)),cex=1.25)

            met=rbind(met,cbind(nazvy[i],nazvy[j]))
        }
    }
    dev.off()

    pocet=ncol(dataM)*(ncol(dataM)-1)/2
    ord=matrix(1:pocet,ncol=1)
    colnames(ord)="Order"
    colnames(met)=c("Met_1","Met_2")
    rownames(met)=paste("A",1:pocet,sep="")

    colnames(cortestp)="P_value"
    cortestpfin=list(ord,met,cortestp)

    writeData(wb,sheet3,cortestpfin,colNames = TRUE, rowNames = TRUE)

    colnames(cortests)="P_value"
    cortestsfin=list(ord,met,cortests)

    writeData(wb,sheet4,cortestsfin,colNames = TRUE, rowNames = TRUE)

    colnames(corelace2)="Cor"
    corelace2fin=list(ord,met,corelace2)

    writeData(wb,sheet5,corelace2fin,colNames = TRUE, rowNames = TRUE)

    colnames(corelaces2)="Cor"
    corelaces2fin=list(ord,met,corelaces2)

    writeData(wb,sheet6,corelaces2fin,colNames = TRUE, rowNames = TRUE)
    saveWorkbook(wb, file = file00,overwrite = TRUE)

}

