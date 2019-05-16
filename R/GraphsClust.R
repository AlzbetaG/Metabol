#' Cluster analysis
#'
#' Display cluster plots and heat maps.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If TRUE (default) quality control samples (QCs) are automatically distinguished. See Details.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @details If quality control samples (QCs) are present in data and QCs=TRUE, versions with QCs and without them are displayed. If QCs=TRUE and QCs are not present in data, this step is automatically skipped.
#' @return Cluster plots and heat maps.
#' @import dendextend
#' @importFrom stats hclust lowess dist as.dendrogram order.dendrogram
#' @importFrom gplots heatmap.2
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsClust(data,name,groupnames)
#' @export
GraphsClust=function(data,name,groupnames,tsf="clr",QCs=TRUE){

    dirout = paste(getwd(),"/",sep = "")

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
        dataS=as.matrix(data)
    }
#################################################################################################
    count=length(groupnames)
    basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown",
                "gray40","darkslateblue","firebrick","darkcyan","darkmagenta", "deeppink1","limegreen","gold2","bisque2",
                "lightcyan3","red","darkolivegreen3")            # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

    groups=NULL
    groupCodes=NULL
    color=NULL
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataS))
        gr=rep(i,length(Gr))
        groups=c(groups,gr)
        gc=rep(LETTERS[i],length(Gr))
        groupCodes=c(groupCodes,gc)
        cl=rep(basecolor[i],length(Gr))
        color=c(color,cl)
    }

    colorCodes=unique(color)
    names(colorCodes)=unique(groupCodes)
    #################################################################################################
    if (count < 6){
      wid=10
    }

    if (count > 5 & count < 9){
      wid=15
    }

    if (count > 8){
      wid=25
    }

##########################################################################################################################
if (QCs==TRUE){
    QC=grep("QC",rownames(dataS))
    if (length(QC)!=0){
        dataall=dataS
        groupsall=groups
        groupCodesall=groupCodes
        colorCodesall=colorCodes

        dataS=dataall[-QC,]
        groups=groupsall[-QC]
        groupCodes=groupCodesall[-QC]
        colorCodes=colorCodesall[unique(groups)]

        PDFQC=paste(dirout,"Cluster_QC_",name,".pdf",sep="")

        pdf((PDFQC),width=wid+5,height=10)
        par(mar = c(6, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

        clusters=hclust(dist(dataall), method="ward.D")
        dend <- as.dendrogram(clusters)
        labels_colors(dend) <- colorCodesall[groupCodesall][order.dendrogram(dend)]
        plot(dend, main = "Cluster Analysis",cex=0.5)
        dev.off()
    }
}
    #################################################################################################
    hei=round((length(colnames(data))/8),digits = 0)

    ##################################################################################
    PDF1=paste(dirout,"Cluster_",name,".pdf",sep="")
    PDF2=paste(dirout,"Heat_hier_",name,".pdf",sep="")
    PDF3=paste(dirout,"Heat_Ward_",name,".pdf",sep="")

    force(groupCodes)
    force(colorCodes)
    rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"), space = "rgb")

    pdf((PDF2),width=wid,height=hei)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    heatmap.2(t(dataS),col=rgb.palette(256), scale="row",margins = c(5.5, 15), key=TRUE, keysize = 1, symkey=FALSE, density.info="none",trace="none", cexRow=0.7,cexCol=0.8)

    dev.off()

    pdf((PDF3),width=wid,height=hei)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    heatmap.2(t(dataS),hclustfun=function(x) hclust(dist(x), method="ward.D"),col=rgb.palette(256), scale="row",margins = c(5.5, 15), key=TRUE, keysize = 1, symkey=FALSE, density.info="none",trace="none", cexRow=0.7,cexCol=0.8)

    dev.off()

    pdf((PDF1),width=wid+5,height=10)
    par(mar = c(6, 4, 4, 4) + 0.1,oma=c(1,1,1,1))

    clusters=hclust(dist(dataS), method="ward.D")
    dend <- as.dendrogram(clusters)
    labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
    plot(dend, main = "Cluster Analysis",cex=0.5)
    dev.off()

}

