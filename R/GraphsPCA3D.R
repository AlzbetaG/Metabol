#' 3D principal component analysis (PCA)
#'
#' Makes 3D principal component analysis (PCA).
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param type A type of plots must be defined by "points" (default) or "names".
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If TRUE (default) quality control samples (QCs) are are left in the graph. If FALSE QCs are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @details If quality control samples (QCs) are present in data and QCs=TRUE, versions with QCs and without them are displayed. If QCs=TRUE and QCs are not present in data, this step is automatically skipped.
#' @return 3D score plot of PCA.
#' @import rgl
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @export
GraphsPCA3D=function(data,name,groupnames,type="points",tsf="clr",QCs=TRUE){

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
##########################################################################################################################
    count=length(groupnames)
    basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown",
                "gray40","darkslateblue","firebrick","darkcyan","darkmagenta", "deeppink1","limegreen","gold2","bisque2",
                "lightcyan3","red","darkolivegreen3")           # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

    groups=NULL
    color=NULL
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataS))
        gr=rep(i,length(Gr))
        groups=c(groups,gr)
        cl=rep(basecolor[i],length(Gr))
        color=c(color,cl)
    }


##########################################################################################################################
    if (QCs==FALSE){
        QC=grep("QC",rownames(dataS))
        if (length(QC)!=0){
            dataS=dataS[-QC,]
            color=color[-QC]
            groups=groups[-QC]
            groupnames=groupnames[unique(groups)]
            count=length(groupnames)
        }
    }

##########################################################################################################################

    L=svd(dataS)$u[,c(1:3)]
    K=diag(svd(dataS)$d[1:3])
    M=svd(dataS)$v[,c(1:3)]
    G=(sqrt(nrow(dataS)-1))*L
    H=(1/sqrt(nrow(dataS)-1))*M%*%K

    rownames(G)=rownames(dataS)
    rownames(H)=colnames(dataS)

    Ge=svd(var(dataS))$u
    Z=dataS%*%Ge
    PC1=var(Z[,1])/sum(apply(Z,2,var))
    PC2=var(Z[,2])/sum(apply(Z,2,var))
    PC3=var(Z[,3])/sum(apply(Z,2,var))

    if (type == "points") {
        PDF=paste(dirout,"Graphs-PCA3D_points - ",name,".pdf",sep="")

    # window settings - optional
    rgl.open()
    par3d(windowRect = 50 + c(0, 0, 640, 640))        # size of the window
    rgl.viewpoint(theta = 45, phi = 30, fov = 50, zoom = 0.95) # view point...

    # optional changes in axis and brackground window colour
    bg3d("lightgrey", top=T)   # colour of backgroung
    rgl.pop("lights")   # styles of light and colour
    light3d(specular = "blue")

    # plotting 3D PCA
    plot3d(G,size=1.5 ,type="s",main = paste("PCA - ",name, sep="")
       ,xlab=paste("PC1 - ", round(PC1*100, digits = 2), "%; Cumulative = ",round(PC1*100, digits = 2)+round(PC2*100, digits = 2)+round(PC3*100, digits = 2),"%")
       ,ylab=paste("PC2 - ", round(PC2*100, digits = 2),"%"),zlab=paste("PC3 - ", round(PC3*100, digits = 2),"%")
       ,col=color,show.plane=T,box=T,axes=T,top=T)
    legend3d("topleft",legend = groupnames, pch = 19, col = unique(color))

    #rgl.postscript(PDF,fmt="pdf",drawText=TRUE) # vytvori pdf
    #rgl.postscript("PCA CML.eps",fmt="eps",drawText=TRUE) # vytvori eps
    }

    if (type == "names") {
        PDF=paste(dirout,"Graphs-PCA3D_names - ",name,".pdf",sep="")

        names=rownames(dataS)

        # window settings - optional
        rgl.open()
        par3d(windowRect = 50 + c(0, 0, 640, 640))        # size of the window
        rgl.viewpoint(theta = 45, phi = 30, fov = 50, zoom = 0.95) # view point...

        # optional changes in axis and brackground window colour
        bg3d("lightgrey", top=T)   # colour of backgroung
        rgl.pop("lights")   # styles of light and colour
        light3d(specular = "blue")

        # plotting 3D PCA
        plot3d(G,size=1.5 ,type="n",main = paste("PCA - ",name, sep="")
               ,xlab=paste("PC1 - ", round(PC1*100, digits = 2), "%; Cumulative = ",round(PC1*100, digits = 2)+round(PC2*100, digits = 2)+round(PC3*100, digits = 2),"%")
               ,ylab=paste("PC2 - ", round(PC2*100, digits = 2),"%"),zlab=paste("PC3 - ", round(PC3*100, digits = 2),"%")
               ,col=color,show.plane=T,box=T,axes=T,top=T)
        text3d(G,texts=names,col=color)
        legend3d("topleft",legend = groupnames, pch = 19, col = unique(color))

        #rgl.postscript(PDF,fmt="pdf",drawText=TRUE) # vytvori pdf
        #rgl.postscript("PCA CML.eps",fmt="eps",drawText=TRUE) # vytvori eps
    }
}

