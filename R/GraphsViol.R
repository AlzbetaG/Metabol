#' Violin plots
#'
#' Display violin plots for all specific groups of data.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return Violin plots.
#' @importFrom robCompositions cenLR
#' @import sm
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsViol(data,name,groupnames)
#' @export
GraphsViol=function(data,name,groupnames,tsf="clr",QCs=FALSE){

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
    ##################################################################################
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

############################################################################################
    if (count < 6){
        wid=10
    }

    if (count > 5 & count < 9){
        wid=15
    }

    if (count > 8){
        wid=20
    }

##################################################################################
names<-colnames(dataM)

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

##########################################################################################################################
    PDF1=paste("Viol_",name,".pdf",sep="")
    PDF2=paste("Viol_P_",name,".pdf",sep="")
############################################################################################
 "viol"<- function (x, range = 1.5, h = NULL, ylim = NULL, names = NULL,
              horizontal = FALSE, col = "magenta", border = "black", lty = 1,
              lwd = 1, rectCol = "black", colMed = "white", pchMed = 19,
              at, add = FALSE, wex = 1, drawRect = TRUE)
    {
     if(length(col)==1) col <- rep(col,n)
        datas <- x
        n <- length(datas)
        if (missing(at))
            at <- 1:n
        upper <- vector(mode = "numeric", length = n)
        lower <- vector(mode = "numeric", length = n)
        q1 <- vector(mode = "numeric", length = n)
        q3 <- vector(mode = "numeric", length = n)
        med <- vector(mode = "numeric", length = n)
        base <- vector(mode = "list", length = n)
        height <- vector(mode = "list", length = n)
        baserange <- c(Inf, -Inf)
        args <- list(display = "none")
        if (!(is.null(h)))
            args <- c(args, h = h)
        for (i in 1:n) {
            data <- datas[[i]]
            data.min <- min(data)
            data.max <- max(data)
            q1[i] <- quantile(data, 0.25)
            q3[i] <- quantile(data, 0.75)
            med[i] <- median(data)
            iqd <- q3[i] - q1[i]
            upper[i] <- min(q3[i] + range * iqd, data.max)
            lower[i] <- max(q1[i] - range * iqd, data.min)
            est.xlim <- c(min(lower[i], data.min), max(upper[i],
                                                       data.max))
            smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
                                             args))
            hscale <- 0.4/max(smout$estimate) * wex
            base[[i]] <- smout$eval.points
            height[[i]] <- smout$estimate * hscale
            t <- range(base[[i]])
            baserange[1] <- min(baserange[1], t[1])
            baserange[2] <- max(baserange[2], t[2])
        }
        if (!add) {
            xlim <- if (n == 1)
                at + c(-0.5, 0.5)
            else range(at) + min(diff(at))/2 * c(-1, 1)
            if (is.null(ylim)) {
                ylim <- baserange
            }
        }
        if (is.null(names)) {
            label <- 1:n
        }
        else {
            label <- names
        }
        boxwidth <- 0.05 * wex
        if (!add)
            plot.new()
        if (!horizontal) {
            if (!add) {
                plot.window(xlim = xlim, ylim = ylim)
                axis(2)
                axis(1, at = at, labels = label)
            }
            box()
            for (i in 1:n) {
                polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
                        c(base[[i]], rev(base[[i]])), col = col[i], border = border,
                        lty = lty, lwd = lwd)
                if (drawRect) {
                    lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
                          lty = lty)
                    rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
                         q3[i], col = rectCol)
                    points(at[i], med[i], pch = pchMed, col = colMed)
                }
            }
        }
        else {
            if (!add) {
                plot.window(xlim = ylim, ylim = xlim)
                axis(1)
                axis(2, at = at, labels = label)
            }
            box()
            for (i in 1:n) {
                polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
                                                        rev(at[i] + height[[i]])), col = col, border = border,
                        lty = lty, lwd = lwd)
                if (drawRect) {
                    lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
                          lty = lty)
                    rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
                             boxwidth/2, col = rectCol)
                    points(med[i], at[i], pch = pchMed, col = colMed)
                }
            }
        }
        invisible(list(upper = upper, lower = lower, median = med,
                       q1 = q1, q3 = q3))
    }
############################################################################################

pdf((PDF1),width=10,height=10)
par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

for(i in 1:ncol(dataM)){
    dataviol=vector("list", count)
    names(dataviol) <- groupnames
    for (j in 1:count){
        dataviol[[j]]=dataM[Group[j,1]:Group[j,2],i]
    }
    viol(dataviol,names=groupnames,col=unique(color))
    title(names[i])
}
dev.off()

pdf((PDF2),width=10,height=10)
par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

for(i in 1:ncol(dataM)){
    dataviol=vector("list", count)
    names(dataviol) <- groupnames
    for (j in 1:count){
        dataviol[[j]]=dataM[Group[j,1]:Group[j,2],i]
    }
    viol(dataviol,names=groupnames,col=c("transparent","transparent"))
    stripchart(dataviol, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
    title(names[i])
}
dev.off()
}


