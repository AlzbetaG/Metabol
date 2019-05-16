#' Partial least squares - discriminant analysis (PLS-DA)
#'
#' Makes partial least squares - discriminant analysis (PLS-DA), displays score plots, loading plots and biplots.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param type A type of plots must be defined by "points" (default), "names" or "both".
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param top How many rays with highest lengths should be in biplot? The default is 20.
#' @param QCs logical. If TRUE (default) quality control samples (QCs) are automatically distinguished. See Details.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @details If quality control samples (QCs) are present in data and QCs=TRUE, versions with QCs and without them are displayed. If QCs=TRUE and QCs are not present in data, this step is automatically skipped.
#' @return Score plot and biplot of PLS-DA. Cross validation plot.
#' @return Excel file with lengths of rays in biplot and mean square error of prediction (MSEP) of data.
#' @import car
#' @import pls
#' @import openxlsx
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsPLSDA(data,name,groupnames)
#' @export
GraphsPLSDA=function(data,name,groupnames,type="points",tsf="clr",top=20,QCs=TRUE){

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
        dataM=data
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
    if (QCs==TRUE){
    QC=grep("QC",rownames(dataM))
    if (length(QC)!=0){
        dataMall=dataM
        colorall=color
        marksall=marks
        groupsall=groups
        groupnamesall=groupnames
        countall=count

        dataM=dataMall[-QC,]
        if (tsf=="pareto"){
            dataSall=scale(dataMall, scale=TRUE, center=TRUE)
        }

        if (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN"|tsf=="lnPQN"){
            dataSall=scale(dataMall, scale=FALSE, center=TRUE)
        }

        if (tsf=="none"){
            dataSall=as.matrix(dataMall)
        }
        color=colorall[-QC]
        marks=marksall[-QC]
        groups=groupsall[-QC]
        groupnames=groupnamesall[unique(groups)]
        count=length(groupnames)

        Groupall=c(rep(NA,countall))
        for (i in 1:countall){
            Gr=grep(groupnamesall[i],rownames(dataSall))
            Groupall[i]=length(Gr)
        }

        if (countall==2) {
            ycon=rep(1,Groupall[1])
            ycon2=rep(0,Groupall[2])
            yd=c(ycon,ycon2)
        }else{
        yd=NULL
        for (k in 1:countall){
            yy=matrix(rep(0,Groupall[k]*countall),ncol=countall)
            yy[1:Groupall[k],k]=1
            yd=rbind(yd,yy)
        }
        }

        dirout = paste(getwd(),"/",sep = "")
        PDFQC=paste(dirout,"PLSDA_QC_",name,".pdf",sep="")

        "biplot.color.QC.P" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- dimnames(x)[[1]]
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                # text(x, xlabs, cex = cex[1], col = col, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "red", ...)
                if (var.axes)
                    ######################
                # Farbe der Variablen-Pfeile
                arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "red",
                       length = arrow.len)
                invisible()
            }


        "biplot.color_transp.QC.P" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- dimnames(x)[[1]]
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                # text(x, xlabs, cex = cex[1], col = col, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "transparent", ...)
                if (var.axes)
                    ######################
                # Farbe der Variablen-Pfeile
                arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "transparent",
                       length = arrow.len)
                invisible()
            }

        "biplot.color.QC.N" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=NULL, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- dimnames(x)[[1]]
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                #  plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                #  cex=pch.cex, ...)
                plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, xlab = xlabMY, ylab = ylabMY, ...)
                ######################
                # kein Text bei Objekten
                text(x, xlabs, cex = cex[1], col = col, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "red", ...)
                if (var.axes)
                    ######################
                # Farbe der Variablen-Pfeile
                arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "red",
                       length = arrow.len)
                invisible()
            }

        "biplot.color_transp.QC.N" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=NULL, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- dimnames(x)[[1]]
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                #  plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                #  cex=pch.cex, ...)
                plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, xlab = xlabMY, ylab = ylabMY, ...)
                ######################
                # kein Text bei Objekten
                text(x, xlabs, cex = cex[1], col = col, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "transparent", ...)
                if (var.axes)
                    ######################
                # Farbe der Variablen-Pfeile
                arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "transparent",
                       length = arrow.len)
                invisible()
            }

        "biplot.color.QC.B" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- rownames(x)
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                text(x, xlabs, cex = cex[1], col = col,pos=3, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "red", ...)
                if (var.axes)
                    ######################
                # Farbe der Variablen-Pfeile
                arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "red",
                       length = arrow.len)
                invisible()
            }


        "biplot.color_transp.QC.B" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- rownames(x)
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                text(x, xlabs, cex = cex[1], col = col,pos=3, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "transparent", ...)
                if (var.axes)
                    ######################
                # Farbe der Variablen-Pfeile
                arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "transparent",
                       length = arrow.len)
                invisible()
            }

        pdf((PDFQC),width=10,height=10)
        ys=as.matrix(scale(yd,scale=FALSE))
        resst <- mvr(ys~dataSall,ncomp=2,method="simpls")

        G=resst$scores
        H=resst$loadings

        rownames(G)=rownames(dataSall)

        #Most discriminating vectors:

        lengthx=matrix(rep(0,nrow(H)),nrow=nrow(H))

        for(i in 1:nrow(H)){
            for(j in 1:nrow(H)){
                a=abs(H[i,1])
                b=abs(H[i,2])
                lengthx[i]=sqrt(a^2+b^2)
                rownames(lengthx)=rownames(H)
            }}

        discrim=as.matrix(lengthx[order(lengthx[,1],decreasing = TRUE),])

        matrixH=cbind(lengthx,H)
        matrixH=matrixH[order(matrixH[,1],decreasing = TRUE),]
        matrixH=matrixH[,-1]

        number=length(which((Groupall==1)==TRUE))
        if (number!=0){
            if (type=="points"){
                biplot.color.QC.P(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color.QC.P(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color_transp.QC.P(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                         ,xlabMY="PC1"
                                         ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                dev.off()
            }
            if (type=="names"){
                biplot.color.QC.N(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color.QC.N(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color_transp.QC.N(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                         ,xlabMY="PC1"
                                         ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                dev.off()
            }
            if (type=="both"){
                biplot.color.QC.B(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color.QC.B(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color_transp.QC.B(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                         ,xlabMY="PC1"
                                         ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                dev.off()
            }
        }else{
            if (type=="points"){
                biplot.color.QC.P(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color.QC.P(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color_transp.QC.P(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                         ,xlabMY="PC1"
                                         ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                gr=as.factor(groupsall)
                dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PLSDA - ",name),
                            xlab="PC1",ylab="PC2")
                legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.95, plot.points=TRUE, cex=1.2,center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PLSDA - ",name),
                            xlab="PC1",ylab="PC2")
                legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)

                dev.off()
            }
            if (type=="names"){
                biplot.color.QC.N(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color.QC.N(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color_transp.QC.N(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                         ,xlabMY="PC1"
                                         ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                gr=as.factor(groupsall)
                dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PLSDA - ",name),
                            xlab="PC1",ylab="PC2")
                legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.95, plot.points=TRUE, cex=1.2,center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PLSDA - ",name),
                            xlab="PC1",ylab="PC2")
                legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)

                dev.off()
            }
            if (type=="both"){
                biplot.color.QC.B(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color.QC.B(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                  ,xlabMY="PC1"
                                  ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                biplot.color_transp.QC.B(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PLSDA - ",name, sep="")
                                         ,xlabMY="PC1"
                                         ,ylabMY="PC2", pch=marksall,lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                gr=as.factor(groupsall)
                dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PLSDA - ",name),
                            xlab="PC1",ylab="PC2")
                legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.95, plot.points=TRUE, cex=1.2,center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PLSDA - ",name),
                            xlab="PC1",ylab="PC2")
                legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)

                dev.off()
            }
        }

    }
    }
##########################################################################################################################
    count=length(groupnames)

    if (tsf=="pareto"){
        dataS=scale(dataM, scale=TRUE, center=TRUE)
    }

    if (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN"|tsf=="lnPQN"){
        dataS=scale(dataM, scale=FALSE, center=TRUE)
    }

    if (tsf=="none"){
        dataS=as.matrix(dataM)
    }

    Group=matrix(rep(NA,count*3),ncol=3)
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(dataS))
        Group[i,1]=min(Gr)
        Group[i,2]=max(Gr)
        Group[i,3]=length(Gr)
    }

    if (count==2) {
        ycon=rep(1,Group[1,3])
        ycon2=rep(0,Group[2,3])
        yd=c(ycon,ycon2)
    }else{
        yd=NULL
        for (k in 1:count){
            yy=matrix(rep(0,Group[k,3]*count),ncol=count)
            yy[1:Group[k,3],k]=1
            yd=rbind(yd,yy)
        }
    }

##########################################################################################################################

    if (type == "names") {
    dirout = paste(getwd(),"/",sep = "")
    #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)
    PDF=paste(dirout,"PLSDA_n_",name,".pdf",sep="")
    PDF2=paste(dirout,"PLSDA_CV_",name,".pdf",sep="")

    wb <- createWorkbook()
    sheet1  <- addWorksheet(wb, sheetName = "Discrim")
    sheet2  <- addWorksheet(wb, sheetName = "MSEP")

    file00=paste(dirout, "PLSDA_",name,".xlsx",sep="")

    "biplot.color.N" <-
        function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                  xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                  arrow.len = 0.1, pch=NULL, pch.cex=1.2, xlabMY, ylabMY, ...)
        {
            n <- nrow(x)
            p <- nrow(y)
            if (is.null(xlabs)) {
                xlabs <- dimnames(x)[[1]]
                if (is.null(xlabs))
                    xlabs <- 1:n
            }
            xlabs <- as.character(xlabs)
            dimnames(x) <- list(xlabs, dimnames(x)[[2]])
            if (missing(ylabs)) {
                ylabs <- dimnames(y)[[1]]
                if (is.null(ylabs))
                    ylabs <- paste("Var", 1:p)
            }
            ylabs <- as.character(ylabs)
            dimnames(y) <- list(ylabs, dimnames(y)[[2]])
            if (length(cex) == 1)
                cex <- c(cex, cex)
            if (missing(col)) {
                col <- par("col")
                if (!is.numeric(col))
                    col <- match(col, palette())
                col <- c(col, col + 1)
            }
            else if (length(col) == 1)
                col <- c(col, col)
            unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
            rangx1 <- unsigned.range(x[, 1])
            rangx2 <- unsigned.range(x[, 2])
            rangy1 <- unsigned.range(y[, 1])
            rangy2 <- unsigned.range(y[, 2])
            if (missing(xlim) && missing(ylim))
                xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
            else if (missing(xlim))
                xlim <- rangx1
            else ylim <- rangx2
            ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
            on.exit(par(oldpar))
            oldpar <- par(pty = "s")
            ######################
            # Plot-Symbol definiert
            #  plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
            #  cex=pch.cex, ...)
            plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, xlab = xlabMY, ylab = ylabMY, ...)
            ######################
            # kein Text bei Objekten
            text(x, xlabs, cex = cex[1], col = col, ...)
            par(new = TRUE)
            plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                     ratio, xlab = "", ylab = "", col = col, ...)
            axis(3, col = 1)
            axis(4, col = 1)
            box(col = 1)
            ######################
            # Text der Variablen
            text(y, labels = ylabs, cex = cex[2], col = "red", ...)
            if (var.axes)
                ######################
            # Farbe der Variablen-Pfeile
            arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "red",
                   length = arrow.len)
            invisible()
        }

    "biplot.color_transp.N" <-
        function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                  xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                  arrow.len = 0.1, pch=NULL, pch.cex=1.2, xlabMY, ylabMY, ...)
        {
            n <- nrow(x)
            p <- nrow(y)
            if (is.null(xlabs)) {
                xlabs <- dimnames(x)[[1]]
                if (is.null(xlabs))
                    xlabs <- 1:n
            }
            xlabs <- as.character(xlabs)
            dimnames(x) <- list(xlabs, dimnames(x)[[2]])
            if (missing(ylabs)) {
                ylabs <- dimnames(y)[[1]]
                if (is.null(ylabs))
                    ylabs <- paste("Var", 1:p)
            }
            ylabs <- as.character(ylabs)
            dimnames(y) <- list(ylabs, dimnames(y)[[2]])
            if (length(cex) == 1)
                cex <- c(cex, cex)
            if (missing(col)) {
                col <- par("col")
                if (!is.numeric(col))
                    col <- match(col, palette())
                col <- c(col, col + 1)
            }
            else if (length(col) == 1)
                col <- c(col, col)
            unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
            rangx1 <- unsigned.range(x[, 1])
            rangx2 <- unsigned.range(x[, 2])
            rangy1 <- unsigned.range(y[, 1])
            rangy2 <- unsigned.range(y[, 2])
            if (missing(xlim) && missing(ylim))
                xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
            else if (missing(xlim))
                xlim <- rangx1
            else ylim <- rangx2
            ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
            on.exit(par(oldpar))
            oldpar <- par(pty = "s")
            ######################
            # Plot-Symbol definiert
            #  plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
            #  cex=pch.cex, ...)
            plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, xlab = xlabMY, ylab = ylabMY, ...)
            ######################
            # kein Text bei Objekten
            text(x, xlabs, cex = cex[1], col = col, ...)
            par(new = TRUE)
            plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                     ratio, xlab = "", ylab = "", col = col, ...)
            axis(3, col = 1)
            axis(4, col = 1)
            box(col = 1)
            ######################
            # Text der Variablen
            text(y, labels = ylabs, cex = cex[2], col = "transparent", ...)
            if (var.axes)
                ######################
            # Farbe der Variablen-Pfeile
            arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "transparent",
                   length = arrow.len)
            invisible()
        }

    #################################################################################

    ys=as.matrix(scale(yd,scale=FALSE))

    resst <- mvr(ys~dataS,ncomp=2,method="simpls")

    G=resst$scores
    H=resst$loadings

    rownames(G)=rownames(dataS)

    # Most discriminating vectors:
    lengthx=matrix(rep(0,nrow(H)),nrow=nrow(H))

    for(i in 1:nrow(H)){
        for(j in 1:nrow(H)){
            a=abs(H[i,1])
            b=abs(H[i,2])
            lengthx[i]=sqrt(a^2+b^2)
            rownames(lengthx)=rownames(H)
        }}

    discrim=as.matrix(lengthx[order(lengthx[,1],decreasing = TRUE),])

    matrixH=cbind(lengthx,H)
    matrixH=matrixH[order(matrixH[,1],decreasing = TRUE),]
    matrixH=matrixH[,-1]

    colnames(discrim)="Length"

    writeData(wb,sheet1,discrim,colNames = TRUE, rowNames = TRUE)

    pdf((PDF),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    number=length(which((Group[,3]==1)==TRUE))
    if (number!=0){
        biplot.color.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                                ,xlabMY="PC1",ylabMY="PC2", pch=groups,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

        biplot.color.N(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=groups,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

        biplot.color_transp.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                            ,xlabMY="PC1",ylabMY="PC2", pch=groups,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

        plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
        abline(h=0,col="gray86")
        abline(v=0,col="gray86")
        text(H[,1],H[,2], labels = rownames(H),cex=0.8)

        dev.off()
    }else {
        if (ncol(dataS)>10){
        if (nrow(ys)<=20){
            ncompp=nrow(ys)-count
            seg=nrow(ys)
        } else {
            ncompp=15
            seg=10
        }

        misclass=mvr(ys~dataS,ncomp=ncompp,method="simpls",validation = "CV",segments=seg)

        MSEPy=MSEP(misclass,estimate="CV",intercept=FALSE)

        if (count == 2) {
            MSEPf=MSEPy$val[1,,1:ncompp]

            writeData(wb,sheet2,MSEPf,colNames = TRUE, rowNames = TRUE)
            }else {
            MSEPf=NULL
            for (i in 1:count){
                MSEPff=MSEPy$val[1, i, 1:ncompp]
                MSEPf=cbind(MSEPf,MSEPff)
            }

            colnames(MSEPf)=paste("Y",1:count,sep="")

            writeData(wb,sheet2,MSEPf,colNames = TRUE, rowNames = TRUE)
            }

        pdf((PDF2),width=10,height=10)
        plot(MSEPy,main="Cross validation - MSEP")
        #axis(1, at=1:15,labels = 1:15)
        dev.off()
        }

    biplot.color.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                 ,xlabMY="PC1",ylabMY="PC2", pch=groups,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

    biplot.color.N(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                 ,xlabMY="PC1",ylabMY="PC2", pch=groups,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

    biplot.color_transp.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                        ,xlabMY="PC1",ylabMY="PC2", pch=groups,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

    gr=as.factor(groups)
    dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=rep(19,length(groupnames)),center.pch=3,main=paste("PLSDA - ",name),
                xlab="PC1",ylab="PC2")
    legend("bottomright",legend = groupnames, pch = 19, col = unique(color), cex=0.6)

    dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.95, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=rep(19,length(groupnames)),center.pch=3,main=paste("PLSDA - ",name),
                xlab="PC1",ylab="PC2")
    legend("bottomright",legend = groupnames, pch = 19, col = unique(color), cex=0.6)

    plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
    abline(h=0,col="gray86")
    abline(v=0,col="gray86")
    text(H[,1],H[,2], labels = rownames(H),cex=0.8)

    dev.off()
    }
   saveWorkbook(wb, file = file00,overwrite = TRUE)
    }
    else if (type == "points") {
        dirout = paste(getwd(),"/",sep = "")
        #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)
        PDF=paste(dirout,"PLSDA_p_",name,".pdf",sep="")
        PDF2=paste(dirout,"PLSDA_CV_",name,".pdf",sep="")

        wb <- createWorkbook()
        sheet1  <- addWorksheet(wb, sheetName = "Discrim")
        sheet2  <- addWorksheet(wb, sheetName = "MSEP")

        file00=paste(dirout, "PLSDA_",name,".xlsx",sep="")


        "biplot.color.P" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- dimnames(x)[[1]]
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                # text(x, xlabs, cex = cex[1], col = col, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "red", ...)
                if (var.axes)
                    ######################
                    # Farbe der Variablen-Pfeile
                    arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "red",
                           length = arrow.len)
                invisible()
            }


        "biplot.color_transp.P" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- dimnames(x)[[1]]
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                # text(x, xlabs, cex = cex[1], col = col, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "transparent", ...)
                if (var.axes)
                    ######################
                    # Farbe der Variablen-Pfeile
                    arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "transparent",
                           length = arrow.len)
                invisible()
            }

        #################################################################################

        ys=as.matrix(scale(yd,scale=FALSE))

        resst <- mvr(ys~dataS,ncomp=2,method="simpls")

        G=resst$scores
        H=resst$loadings

        rownames(G)=rownames(dataS)

        #Most discriminating vectors:
        lengthx=matrix(rep(0,nrow(H)),nrow=nrow(H))

        for(i in 1:nrow(H)){
            for(j in 1:nrow(H)){
                a=abs(H[i,1])
                b=abs(H[i,2])
                lengthx[i]=sqrt(a^2+b^2)
                rownames(lengthx)=rownames(H)
            }}

        discrim=as.matrix(lengthx[order(lengthx[,1],decreasing = TRUE),])

        matrixH=cbind(lengthx,H)
        matrixH=matrixH[order(matrixH[,1],decreasing = TRUE),]
        matrixH=matrixH[,-1]

        colnames(discrim)="Length"

        writeData(wb,sheet1,discrim,colNames = TRUE, rowNames = TRUE)

        pdf((PDF),width=10,height=10)
        par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

        number=length(which((Group[,3]==1)==TRUE))
        if (number!=0){
            biplot.color.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                         ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            biplot.color.P(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                         ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            biplot.color_transp.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                                ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
            abline(h=0,col="gray86")
            abline(v=0,col="gray86")
            text(H[,1],H[,2], labels = rownames(H),cex=0.8)

            dev.off()
        }else{
            if (ncol(dataS)>10){
            if (nrow(ys)<=20){
                ncompp=nrow(ys)-count
                seg=nrow(ys)
            } else {
                ncompp=15
                seg=10
            }

            misclass=mvr(ys~dataS,ncomp=ncompp,method="simpls",validation = "CV",segments=seg)

            MSEPy=MSEP(misclass,estimate="CV",intercept=FALSE)

            if (count == 2) {
                MSEPf=MSEPy$val[1,,1:ncompp]

                writeData(wb,sheet2,MSEPf,colNames = TRUE, rowNames = TRUE)
                }else{

                MSEPf=NULL
                for (i in 1:count){
                    MSEPff=MSEPy$val[1, i, 1:ncompp]
                    MSEPf=cbind(MSEPf,MSEPff)
                }

                colnames(MSEPf)=paste("Y",1:count,sep="")

                writeData(wb,sheet2,MSEPf,colNames = TRUE, rowNames = TRUE)
                }

            pdf((PDF2),width=10,height=10)
            plot(MSEPy,main="Cross validation - MSEP")
            #axis(1, at=1:15,labels = 1:15)
            dev.off()
            }

        biplot.color.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color.P(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color_transp.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                            ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        gr=as.factor(groups)
        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PLSDA - ",name),
                    xlab="PC1",ylab="PC2")
        legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.95, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PLSDA - ",name),
                    xlab="PC1",ylab="PC2")
        legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

        plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
        abline(h=0,col="gray86")
        abline(v=0,col="gray86")
        text(H[,1],H[,2], labels = rownames(H),cex=0.8)

        dev.off()
        }
       saveWorkbook(wb, file = file00,overwrite = TRUE)
    }
    else if (type == "both") {
        dirout = paste(getwd(),"/",sep = "")
        #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)
        PDF=paste(dirout,"PLSDA_b_",name,".pdf",sep="")
        PDF2=paste(dirout,"PLSDA_CV_",name,".pdf",sep="")

        wb <- createWorkbook()
        sheet1  <- addWorksheet(wb, sheetName = "Discrim")
        sheet2  <- addWorksheet(wb, sheetName = "MSEP")

        file00=paste(dirout, "PLSDA_",name,".xlsx",sep="")

        "biplot.color.B" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- rownames(x)
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                text(x, xlabs, cex = cex[1], col = col,pos=3, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "red", ...)
                if (var.axes)
                    ######################
                    # Farbe der Variablen-Pfeile
                    arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "red",
                           length = arrow.len)
                invisible()
            }


        "biplot.color_transp.B" <-
            function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                      xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                      arrow.len = 0.1, pch=marks, pch.cex=1.2, xlabMY, ylabMY, ...)
            {
                n <- nrow(x)
                p <- nrow(y)
                if (is.null(xlabs)) {
                    xlabs <- rownames(x)
                    if (is.null(xlabs))
                        xlabs <- 1:n
                }
                xlabs <- as.character(xlabs)
                dimnames(x) <- list(xlabs, dimnames(x)[[2]])
                if (missing(ylabs)) {
                    ylabs <- dimnames(y)[[1]]
                    if (is.null(ylabs))
                        ylabs <- paste("Var", 1:p)
                }
                ylabs <- as.character(ylabs)
                dimnames(y) <- list(ylabs, dimnames(y)[[2]])
                if (length(cex) == 1)
                    cex <- c(cex, cex)
                if (missing(col)) {
                    col <- par("col")
                    if (!is.numeric(col))
                        col <- match(col, palette())
                    col <- c(col, col + 1)
                }
                else if (length(col) == 1)
                    col <- c(col, col)
                unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
                rangx1 <- unsigned.range(x[, 1])
                rangx2 <- unsigned.range(x[, 2])
                rangy1 <- unsigned.range(y[, 1])
                rangy2 <- unsigned.range(y[, 2])
                if (missing(xlim) && missing(ylim))
                    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
                else if (missing(xlim))
                    xlim <- rangx1
                else ylim <- rangx2
                ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
                on.exit(par(oldpar))
                oldpar <- par(pty = "s")
                ######################
                # Plot-Symbol definiert
                plot(x, type = "p", xlim = xlim, ylim = ylim, col = col,pch=pch,
                     cex=pch.cex, xlab = xlabMY, ylab = ylabMY, ...)
                # plot(x, type = "n", xlim = xlim, ylim = ylim, col = col, ...)
                ######################
                # kein Text bei Objekten
                text(x, xlabs, cex = cex[1], col = col,pos=3, ...)
                par(new = TRUE)
                plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
                         ratio, xlab = "", ylab = "", col = col, ...)
                axis(3, col = 1)
                axis(4, col = 1)
                box(col = 1)
                ######################
                # Text der Variablen
                text(y, labels = ylabs, cex = cex[2], col = "transparent", ...)
                if (var.axes)
                    ######################
                    # Farbe der Variablen-Pfeile
                    arrows(0, 0, y[, 1] * 0.8, y[, 2] * 0.8, col = "transparent",
                           length = arrow.len)
                invisible()
            }

        #################################################################################

        ys=as.matrix(scale(yd,scale=FALSE))

        resst <- mvr(ys~dataS,ncomp=2,method="simpls")

        G=resst$scores
        H=resst$loadings

        rownames(G)=rownames(dataS)

        #Mosr discriminating vectors:
        lengthx=matrix(rep(0,nrow(H)),nrow=nrow(H))

        for(i in 1:nrow(H)){
            for(j in 1:nrow(H)){
                a=abs(H[i,1])
                b=abs(H[i,2])
                lengthx[i]=sqrt(a^2+b^2)
                rownames(lengthx)=rownames(H)
            }}

        discrim=as.matrix(lengthx[order(lengthx[,1],decreasing = TRUE),])

        matrixH=cbind(lengthx,H)
        matrixH=matrixH[order(matrixH[,1],decreasing = TRUE),]
        matrixH=matrixH[,-1]

        colnames(discrim)="Length"

        writeData(wb,sheet1,discrim,colNames = TRUE, rowNames = TRUE)

        pdf((PDF),width=10,height=10)
        par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

        number=length(which((Group[,3]==1)==TRUE))
        if (number!=0){
            biplot.color.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                         ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            biplot.color.B(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                         ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            biplot.color_transp.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                                ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
            abline(h=0,col="gray86")
            abline(v=0,col="gray86")
            text(H[,1],H[,2], labels = rownames(H),cex=0.8)

            dev.off()
        }else{

            if (ncol(dataS)>10){
            if (nrow(ys)<=20){
                ncompp=nrow(ys)-count
                seg=nrow(ys)
            } else {
                ncompp=15
                seg=10
            }

            misclass=mvr(ys~dataS,ncomp=ncompp,method="simpls",validation = "CV",segments=seg)

            MSEPy=MSEP(misclass,estimate="CV",intercept=FALSE)

            if (count == 2) {
                MSEPf=MSEPy$val[1,,1:ncompp]
                writeData(wb,sheet2,MSEPf,colNames = TRUE, rowNames = TRUE)
                }else {
                MSEPf=NULL
                for (i in 1:count){
                    MSEPff=MSEPy$val[1, i, 1:ncompp]
                    MSEPf=cbind(MSEPf,MSEPff)
                }

                colnames(MSEPf)=paste("Y",1:count,sep="")

                writeData(wb,sheet2,MSEPf,colNames = TRUE, rowNames = TRUE)
                }

            pdf((PDF2),width=10,height=10)
            plot(MSEPy,main="Cross validation - MSEP")
            #axis(1, at=1:15,labels = 1:15)
            dev.off()
            }

        biplot.color.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color.B(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color_transp.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PLSDA - ",name, sep="")
                            ,xlabMY="PC1",ylabMY="PC2", pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        gr=as.factor(groups)
        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PLSDA - ",name),
                    xlab="PC1",ylab="PC2")
        legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.95, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PLSDA - ",name),
                    xlab="PC1",ylab="PC2")
        legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

        plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),xlab="PC1",ylab="PC2")
        abline(h=0,col="gray86")
        abline(v=0,col="gray86")
        text(H[,1],H[,2], labels = rownames(H),cex=0.8)

        dev.off()
        }
       saveWorkbook(wb, file = file00,overwrite = TRUE)
    }

}

