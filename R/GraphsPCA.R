#' Principal component analysis (PCA)
#'
#' Makes principal component analysis (PCA), displays score plots, loading plots, scree plots and biplots.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param type A type of plots must be defined by "points" (default), "names" or "both".
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param top How many rays with highest lengths should be visualised in biplot? The default is 20.
#' @param QCs logical. If TRUE (default) quality control samples (QCs) are automatically distinguished. See Details.
#' @param pairs logical. If TRUE (default) PCAs of all combinations of pairs of specific groups are displayed.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @details If quality control samples (QCs) are present in data and QCs=TRUE, versions with QCs and without them are displayed. If QCs=TRUE and QCs are not present in data, this step is automatically skipped.
#' @return Score plot, scree plot and biplot of PCA.
#' @return Excel file with lengths of rays in biplot and degree of class separation (DCS) of all possible pairs of groups in data.
#' @import car
#' @import openxlsx
#' @importFrom stats var
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @references Pierce, K.M. et al. (2005) Classification of gasoline data obtained by gas chromatography using a piecewise alignment algorithm combined with feature selection and principal component analysis, J CHROMATOGR A 1096, p. 101-110.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsPCA(data,name,groupnames)
#' @export
GraphsPCA=function(data,name,groupnames,type="points",tsf="clr",top=20,QCs=TRUE,pairs=TRUE){

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

        groupss=Group[,c(1,2)]
##########################################################################################################################
    if (QCs==TRUE){
    QC=grep("QC",rownames(dataM))
    if (length(QC)!=0){
         dataMall=dataM
         colorall=color
         marksall=marks
         groupsall=groups
         groupnamesall=groupnames

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
         Group=matrix(rep(NA,count*3),ncol=3)
         for (i in 1:count){
             Gr=grep(groupnames[i],rownames(dataM))
             Group[i,1]=min(Gr)
             Group[i,2]=max(Gr)
             Group[i,3]=length(Gr)
         }
         groupss=Group[,c(1,2)]

         dirout = paste(getwd(),"/",sep = "")
         PDFQC=paste(dirout,"PCA_QC_",name,".pdf",sep="")

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
         par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
         #hist(data, main = paste("Histogram - ",name, sep=""))
         L=svd(dataSall)$u[,c(1,2)]
         K=diag(svd(dataSall)$d[1:2])
         M=svd(dataSall)$v[,c(1,2)]
         G=(sqrt(nrow(dataSall)-1))*L
         H=(1/sqrt(nrow(dataSall)-1))*M%*%K

         rownames(G)=rownames(dataSall)
         rownames(H)=colnames(dataSall)

         Ge=svd(var(dataSall))$u
         Z=dataSall%*%Ge
         PC1=100*var(Z[,1])/sum(apply(Z,2,var))
         PC2=100*var(Z[,2])/sum(apply(Z,2,var))

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

         number=length(which((Group[,3]==1)==TRUE))
         if (number!=0){
             if (type=="points"){
             biplot.color.QC.P(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
             legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

             biplot.color.QC.P(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
             legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

             biplot.color_transp.QC.P(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                 ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                 ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
             legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

             plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                  ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                  ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
             abline(h=0,col="gray86")
             abline(v=0,col="gray86")
             text(H[,1],H[,2], labels = rownames(H),cex=0.8)
             dev.off()
             }
             if (type=="names"){
                 biplot.color.QC.N(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color.QC.N(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color_transp.QC.N(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),
                      ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                      ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 abline(h=0,col="gray86")
                 abline(v=0,col="gray86")
                 text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                 dev.off()
             }
             if (type=="both"){
                 biplot.color.QC.B(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color.QC.B(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color_transp.QC.B(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep=""),
                      ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                      ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 abline(h=0,col="gray86")
                 abline(v=0,col="gray86")
                 text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                 dev.off()
             }
         }else{
             if (type=="points"){
         biplot.color.QC.P(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                      ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                      ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
         legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

         biplot.color.QC.P(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                      ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                      ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
         legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

         biplot.color_transp.QC.P(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                             ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                             ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
         legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

         gr=as.factor(groupsall)
         dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PCA - ",name),
                     xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
         legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

         dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.95, plot.points=TRUE, cex=1.2,center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PCA - ",name),
                     xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
         legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

         plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
              ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
              ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
         abline(h=0,col="gray86")
         abline(v=0,col="gray86")
         text(H[,1],H[,2], labels = rownames(H),cex=0.8)
         dev.off()
             }
             if (type=="names"){
                 biplot.color.QC.N(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color.QC.N(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color_transp.QC.N(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 gr=as.factor(groupsall)
                 dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PCA - ",name),
                             xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                 dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.95, plot.points=TRUE, cex=1.2,center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PCA - ",name),
                             xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                 plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                      ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                      ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 abline(h=0,col="gray86")
                 abline(v=0,col="gray86")
                 text(H[,1],H[,2], labels = rownames(H),cex=0.8)

                 dev.off()
             }
             if (type=="both"){
                 biplot.color.QC.B(G,H,cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color.QC.B(G,matrixH[1:top,],cex=c(1.2,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 biplot.color_transp.QC.B(G,H,cex=c(0.8,0.6),col=c(colorall),main = paste("PCA - ",name, sep="")
                                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marksall,lwd=0.00000001,arrow.len=0.05)
                 legend("topleft",legend = groupnamesall, pch = unique(marksall), col = unique(colorall),cex=0.6)

                 gr=as.factor(groupsall)
                 dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PCA - ",name),
                             xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                 dataEllipse(G[,1],G[,2],gr,group.labels=groupnamesall, levels=0.95, plot.points=TRUE, cex=1.2,center.cex=0.2, col=unique(colorall), pch=unique(marksall),center.pch=3,main=paste("PCA - ",name),
                             xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 legend("bottomright",legend = groupnamesall, pch = unique(marksall), col = unique(colorall), cex=0.6)

                 plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                      ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                      ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                 abline(h=0,col="gray86")
                 abline(v=0,col="gray86")
                 text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                 dev.off()
             }
         }

    }
    }

        if (tsf=="pareto"){
         dataS=scale(dataM, scale=TRUE, center=TRUE)
        }

        if (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN"|tsf=="lnPQN"){
            dataS=scale(dataM, scale=FALSE, center=TRUE)
        }

        if (tsf=="none"){
            dataS=as.matrix(dataM)
        }

    if (type == "names") {
    dirout = paste(getwd(),"/",sep = "")
    #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)
    PDF=paste(dirout,"PCA_n_",name,".pdf",sep="")
    PDF2=paste(dirout,"PCA_scree_",name,".pdf",sep="")

    wb <- createWorkbook()
    sheet1  <- addWorksheet(wb, sheetName = "Discrim")
    sheet2  <- addWorksheet(wb, sheetName = "DCS")

    file00=paste(dirout, "PCA_",name,".xlsx",sep="")

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

    pdf((PDF),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
    #hist(dataS, main = paste("Histogram - ",name, sep=""))
    L=svd(dataS)$u[,c(1,2)]
    K=diag(svd(dataS)$d[1:2])
    M=svd(dataS)$v[,c(1,2)]
    G=(sqrt(nrow(dataS)-1))*L
    H=(1/sqrt(nrow(dataS)-1))*M%*%K

    rownames(G)=rownames(dataS)
    rownames(H)=colnames(dataS)

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

    #cloud=NULL
    #for(i in 1:length(groupnames)){
    #    cloud1=matrix(rep(0,2),ncol=2)
    #    cloud1[1,1]=median(G[(groupss[i,1]:groupss[i,2]),1])
    #    cloud1[1,2]=median(G[(groupss[i,1]:groupss[i,2]),2])
    #    cloud=rbind(cloud,cloud1)
    #}

    #rownames(cloud)=c(paste("Angle_",groupnames,sep=""))

    #angle=matrix(rep(0,nrow(cloud)*nrow(H)),ncol=nrow(cloud))
    #theta=matrix(rep(0,nrow(cloud)*nrow(H)),ncol=nrow(cloud))
    #colnames(angle)=paste("Angle_",groupnames,sep="")
    #rownames(angle)=rownames(H)

    #for(i in 1:nrow(cloud)){
    #    for(j in 1:nrow(H)){
    #        a=cloud[i,]
    #        b=H[j,1:2]
    #        theta[j,i]= acos(sum(a*b) / (sqrt(sum(a * a)) * sqrt(sum(b * b))))
    #        angle[j,i]=theta[j,i]*180/pi}}    #prevod na stupne
    # angle

    colnames(discrim)="Length"
    #datafile=cbind(discrim,angle)

    writeData(wb,sheet1,discrim,colNames = TRUE, rowNames = TRUE)

    # Scree plot
    Ge=svd(var(dataS))$u
    Z=dataS%*%Ge

    PC1=100*var(Z[,1])/sum(apply(Z,2,var))
    PC2=100*var(Z[,2])/sum(apply(Z,2,var))

    if (ncol(dataS)>10){
    PC3=100*var(Z[,3])/sum(apply(Z,2,var))
    PC4=100*var(Z[,4])/sum(apply(Z,2,var))
    PC5=100*var(Z[,5])/sum(apply(Z,2,var))
    PC6=100*var(Z[,6])/sum(apply(Z,2,var))
    PC7=100*var(Z[,7])/sum(apply(Z,2,var))
    PC8=100*var(Z[,8])/sum(apply(Z,2,var))
    PC9=100*var(Z[,9])/sum(apply(Z,2,var))
    PC10=100*var(Z[,10])/sum(apply(Z,2,var))

    PC=c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
    PC

    cumulative=c(PC1,PC1+PC2,PC1+PC2+PC3,PC1+PC2+PC3+PC4,PC1+PC2+PC3+PC4+PC5,PC1+PC2+PC3+PC4+PC5+PC6,PC1+PC2+PC3+PC4+PC5+PC6+PC7,
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9,sum(PC))
    cumulative

    pdf(PDF2,width=10,height=10)
    b=barplot(PC,width=1,ylim=c(0,100),names.arg=1:10,col="grey80",xlab="Principal component",ylab="Variance explained (%)",cex.lab=1.5,
              cex.axis=1.5,cex.names=1.5,main=paste("Scree plot - ",name))
    lines(b,cumulative,col="black",lwd=1.5,type="b",pch=16)
    abline(h=100)
    abline(h=0)
    abline(v=12.45)
    text(b[1],(PC1+2.5),round(PC1,1))
    text(b[2],(PC2+2.5),round(PC2,1))
    text(b[3],(PC3+2.5),round(PC3,1))
    text(b[4],(PC4+2.5),round(PC4,1))
    text(b[5],(PC5+2.5),round(PC5,1))
    text(b[6],(PC6+2.5),round(PC6,1))
    text(b[7],(PC7+2.5),round(PC7,1))
    text(b[8],(PC8+2.5),round(PC8,1))
    text(b[9],(PC9+2.5),round(PC9,1))
    text(b[10],(PC10+2.5),round(PC10,1))
    legend("topleft", legend = c("Absolute","Cumulative"),col = c("black", "black"),lwd=c(NA,1.5),lty = c(NA, 1),pch = c(22, 16),
           pt.bg = c("grey80",NA),pt.cex = c(2,1))
    dev.off()
    }

    number=length(which((Group[,3]==1)==TRUE))
    if (number!=0){
        biplot.color.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                     ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=groups,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

        biplot.color.N(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                     ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=groups,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

        biplot.color_transp.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                            ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                            ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=groups,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

        plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
             ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
             ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
        abline(h=0,col="gray86")
        abline(v=0,col="gray86")
        text(H[,1],H[,2], labels = rownames(H),cex=0.8)
        dev.off()
    }else{
    biplot.color.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                 ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                 ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=groups,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

    biplot.color.N(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                 ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                 ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=groups,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

    biplot.color_transp.N(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                        ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                        ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=groups,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = 19, col = unique(color),cex=0.6)

    gr=as.factor(groups)
    dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=rep(19,length(groupnames)),center.pch=3,main=paste("PCA - ",name),
                xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
    legend("bottomright",legend = groupnames, pch = 19, col = unique(color), cex=0.6)

    dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.95, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=rep(19,length(groupnames)),center.pch=3,main=paste("PCA - ",name),
                xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
    legend("bottomright",legend = groupnames, pch = 19, col = unique(color), cex=0.6)

    plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
         ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
         ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
    abline(h=0,col="gray86")
    abline(v=0,col="gray86")
    text(H[,1],H[,2], labels = rownames(H),cex=0.8)
    dev.off()

    ####################################################################################################
    count=length(unique(color))

    mmean=matrix(c(rep(0,count*2)),nrow=count)
    dev=matrix(c(rep(0,count)),nrow=count)
    for(i in 1:count){
        scores=G[(groupss[i,1]:groupss[i,2]),]
        mmean[i,1]=mean(scores[,1])
        mmean[i,2]=mean(scores[,2])
        dev[i,]=var(dist(rbind(mmean[i,],scores),method = "euclidean")[1:nrow(scores)])
    }

    DCS=NULL
    for(i in 1:(count-1)){
        for(j in (i+1):count){
            d=dist(rbind(mmean[i,],mmean[j,]),method = "euclidean")[1]
            DCS2=d/(sqrt(dev[i,]+dev[j,]))
            DCS=cbind(DCS,DCS2)
        }
    }

    nazvy=NULL
    for(i in 1:(count-1)){
        for(j in (i+1):count){
            col=paste(groupnames[i],"_",groupnames[j],sep="")
            nazvy=c(nazvy,col)
        }
    }

    colnames(DCS)=nazvy

    writeData(wb,sheet2,DCS,colNames = TRUE, rowNames = TRUE)
    }
    saveWorkbook(wb,file = file00,overwrite = TRUE)

    if(pairs == TRUE & count > 2){
        dirout2 = paste(dirout,"/","PCA",sep = "")
        dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
        setwd(dirout2)
        for (k in 1:(count-1)){
            for (l in (k+1):count){
                PDFpair=paste("PCA_",name,"_",groupnames[k],"_",groupnames[l],".pdf",sep="")
                vyber=c(Group[k,1]:Group[k,2],Group[l,1]:Group[l,2])

                datavyb0=dataM[vyber,]
                if (tsf=="pareto"){
                    datavyb=scale(datavyb0, scale=TRUE, center=TRUE)
                }

                if (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN"|tsf=="lnPQN"){
                    datavyb=scale(datavyb0, scale=FALSE, center=TRUE)
                }

                if (tsf=="none"){
                    datavyb=as.matrix(datavyb0)
                }
                L=svd(datavyb)$u[,c(1,2)]
                K=diag(svd(datavyb)$d[1:2])
                M=svd(datavyb)$v[,c(1,2)]
                G=(sqrt(nrow(datavyb)-1))*L
                H=(1/sqrt(nrow(datavyb)-1))*M%*%K

                rownames(G)=rownames(datavyb)
                rownames(H)=colnames(datavyb)

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

                Ge=svd(var(datavyb))$u
                Z=datavyb%*%Ge
                PC1=100*var(Z[,1])/sum(apply(Z,2,var))
                PC2=100*var(Z[,2])/sum(apply(Z,2,var))

                pdf((PDFpair),width=10,height=10)
                biplot.color.N(G,H,cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                             ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                             ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                biplot.color.N(G,matrixH[1:top,],cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                             ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                             ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                biplot.color_transp.N(G,H,cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                                    ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                    ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                     ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)

                dev.off()

            }
        }
        setwd(dirout)
    }

    }
    else if (type == "points") {
    dirout = paste(getwd(),"/",sep = "")
    #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)
    PDF=paste(dirout,"PCA_p_",name,".pdf",sep="")
    PDF2=paste(dirout,"PCA_scree_",name,".pdf",sep="")

    wb <- createWorkbook()
    sheet1  <- addWorksheet(wb, sheetName = "Discrim")
    sheet2  <- addWorksheet(wb, sheetName = "DCS")

    file00=paste(dirout, "PCA_",name,".xlsx",sep="")


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

    pdf((PDF),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
    #hist(dataS, main = paste("Histogram - ",name, sep=""))
    L=svd(dataS)$u[,c(1,2)]
    K=diag(svd(dataS)$d[1:2])
    M=svd(dataS)$v[,c(1,2)]
    G=(sqrt(nrow(dataS)-1))*L
    H=(1/sqrt(nrow(dataS)-1))*M%*%K

    rownames(G)=rownames(dataS)
    rownames(H)=colnames(dataS)


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

    # Scree plot
    Ge=svd(var(dataS))$u
    Z=dataS%*%Ge

    PC1=100*var(Z[,1])/sum(apply(Z,2,var))
    PC2=100*var(Z[,2])/sum(apply(Z,2,var))

    if (ncol(dataS)>10){
    PC3=100*var(Z[,3])/sum(apply(Z,2,var))
    PC4=100*var(Z[,4])/sum(apply(Z,2,var))
    PC5=100*var(Z[,5])/sum(apply(Z,2,var))
    PC6=100*var(Z[,6])/sum(apply(Z,2,var))
    PC7=100*var(Z[,7])/sum(apply(Z,2,var))
    PC8=100*var(Z[,8])/sum(apply(Z,2,var))
    PC9=100*var(Z[,9])/sum(apply(Z,2,var))
    PC10=100*var(Z[,10])/sum(apply(Z,2,var))

    PC=c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
    PC

    cumulative=c(PC1,PC1+PC2,PC1+PC2+PC3,PC1+PC2+PC3+PC4,PC1+PC2+PC3+PC4+PC5,PC1+PC2+PC3+PC4+PC5+PC6,PC1+PC2+PC3+PC4+PC5+PC6+PC7,
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9,sum(PC))
    cumulative

    colnames(discrim)="Length"

    writeData(wb,sheet1,discrim,colNames = TRUE, rowNames = TRUE)

    pdf(PDF2,width=10,height=10)
    b=barplot(PC,width=1,ylim=c(0,100),names.arg=1:10,col="grey80",xlab="Principal component",ylab="Variance explained (%)",cex.lab=1.5,
              cex.axis=1.5,cex.names=1.5,main=paste("Scree plot - ",name))
    lines(b,cumulative,col="black",lwd=1.5,type="b",pch=16)
    abline(h=100)
    abline(h=0)
    abline(v=12.45)
    text(b[1],(PC1+2.5),round(PC1,1))
    text(b[2],(PC2+2.5),round(PC2,1))
    text(b[3],(PC3+2.5),round(PC3,1))
    text(b[4],(PC4+2.5),round(PC4,1))
    text(b[5],(PC5+2.5),round(PC5,1))
    text(b[6],(PC6+2.5),round(PC6,1))
    text(b[7],(PC7+2.5),round(PC7,1))
    text(b[8],(PC8+2.5),round(PC8,1))
    text(b[9],(PC9+2.5),round(PC9,1))
    text(b[10],(PC10+2.5),round(PC10,1))
    legend("topleft", legend = c("Absolute","Cumulative"),col = c("black", "black"),lty = c(NA, 1),pch = c(22, 16),
           pt.bg = c("grey80",NA),pt.cex = c(2,1))
    dev.off()
    }

    number=length(which((Group[,3]==1)==TRUE))
    if (number!=0){
        biplot.color.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                     ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color.P(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                     ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

         biplot.color_transp.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                            ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                            ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
             ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
             ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
        abline(h=0,col="gray86")
        abline(v=0,col="gray86")
        text(H[,1],H[,2], labels = rownames(H),cex=0.8)

        dev.off()
    } else{
    biplot.color.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                 ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                 ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

    biplot.color.P(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                 ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                 ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

    biplot.color_transp.P(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                        ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                        ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
    legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

    gr=as.factor(groups)
    dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PCA - ",name),
                xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
    legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

    dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.95, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PCA - ",name),
                xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
    legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

    plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
         ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
         ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
    abline(h=0,col="gray86")
    abline(v=0,col="gray86")
    text(H[,1],H[,2], labels = rownames(H),cex=0.8)
    dev.off()

    ####################################################################################################
    count=length(unique(color))

    mmean=matrix(c(rep(0,count*2)),nrow=count)
    dev=matrix(c(rep(0,count)),nrow=count)
    for(i in 1:count){
        scores=G[(groupss[i,1]:groupss[i,2]),]
        mmean[i,1]=mean(scores[,1])
        mmean[i,2]=mean(scores[,2])
        dev[i,]=var(dist(rbind(mmean[i,],scores),method = "euclidean")[1:nrow(scores)])
    }

    DCS=NULL
    for(i in 1:(count-1)){
        for(j in (i+1):count){
            d=dist(rbind(mmean[i,],mmean[j,]),method = "euclidean")[1]
            DCS2=d/(sqrt(dev[i,]+dev[j,]))
            DCS=cbind(DCS,DCS2)
        }
    }

    nazvy=NULL
    for(i in 1:(count-1)){
        for(j in (i+1):count){
            col=paste(groupnames[i],"_",groupnames[j],sep="")
            nazvy=c(nazvy,col)
        }
    }

    colnames(DCS)=nazvy

    writeData(wb,sheet2,DCS,colNames = TRUE, rowNames = TRUE)
    }
    saveWorkbook(wb,file = file00,overwrite = TRUE)

    if(pairs == TRUE & count > 2){
        dirout2 = paste(dirout,"/","PCA",sep = "")
        dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
        setwd(dirout2)
        for (k in 1:(count-1)){
            for (l in (k+1):count){
                PDFpair=paste("PCA_",name,"_",groupnames[k],"_",groupnames[l],".pdf",sep="")
                vyber=c(Group[k,1]:Group[k,2],Group[l,1]:Group[l,2])

                datavyb0=dataM[vyber,]
                if (tsf=="pareto"){
                    datavyb=scale(datavyb0, scale=TRUE, center=TRUE)
                }

                if (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN"|tsf=="lnPQN"){
                    datavyb=scale(datavyb0, scale=FALSE, center=TRUE)
                }

                if (tsf=="none"){
                    datavyb=as.matrix(datavyb0)
                }
                L=svd(datavyb)$u[,c(1,2)]
                K=diag(svd(datavyb)$d[1:2])
                M=svd(datavyb)$v[,c(1,2)]
                G=(sqrt(nrow(datavyb)-1))*L
                H=(1/sqrt(nrow(datavyb)-1))*M%*%K

                rownames(G)=rownames(datavyb)
                rownames(H)=colnames(datavyb)

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

                Ge=svd(var(datavyb))$u
                Z=datavyb%*%Ge
                PC1=100*var(Z[,1])/sum(apply(Z,2,var))
                PC2=100*var(Z[,2])/sum(apply(Z,2,var))

                pdf((PDFpair),width=10,height=10)
                biplot.color.P(G,H,cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                               ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                               ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                biplot.color.P(G,matrixH[1:top,],cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                               ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                               ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                biplot.color_transp.P(G,H,cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                                      ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                      ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                     ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                abline(h=0,col="gray86")
                abline(v=0,col="gray86")
                text(H[,1],H[,2], labels = rownames(H),cex=0.8)
                dev.off()

            }
        }
        setwd(dirout)
    }
}

    if (type == "both") {
        dirout = paste(getwd(),"/",sep = "")
        #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)
        PDF=paste(dirout,"PCA_b_",name,".pdf",sep="")
        PDF2=paste(dirout,"PCA_scree_",name,".pdf",sep="")

        wb <- createWorkbook()
        sheet1  <- addWorksheet(wb, sheetName = "Discrim")
        sheet2  <- addWorksheet(wb, sheetName = "DCS")

        file00=paste(dirout, "PCA_",name,".xlsx",sep="")


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

        pdf((PDF),width=10,height=10)
        par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
        #hist(dataS, main = paste("Histogram - ",name, sep=""))
        L=svd(dataS)$u[,c(1,2)]
        K=diag(svd(dataS)$d[1:2])
        M=svd(dataS)$v[,c(1,2)]
        G=(sqrt(nrow(dataS)-1))*L
        H=(1/sqrt(nrow(dataS)-1))*M%*%K

        rownames(G)=rownames(dataS)
        rownames(H)=colnames(dataS)

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

        # Scree plot
        Ge=svd(var(dataS))$u
        Z=dataS%*%Ge

        PC1=100*var(Z[,1])/sum(apply(Z,2,var))
        PC2=100*var(Z[,2])/sum(apply(Z,2,var))

        if (ncol(dataS)>10){
        PC3=100*var(Z[,3])/sum(apply(Z,2,var))
        PC4=100*var(Z[,4])/sum(apply(Z,2,var))
        PC5=100*var(Z[,5])/sum(apply(Z,2,var))
        PC6=100*var(Z[,6])/sum(apply(Z,2,var))
        PC7=100*var(Z[,7])/sum(apply(Z,2,var))
        PC8=100*var(Z[,8])/sum(apply(Z,2,var))
        PC9=100*var(Z[,9])/sum(apply(Z,2,var))
        PC10=100*var(Z[,10])/sum(apply(Z,2,var))

        PC=c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
        PC

        cumulative=c(PC1,PC1+PC2,PC1+PC2+PC3,PC1+PC2+PC3+PC4,PC1+PC2+PC3+PC4+PC5,PC1+PC2+PC3+PC4+PC5+PC6,PC1+PC2+PC3+PC4+PC5+PC6+PC7,
                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9,sum(PC))
        cumulative

        pdf(PDF2,width=10,height=10)
        b=barplot(PC,width=1,ylim=c(0,100),names.arg=1:10,col="grey80",xlab="Principal component",ylab="Variance explained (%)",cex.lab=1.5,
                  cex.axis=1.5,cex.names=1.5,main=paste("Scree plot - ",name))
        lines(b,cumulative,col="black",lwd=1.5,type="b",pch=16)
        abline(h=100)
        abline(h=0)
        abline(v=12.45)
        text(b[1],(PC1+2.5),round(PC1,1))
        text(b[2],(PC2+2.5),round(PC2,1))
        text(b[3],(PC3+2.5),round(PC3,1))
        text(b[4],(PC4+2.5),round(PC4,1))
        text(b[5],(PC5+2.5),round(PC5,1))
        text(b[6],(PC6+2.5),round(PC6,1))
        text(b[7],(PC7+2.5),round(PC7,1))
        text(b[8],(PC8+2.5),round(PC8,1))
        text(b[9],(PC9+2.5),round(PC9,1))
        text(b[10],(PC10+2.5),round(PC10,1))
        legend("topleft", legend = c("Absolute","Cumulative"),col = c("black", "black"),lty = c(NA, 1),pch = c(22, 16),
               pt.bg = c("grey80",NA),pt.cex = c(2,1))
        dev.off()
        }

        colnames(discrim)="Length"

        writeData(wb,sheet1,discrim,colNames = TRUE, rowNames = TRUE)

        number=length(which((Group[,3]==1)==TRUE))
        if (number!=0){
            biplot.color.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                         ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                         ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            biplot.color.B(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                         ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                         ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            biplot.color_transp.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                                ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
            legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

            plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                 ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                 ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
            abline(h=0,col="gray86")
            abline(v=0,col="gray86")
            text(H[,1],H[,2], labels = rownames(H),cex=0.8)

            dev.off()
        }else{
        biplot.color.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                     ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color.B(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                     ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                     ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        biplot.color_transp.B(G,H,cex=c(0.8,0.6),col=c(color),main = paste("PCA - ",name, sep="")
                            ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                            ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames, pch = unique(marks), col = unique(color),cex=0.6)

        gr=as.factor(groups)
        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.75, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PCA - ",name),
                    xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
        legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames, levels=0.95, plot.points=TRUE,cex=1.2, center.cex=0.2, col=unique(color), pch=unique(marks),center.pch=3,main=paste("PCA - ",name),
                    xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%"),ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
        legend("bottomright",legend = groupnames, pch = unique(marks), col = unique(color), cex=0.6)

        plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
             ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
             ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
        abline(h=0,col="gray86")
        abline(v=0,col="gray86")
        text(H[,1],H[,2], labels = rownames(H),cex=0.8)

        dev.off()

        ####################################################################################################
        count=length(unique(color))

        mmean=matrix(c(rep(0,count*2)),nrow=count)
        dev=matrix(c(rep(0,count)),nrow=count)
        for(i in 1:count){
            scores=G[(groupss[i,1]:groupss[i,2]),]
            mmean[i,1]=mean(scores[,1])
            mmean[i,2]=mean(scores[,2])
            dev[i,]=var(dist(rbind(mmean[i,],scores),method = "euclidean")[1:nrow(scores)])
        }

        DCS=NULL
        for(i in 1:(count-1)){
            for(j in (i+1):count){
                d=dist(rbind(mmean[i,],mmean[j,]),method = "euclidean")[1]
                DCS2=d/(sqrt(dev[i,]+dev[j,]))
                DCS=cbind(DCS,DCS2)
            }
        }

        nazvy=NULL
        for(i in 1:(count-1)){
            for(j in (i+1):count){
                col=paste(groupnames[i],"_",groupnames[j],sep="")
                nazvy=c(nazvy,col)
            }
        }

        colnames(DCS)=nazvy

        writeData(wb,sheet2,DCS,colNames = TRUE, rowNames = TRUE)
        }
        saveWorkbook(wb,file = file00,overwrite = TRUE)

        if(pairs == TRUE & count > 2){
            dirout2 = paste(dirout,"/","PCA",sep = "")
            dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
            setwd(dirout2)
            for (k in 1:(count-1)){
                for (l in (k+1):count){
                    PDFpair=paste("PCA_",name,"_",groupnames[k],"_",groupnames[l],".pdf",sep="")
                    vyber=c(Group[k,1]:Group[k,2],Group[l,1]:Group[l,2])

                    datavyb0=dataM[vyber,]
                    if (tsf=="pareto"){
                        datavyb=scale(datavyb0, scale=TRUE, center=TRUE)
                    }

                    if (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN"|tsf=="lnPQN"){
                        datavyb=scale(datavyb0, scale=FALSE, center=TRUE)
                    }

                    if (tsf=="none"){
                        datavyb=as.matrix(datavyb0)
                    }
                    L=svd(datavyb)$u[,c(1,2)]
                    K=diag(svd(datavyb)$d[1:2])
                    M=svd(datavyb)$v[,c(1,2)]
                    G=(sqrt(nrow(datavyb)-1))*L
                    H=(1/sqrt(nrow(datavyb)-1))*M%*%K

                    rownames(G)=rownames(datavyb)
                    rownames(H)=colnames(datavyb)

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

                    Ge=svd(var(datavyb))$u
                    Z=datavyb%*%Ge
                    PC1=100*var(Z[,1])/sum(apply(Z,2,var))
                    PC2=100*var(Z[,2])/sum(apply(Z,2,var))

                    pdf((PDFpair),width=10,height=10)
                    biplot.color.B(G,H,cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                    legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                    biplot.color.B(G,matrixH[1:top,],cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                                   ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                   ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                    legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                    biplot.color_transp.B(G,H,cex=c(0.8,0.6),col=color[vyber],main = paste("PCA - ",name, sep="")
                                          ,xlabMY=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                                          ,ylabMY=paste("PC2 - ", round(PC2, digits = 2),"%"), pch=marks[vyber],lwd=0.00000001,arrow.len=0.05)
                    legend("topleft",legend = groupnames[c(k,l)], pch = unique(marks[vyber]), col = unique(color[vyber]),cex=0.6)

                    plot(H[,1],H[,2],col="black",type="n",main = paste("Loading plot - ",name, sep="")
                         ,xlab=paste("PC1 - ", round(PC1, digits = 2), "%; Cumulative = ",round(PC1, digits = 2)+round(PC2, digits = 2),"%")
                         ,ylab=paste("PC2 - ", round(PC2, digits = 2),"%"))
                    abline(h=0,col="gray86")
                    abline(v=0,col="gray86")
                    text(H[,1],H[,2], labels = rownames(H),cex=0.8)

                    dev.off()

                }
            }
            setwd(dirout)
        }
    }

}

