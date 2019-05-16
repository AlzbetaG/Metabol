#' Variable importance plots (VIP)
#'
#' Makes variable importance plots (VIP) of partial least squares - discriminant analysis (PLS-DA), also displays score plots and biplots.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param top How many most important variables should be in zoomed VIP plot and biplot? The default is 30.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details VIP plot can be used only for comparison of two groups. If there is more groups in data, all possible combinations of pairs are evaluated.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return VIP plots, score plots and biplots of PLS-DA.
#' @return Excel file with mean, sd and angles of rays from group center for every variable.
#' @import car
#' @importFrom openxlsx write.xlsx
#' @importFrom stats sd quantile
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @references Wold, S. et al. (2001) PLS-regression: a basic tool of chemometrics, CHEMOLAB 58 (2), p. 109-130.
#' @references Chong, I.G., Jun, C.H. (2005) Performance of some variable selection methods when multicollinearity is present, CHEMOLAB 78, p. 103-112.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsVIP(data,name,groupnames)
#' @export
GraphsVIP=function(data,name,groupnames,tsf="clr",top=30,QCs=FALSE){

    ############################################################################################
    dirout=getwd()
    dirout2 = paste(dirout,"/","PLSDA",sep = "")
    dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)        # vytvori v zadane ceste slozku s nazvem "note"
    setwd(dirout2)

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
    ############################################################################################
    #VIP

    VIP <- function(object) {
        if (object$method != "oscorespls")
            stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
        if (nrow(object$Yloadings) > 1)
            stop("Only implemented for single-response models")

        SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
        Wnorm2 <- colSums(object$loading.weights^2)
        SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
        sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
    }

    ############################################################################################
    "biplot.color" <-
        function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                  xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                  arrow.len = 0.1, pch=marks, pch.cex=1, xlabMY, ylabMY, ...)
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


    "biplot.color_transp" <-
        function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                  xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                  arrow.len = 0.1, pch=marks, pch.cex=1, xlabMY, ylabMY, ...)
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

    ############################################################################################

    for (i in 1:(count-1)){
        for (j in (i+1):count){
            PDF2=paste("VIP_all_",name,"_",groupnames[i],"_",groupnames[j],".pdf",sep="")
            PDF3=paste("VIP_top_",name,"_",groupnames[i],"_",groupnames[j],".pdf",sep="")
            PDF4=paste("PLSDA_",name,"_",groupnames[i],"_",groupnames[j],".pdf",sep="")

        vyber=c(Group[i,1]:Group[i,2],Group[j,1]:Group[j,2])
        groups2=c(rep(0,Group[i,3]),rep(1,Group[j,3]))
        lc=Group[i,3]
        lp=Group[j,3]
        datanew=dataM[vyber,]
        dataSnew=scale(datanew, scale=FALSE, center=TRUE)
        y1=rep(1,lc)
        y2=rep(0,lp)
        yyd=c(y1,y2)
        color2=color[vyber]
        marks2=marks[vyber]
        groups2=c(rep(1,lc),rep(2,lp))
        groupnames2=groupnames[c(i,j)]

        ############################################################################################
        #The distance of rays from the center:
        ys=as.matrix(scale(yyd,scale=FALSE))
        resst2 <- mvr(ys~dataSnew,ncomp=2,method="simpls")
        G=resst2$scores
        H=resst2$loadings
        rownames(G)=rownames(dataSnew)

        lengthx=matrix(rep(0,nrow(H)),nrow=nrow(H))

        for(k in 1:nrow(H)){
            for(l in 1:nrow(H)){
                a=abs(H[k,1])
                b=abs(H[k,2])
                lengthx[k]=sqrt(a^2+b^2)
                rownames(lengthx)=rownames(H)
            }}

        discrim=as.matrix(lengthx[order(lengthx[,1],decreasing = TRUE),])

        matrixH=cbind(lengthx,H)
        matrixH=matrixH[order(matrixH[,1],decreasing = TRUE),]
        matrixH=matrixH[,-1]

        sc1=quantile(G[1:lc,1],0.5)
        sc2=quantile(G[1:lc,2],0.5)
        sc=c(sc1,sc2)
        sp1=quantile(G[(lc+1):(lc+lp),1],0.5)
        sp2=quantile(G[(lc+1):(lc+lp),2],0.5)
        sp=c(sp1,sp2)

        colorc=unique(color2)[1]
        colorp=unique(color2)[2]

        theta <- function(x1, x2) acos( sum(x1*x2) / ( sqrt(sum(x1 * x1)) * sqrt(sum(x2 * x2)) ) )

        uhel=matrix(rep(0,nrow(H)*2),ncol=2)
        rownames(uhel)=rownames(H)
        colnames(uhel)=c(groupnames2[1],groupnames2[2])

        for(m in 1:nrow(H)){
            uhel[m,1]=theta(sc,H[m,])
            uhel[m,2]=theta(sp,H[m,])
        }

        ###############################################################

        ys=as.matrix(scale(yyd,scale=FALSE))

        if (length(ys)<=10){ncompp=length(ys)-1
        } else if (ncol(dataSnew)<10) {ncompp=ncol(dataSnew)-1
        } else {ncompp=10
        }

        resst <- mvr(ys~dataSnew,ncomp=ncompp,method="oscorespls")

        msep=MSEP(resst,intercept=FALSE)
        comp=which.min(msep$val[1,,])

        if (comp>(ncompp-2)){comp=ncompp-2
        } else {
            comp=comp
        }

        ###############################################################

        vip=VIP(resst)[comp,]

        vipmean=matrix(vip,ncol=ncol(dataSnew))
        colnames(vipmean)=colnames(dataSnew)

        ycon=rep(1,length(which(yyd==1)))
        ycon2=rep(0,length(which(yyd==0)))

        BE=matrix(rep(NA,length(ycon)*length(ycon2)*ncol(dataSnew)),ncol=ncol(dataSnew))
        #jackknife for standard deviation estimation
        k=1
        for(m in 1:length(ycon)){
            for(n in 1:length(ycon2)){
                testX=dataSnew[-c(m,length(ycon)+n),]
                testy=c(ycon[-m],ycon2[-n])
                testys=as.matrix(scale(testy,scale=FALSE))
                resst <- mvr(testys~testX,ncomp=comp,method="oscorespls")
                BE[k,]=VIP(resst)[comp,]
                k=k+1
            }
        }

        sdbeta=matrix(apply(BE,2,sd),ncol=ncol(dataSnew))
        matice=rbind(vipmean,sdbeta,t(uhel))

        matice=matice[,order(matice[1,],decreasing =TRUE)]

        vipmean=matice[1,]
        sdbeta=matice[2,]
        uhelcon=matice[3,]
        uhelpac=matice[4,]

        maxy=max(vipmean)+1.75*max(sdbeta)

        ccolor=rep(0,ncol(matice))
        for(l in 1:ncol(matice)){
            if (uhelcon[l]<uhelpac[l]){ccolor[l]=colorc}
            else ccolor[l]=colorp
        }

        pdf((PDF2),height=10,width=round((length(colnames(data))/8),digits = 0))
        bar=barplot(vipmean,ylim=c(0,maxy),main=paste("VIP plot - ",name, sep=""), ylab="VIP",col=ccolor,cex.names=0.75,las=2)
        segments(bar, vipmean - sdbeta, bar, vipmean + sdbeta, lwd=1)
        segments(bar - 0.1, vipmean - sdbeta, bar + 0.1, vipmean - sdbeta, lwd=1)
        segments(bar - 0.1, vipmean + sdbeta, bar + 0.1, vipmean + sdbeta, lwd=1)
        legend("topright",legend = groupnames2, pch = 16, col = unique(color2))
        dev.off()

        rownames(matice)=c("Mean","Sd",paste("Angle_",groupnames[i],sep=""),paste("Angle_",groupnames[j],sep=""))

        write.xlsx(t(matice),file = paste("VIP_",name,"_",groupnames[i],"_",groupnames[j],".xlsx", sep=""),
                   sheetName=paste(groupnames[i],"_",groupnames[j],sep=""),
                   col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)

        vipmean=vipmean[1:top]
        sdbeta=sdbeta[1:top]

        pdf((PDF3),width=15,height=10)
        par(mar=c(11,4.5,1,0))
        bar=barplot(vipmean,width=rep(0.25,8),xlim=c(0,10),ylim=c(0,maxy),ylab="VIP",col=ccolor[1:top],cex.names=1.2,las=2,main="VIP score",cex.axis=1.5,cex=1.5,cex.lab=1.5)
        segments(bar, vipmean - sdbeta, bar, vipmean + sdbeta, lwd=1)
        segments(bar - 0.1, vipmean - sdbeta, bar + 0.1, vipmean - sdbeta, lwd=1)
        segments(bar - 0.1, vipmean + sdbeta, bar + 0.1, vipmean + sdbeta, lwd=1)
        legend("topright",legend = groupnames2, pch = 16, col = unique(color2))
        dev.off()

        pdf((PDF4),width=10,height=10)
        biplot.color(G,H,cex=c(0.8,0.6),col=c(color2),main = paste("PLSDA - ",name,"_",groupnames[i],"_",groupnames[j], sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=marks2,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames2, pch = unique(marks2), col = unique(color2),cex=0.6)

        biplot.color(G,matrixH[1:top,],cex=c(0.8,0.6),col=c(color2),main = paste("PLSDA - ",name,"_",groupnames[i],"_",groupnames[j],sep="")
                     ,xlabMY="PC1",ylabMY="PC2", pch=marks2,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames2, pch = unique(marks2), col = unique(color2),cex=0.6)

        biplot.color_transp(G,H,cex=c(0.8,0.6),col=c(color2),main = paste("PLSDA - ",name,"_",groupnames[i],"_",groupnames[j],sep="")
                            ,xlabMY="PC1",ylabMY="PC2", pch=marks2,lwd=0.00000001,arrow.len=0.05)
        legend("topleft",legend = groupnames2, pch = unique(marks2), col = unique(color2),cex=0.6)

        gr=as.factor(groups2)
        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames2, levels=0.75, plot.points=TRUE, center.cex=0.2, col=unique(color2),
                    pch=unique(marks2),center.pch=3,main=paste("PLSDA - ",name,"_",groupnames[i],"_",groupnames[j],sep=""),xlab="PC1",ylab="PC2")
        legend("bottomright",legend = groupnames2, pch = unique(marks2), col = unique(color2), cex=0.6)

        dataEllipse(G[,1],G[,2],gr,group.labels=groupnames2, levels=0.95, plot.points=TRUE, center.cex=0.2, col=unique(color2),
                    pch=unique(marks2),center.pch=3,main=paste("PLSDA - ",name,"_",groupnames[i],"_",groupnames[j],sep=""),
                    xlab="PC1",ylab="PC2")
        legend("bottomright",legend = groupnames2, pch = unique(marks2), col = unique(color2), cex=0.6)

        dev.off()

        }
    }

    setwd(dirout)
}

