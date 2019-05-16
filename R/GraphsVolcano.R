#' Logratio volcano plots
#'
#' Display logratio volcano plot comparing two specific groups.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param type Parametric ("par" - default) or nonparametric ("nonpar") test must be chosen.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param pair logical. If TRUE then the paired tests are done. Default is FALSE.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @param hy Important points are evaluted on y-axis with this limit. The default is alpha=0.05 with Bonferroni correction.
#' @param hx logical. If TRUE important points are evaluated from both, x-axis and y-axis. Default is FALSE.
#' @param hhx If hx=TRUE, set the absolute value of the limit for most important differences of medians (deafult is 1.5).
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Logratio volcano plot can be used only for comparison of two groups. If there is more groups in data, all possible combinations of pairs are evaluated.
#' @details Logratio volcano plot is a version of standard volcano plot with difference of group medians on x-axis (only for tsf=c("clr","log")) and -log10 of p-value from t-test or Wilcoxon test on y-axis. For "pareto" transformation standard volcano plot is done with log2(fold change) on x-axis.
#' @return Volcano plots.
#' @return Excel sheet with differences of medians and p-values from chosen test.
#' @importFrom stats t.test wilcox.test
#' @importFrom robCompositions cenLR
#' @import openxlsx
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsVolcano(data,name,groupnames)
#' @export
GraphsVolcano=function(data,name,groupnames,type="par",tsf="clr",pair=FALSE,QCs=FALSE,hy=(0.05/length(colnames(data))),hx=FALSE,hhx=1.5){

    ##########################################################################################################################
    count=length(groupnames)
    groupss=matrix(rep(NA,count*2),ncol=2)
    colnames(groupss)=c("min","max")
    rownames(groupss)=groupnames
    groups=NULL
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(data))
        groupss[i,1]=min(Gr)
        groupss[i,2]=max(Gr)
        gr=rep(i,length(Gr))
        groups=c(groups,gr)
    }

####################################################################################################
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
if (QCs==FALSE){
    QC=grep("QC",rownames(dataM))
    if (length(QC)!=0){
        dataM=dataM[-QC,]
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
    ####################################################################################################
    if (type == "par" & (tsf == "clr" | tsf=="log"| tsf=="log10"|tsf=="PQN")) {
    dirout = paste(getwd(),"/",sep = "")
    PDF1=paste(dirout,"Volcano_par_all_",name,".pdf",sep="")
    PDF2=paste(dirout,"Volcano_par_",name,".pdf",sep="")

    wb <- createWorkbook()
    sheet1  <- addWorksheet(wb, sheetName = "T-test")
    sheet2  <- addWorksheet(wb, sheetName = "Diff of med")

    file00=paste(dirout, "Volcano_par_",name,".xlsx",sep="")

    ####################################################################################################
    names<-colnames(dataM)

    if (pair == TRUE) {
        #t-test:

        pvaluew=NULL
        for(k in 1:(ncol(dataM))){
            pvaluew2=NULL
            for(i in 1:(nrow(groupss)-1)){
                A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                for(j in (i+1):nrow(groupss)){
                    B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                    c=t.test(A,B,alternative = "two.sided",paired=TRUE)$p.value
                    pvaluew2=cbind(pvaluew2,c)
                }
            }
            pvaluew=rbind(pvaluew,pvaluew2)
        }

        nazvy=NULL
        for(i in 1:(nrow(groupss)-1)){
            for(j in (i+1):nrow(groupss)){
                col=paste(groupnames[i],"_",groupnames[j],sep="")
                nazvy=c(nazvy,col)
            }
        }

        colnames(pvaluew)=nazvy
        row.names(pvaluew)=names

        writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
            }else{

    #t-test:

    pvaluew=NULL
    for(k in 1:(ncol(dataM))){
        pvaluew2=NULL
        for(i in 1:(nrow(groupss)-1)){
            A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
            for(j in (i+1):nrow(groupss)){
                B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                c=t.test(A,B,alternative = "two.sided")$p.value
                pvaluew2=cbind(pvaluew2,c)
            }
        }
        pvaluew=rbind(pvaluew,pvaluew2)
    }

    nazvy=NULL
    for(i in 1:(nrow(groupss)-1)){
        for(j in (i+1):nrow(groupss)){
            col=paste(groupnames[i],"_",groupnames[j],sep="")
            nazvy=c(nazvy,col)
        }
    }

    colnames(pvaluew)=nazvy
    row.names(pvaluew)=names

    writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
    }
    ####################################################################################################
    #Difference of medians:
    med=matrix(c(rep(0,ncol(dataM)*length(unique(groups)))),nrow=ncol(dataM))
    rownames(med)=colnames(dataM)
    colnames(med)=groupnames

    for(i in 1:ncol(dataM)){
        dataS=dataM[,i]
        groups=groups
        med[i,]=tapply(dataS, groups, median)
        }

    rozdily=NULL
    for(i in 1:(ncol(med)-1)){
        a=med[,i]
        for(j in (i+1):ncol(med)){
            b=med[,j]
            c=a-b
            rozdily=cbind(rozdily,c)
        }
    }

    colnames(rozdily)=nazvy

    writeData(wb,sheet2,rozdily,colNames = TRUE, rowNames = TRUE)
    saveWorkbook(wb, file = file00,overwrite = TRUE)

    if (hx == TRUE) {

    hranicex=hhx
    hranicey=hy

    pdf((PDF1),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
    for(i in 1:ncol(pvaluew)){
        plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
        abline(v=0,col="gray")
        text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
        barevne=which(((abs(rozdily[,i])) > hranicex) & (pvaluew[,i] < hranicey))
        points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
        #text(rozdily[barevne,i],-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
        }
    dev.off()

    pdf((PDF2),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
    for(i in 1:ncol(pvaluew)){
        plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
        abline(v=0,col="gray")
        #text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
        barevne=which(((abs(rozdily[,i])) > hranicex) & (pvaluew[,i] < hranicey))
        points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
        if(length(barevne)!=0){
        text(rozdily[barevne,i],-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
        }
    }
    dev.off()
    }

    if (hx == FALSE) {

        hranicey=hy

        pdf((PDF1),width=10,height=10)
        par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
        for(i in 1:ncol(pvaluew)){
            plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
            abline(v=0,col="gray")
            text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
            barevne=which(pvaluew[,i] < hranicey)
            points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
            #text(rozdily[barevne,i],-log(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
        }
        dev.off()

        pdf((PDF2),width=10,height=10)
        par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
        for(i in 1:ncol(pvaluew)){
            plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
            abline(v=0,col="gray")
            #text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
            barevne=which(pvaluew[,i] < hranicey)
            points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
            if(length(barevne)!=0){
            text(rozdily[barevne,i],-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
            }
        }
        dev.off()
    }

    }
####################################################################################################
####################################################################################################
    if (type == "par" & (tsf == "pareto"|tsf=="none")) {
        dirout = paste(getwd(),"/",sep = "")
        PDF1=paste(dirout,"Volcano_par_all_",name,".pdf",sep="")
        PDF2=paste(dirout,"Volcano_par_",name,".pdf",sep="")

        wb <- createWorkbook()
        sheet1  <- addWorksheet(wb, sheetName = "T-test")
        sheet2  <- addWorksheet(wb, sheetName = "Diff of med")

        file00=paste(dirout, "Volcano_par_",name,".xlsx",sep="")

        ####################################################################################################
        names<-colnames(dataM)

        if (pair == TRUE) {
            #t-test:

            pvaluew=NULL
            for(k in 1:(ncol(dataM))){
                pvaluew2=NULL
                for(i in 1:(nrow(groupss)-1)){
                    A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                    for(j in (i+1):nrow(groupss)){
                        B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                        c=t.test(A,B,alternative = "two.sided",paired=TRUE)$p.value
                        pvaluew2=cbind(pvaluew2,c)
                    }
                }
                pvaluew=rbind(pvaluew,pvaluew2)
            }

            nazvy=NULL
            for(i in 1:(nrow(groupss)-1)){
                for(j in (i+1):nrow(groupss)){
                    col=paste(groupnames[i],"_",groupnames[j],sep="")
                    nazvy=c(nazvy,col)
                }
            }

            colnames(pvaluew)=nazvy
            row.names(pvaluew)=names

            writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
            }else{
        #t-test:

        pvaluew=NULL
        for(k in 1:(ncol(dataM))){
            pvaluew2=NULL
            for(i in 1:(nrow(groupss)-1)){
                A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                for(j in (i+1):nrow(groupss)){
                    B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                    c=t.test(A,B,alternative = "two.sided")$p.value
                    pvaluew2=cbind(pvaluew2,c)
                }
            }
            pvaluew=rbind(pvaluew,pvaluew2)
        }

        nazvy=NULL
        for(i in 1:(nrow(groupss)-1)){
            for(j in (i+1):nrow(groupss)){
                col=paste(groupnames[i],"_",groupnames[j],sep="")
                nazvy=c(nazvy,col)
            }
        }

        colnames(pvaluew)=nazvy
        row.names(pvaluew)=names

        writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
        }
        ####################################################################################################
        #Log2 fold change:
        med=matrix(c(rep(0,ncol(dataM)*length(unique(groups)))),nrow=ncol(dataM))
        rownames(med)=colnames(dataM)
        colnames(med)=groupnames

        for(i in 1:ncol(dataM)){
            dataS=dataM[,i]
            groups=groups
            med[i,]=tapply(dataS, groups, median)
        }

        podily=NULL
        for(i in 1:(ncol(med)-1)){
            a=med[,i]
            for(j in (i+1):ncol(med)){
                b=med[,j]
                c=a/b
                podily=cbind(podily,c)
            }
        }

        colnames(podily)=nazvy

        writeData(wb,sheet2,podily,colNames = TRUE, rowNames = TRUE)
        saveWorkbook(wb, file = file00,overwrite = TRUE)

        if (hx == TRUE) {

            hranicex=hhx
            hranicey=hy

            pdf((PDF1),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(((abs(log2(podily[,i]))) > hranicex) & (pvaluew[,i] < hranicey))
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
            }
            dev.off()

            pdf((PDF2),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                #text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(((abs(log2(podily[,i]))) > hranicex) & (pvaluew[,i] < hranicey))
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                if(length(barevne)!=0){
                text(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],col="red",pos=4,cex=0.7)
                }
            }
            dev.off()
        }

        if (hx == FALSE) {

            hranicey=hy

            pdf((PDF1),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(pvaluew[,i] < hranicey)
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
            }
            dev.off()

            pdf((PDF2),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                #text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(pvaluew[,i] < hranicey)
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                if(length(barevne)!=0){
                text(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],col="red",pos=4,cex=0.7)
                }
            }
            dev.off()
        }

    }
    ####################################################################################################
    ####################################################################################################

    if (type == "nonpar" & (tsf == "clr" | tsf=="log" | tsf=="log10"|tsf=="PQN")) {
        dirout = paste(getwd(),"/",sep = "")
        PDF1=paste(dirout,"Volcano_nonpar_all_",name,".pdf",sep="")
        PDF2=paste(dirout,"Volcano_nonpar_",name,".pdf",sep="")

        wb <- createWorkbook()
        sheet1  <- addWorksheet(wb, sheetName = "Wilcoxon")
        sheet2  <- addWorksheet(wb, sheetName = "Diff of med")

        file00=paste(dirout, "Volcano_nonpar_",name,".xlsx",sep="")

        ####################################################################################################
        names<-colnames(dataM)

        if (pair == TRUE) {
            #Wilcoxon test:

            pvaluew=NULL
            for(k in 1:(ncol(dataM))){
                pvaluew2=NULL
                for(i in 1:(nrow(groupss)-1)){
                    A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                    for(j in (i+1):nrow(groupss)){
                        B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                        c=wilcox.test(A,B,alternative = "two.sided",paired=TRUE)$p.value
                        pvaluew2=cbind(pvaluew2,c)
                    }
                }
                pvaluew=rbind(pvaluew,pvaluew2)
            }

            nazvy=NULL
            for(i in 1:(nrow(groupss)-1)){
                for(j in (i+1):nrow(groupss)){
                    col=paste(groupnames[i],"_",groupnames[j],sep="")
                    nazvy=c(nazvy,col)
                }
            }

            colnames(pvaluew)=nazvy
            row.names(pvaluew)=names

            writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)

        }else{
        #Wilcoxon test:

        pvaluew=NULL
        for(k in 1:(ncol(dataM))){
            pvaluew2=NULL
            for(i in 1:(nrow(groupss)-1)){
                A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                for(j in (i+1):nrow(groupss)){
                    B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                    c=wilcox.test(A,B,alternative = "two.sided")$p.value
                    pvaluew2=cbind(pvaluew2,c)
                }
            }
            pvaluew=rbind(pvaluew,pvaluew2)
        }

        nazvy=NULL
        for(i in 1:(nrow(groupss)-1)){
            for(j in (i+1):nrow(groupss)){
                col=paste(groupnames[i],"_",groupnames[j],sep="")
                nazvy=c(nazvy,col)
            }
        }

        colnames(pvaluew)=nazvy
        row.names(pvaluew)=names

        writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
                }
        ####################################################################################################
        #Difference of medians:
        med=matrix(c(rep(0,ncol(dataM)*length(unique(groups)))),nrow=ncol(dataM))
        rownames(med)=colnames(dataM)
        colnames(med)=groupnames

        for(i in 1:ncol(dataM)){
            dataS=dataM[,i]
            groups=groups
            med[i,]=tapply(dataS, groups, median)
        }


        rozdily=NULL
        for(i in 1:(ncol(med)-1)){
            a=med[,i]
            for(j in (i+1):ncol(med)){
                b=med[,j]
                c=a-b
                rozdily=cbind(rozdily,c)
            }
        }

        colnames(rozdily)=nazvy

        writeData(wb,sheet2,rozdily,colNames = TRUE, rowNames = TRUE)
        saveWorkbook(wb, file = file00,overwrite = TRUE)

        hranicex=hhx
        hranicey=hy

        if (hx == TRUE) {

            hranicex=hx
            hranicey=hy

            pdf((PDF1),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(((abs(rozdily[,i])) > hranicex) & (pvaluew[,i] < hranicey))
                points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                #text(rozdily[barevne,i],-log(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
            }
            dev.off()

            pdf((PDF2),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                #text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(((abs(rozdily[,i])) > hranicex) & (pvaluew[,i] < hranicey))
                points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                if(length(barevne)!=0){
                    text(rozdily[barevne,i],-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
                    }
            }
            dev.off()
        }

        if (hx == FALSE) {

            hranicey=hy

            pdf((PDF1),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(pvaluew[,i] < hranicey)
                points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                #text(rozdily[barevne,i],-log(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
            }
            dev.off()

            pdf((PDF2),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(rozdily[,i],-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="Difference of medians",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                #text(rozdily[,i],-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(pvaluew[,i] < hranicey)
                points(rozdily[barevne,i],-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                if(length(barevne)!=0){
                    text(rozdily[barevne,i],-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
                }
            }
            dev.off()
        }

    }

    ####################################################################################################
    ####################################################################################################
    if (type == "nonpar" & (tsf == "pareto"|tsf=="none")) {
        dirout = paste(getwd(),"/",sep = "")
        PDF1=paste(dirout,"Volcano_nonpar_all_",name,".pdf",sep="")
        PDF2=paste(dirout,"Volcano_nonpar_",name,".pdf",sep="")

        wb <- createWorkbook()
        sheet1  <- addWorksheet(wb, sheetName = "T-test")
        sheet2  <- addWorksheet(wb, sheetName = "Diff of med")

        file00=paste(dirout, "Volcano_nonpar_",name,".xlsx",sep="")

        ####################################################################################################
        names<-colnames(dataM)

        if (pair == TRUE) {
            #Wilcoxon test:

            pvaluew=NULL
            for(k in 1:(ncol(dataM))){
                pvaluew2=NULL
                for(i in 1:(nrow(groupss)-1)){
                    A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                    for(j in (i+1):nrow(groupss)){
                        B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                        c=wilcox.test(A,B,alternative = "two.sided",paired = TRUE)$p.value
                        pvaluew2=cbind(pvaluew2,c)
                    }
                }
                pvaluew=rbind(pvaluew,pvaluew2)
            }

            nazvy=NULL
            for(i in 1:(nrow(groupss)-1)){
                for(j in (i+1):nrow(groupss)){
                    col=paste(groupnames[i],"_",groupnames[j],sep="")
                    nazvy=c(nazvy,col)
                }
            }

            colnames(pvaluew)=nazvy
            row.names(pvaluew)=names

            writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
            }else{
        #Wilcoxon test:

        pvaluew=NULL
        for(k in 1:(ncol(dataM))){
            pvaluew2=NULL
            for(i in 1:(nrow(groupss)-1)){
                A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
                for(j in (i+1):nrow(groupss)){
                    B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                    c=wilcox.test(A,B,alternative = "two.sided")$p.value
                    pvaluew2=cbind(pvaluew2,c)
                }
            }
            pvaluew=rbind(pvaluew,pvaluew2)
        }

        nazvy=NULL
        for(i in 1:(nrow(groupss)-1)){
            for(j in (i+1):nrow(groupss)){
                col=paste(groupnames[i],"_",groupnames[j],sep="")
                nazvy=c(nazvy,col)
            }
        }

        colnames(pvaluew)=nazvy
        row.names(pvaluew)=names

        writeData(wb,sheet1,pvaluew,colNames = TRUE, rowNames = TRUE)
        }
        ####################################################################################################
        #Log2 fold change:
        med=matrix(c(rep(0,ncol(dataM)*length(unique(groups)))),nrow=ncol(dataM))
        rownames(med)=colnames(dataM)
        colnames(med)=groupnames

        for(i in 1:ncol(dataM)){
            dataS=dataM[,i]
            groups=groups
            med[i,]=tapply(dataS, groups, median)
        }

        podily=NULL
        for(i in 1:(ncol(med)-1)){
            a=med[,i]
            for(j in (i+1):ncol(med)){
                b=med[,j]
                c=a/b
                podily=cbind(podily,c)
            }
        }

        colnames(podily)=nazvy

        writeData(wb,sheet2,podily,colNames = TRUE, rowNames = TRUE)
        saveWorkbook(wb, file = file00,overwrite = TRUE)

        if (hx == TRUE) {

            hranicex=hhx
            hranicey=hy

            pdf((PDF1),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(((abs(log2(podily[,i]))) > hranicex) & (pvaluew[,i] < hranicey))
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
            }
            dev.off()

            pdf((PDF2),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                #text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(((abs(log2(podily[,i]))) > hranicex) & (pvaluew[,i] < hranicey))
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                if(length(barevne)!=0){
                    text(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
                }
            }
            dev.off()
        }

        if (hx == FALSE) {

            hranicey=hy

            pdf((PDF1),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(pvaluew[,i] < hranicey)
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
            }
            dev.off()

            pdf((PDF2),width=10,height=10)
            par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))
            for(i in 1:ncol(pvaluew)){
                plot(log2(podily[,i]),-log10(pvaluew[,i]),pch=16,main=colnames(pvaluew)[i],xlab="log2(fold change)",ylab="-log10(p-value)")
                abline(v=0,col="gray")
                #text(log2(podily[,i]),-log10(pvaluew[,i]),labels=rownames(pvaluew),pos=4,cex=0.7)
                barevne=which(pvaluew[,i] < hranicey)
                points(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),col="red",pch=17,cex=1.25)
                if(length(barevne)!=0){
                    text(log2(podily[barevne,i]),-log10(pvaluew[barevne,i]),labels=rownames(pvaluew)[barevne],pos=4,cex=0.7,col="red")
                }
            }
            dev.off()
        }

    }

}

