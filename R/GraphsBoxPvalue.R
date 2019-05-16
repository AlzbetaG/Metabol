#' Boxplots with p-values
#'
#' Display 4 variants of boxplots with p-values from t-test or Wilcoxon test in PDF files for every variable (metabolite) in data table separately.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every ouput.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param type Parametric ("par" - default) or nonparametric ("nonpar") test must be chosen.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to nine different groups can be distinguished in data.
#' @return Boxplot graphs with p-values and stars representing p-values.
#' @return Excel file with medians and differences of these medians of all groups in data.
#' @import graphics
#' @import grDevices
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' GraphsBoxPvalue(data,name,groupnames)
#' @export
GraphsBoxPvalue=function(data,name,groupnames,type="par",tsf="clr",QCs=FALSE){

    ##################################################################################
    dirout=getwd()
    dirout2 = paste(dirout,"/","Boxplots",sep = "")
    dir.create(dirout2, recursive=TRUE, showWarnings = FALSE)
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
####################################################################################################
    count=length(groupnames)
    basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown","gray40","darkslateblue","firebrick","darkcyan","darkmagenta",
                "deeppink1","limegreen")           # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    basemarks=c(15,17,18,8,11,2,0,16,5,6,4,10,3,7,9)
    groupss=matrix(rep(NA,count*2),ncol=2)
    colnames(groupss)=c("min","max")
    rownames(groupss)=groupnames
    groups=NULL
    marks=NULL
    color=NULL
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
        groups=groups[-QC]
        color=color[-QC]
        marks=marks[-QC]
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
##################################################################################

# Boxplots

names<-colnames(dataM)

#names<-abbreviate(colnames(dataM),minlength = 20,dot = FALSE, strict = FALSE,method = c("left.kept", "both.sides"))
##################################################################################
if (type == "par"){
pvalue=NULL
for(k in 1:(ncol(dataM))){
    pvalue2=NULL
    for(i in 1:(nrow(groupss)-1)){
        A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
        for(j in (i+1):nrow(groupss)){
            B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
            c=t.test(A,B,alternative = "two.sided",paired =FALSE)$p.value
            pvalue2=cbind(pvalue2,c)
        }
    }
    pvalue=rbind(pvalue,pvalue2)
}

tip="par"
}
##################################################################################

if (type == "nonpar"){
    pvalue=NULL
    for(k in 1:(ncol(dataM))){
        pvalue2=NULL
        for(i in 1:(nrow(groupss)-1)){
            A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
            for(j in (i+1):nrow(groupss)){
                B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                c=wilcox.test(A,B,alternative = "two.sided",paired =FALSE)$p.value
                pvalue2=cbind(pvalue2,c)
            }
        }
        pvalue=rbind(pvalue,pvalue2)
    }

tip="nonpar"
}

##################################################################################
if (count == 2){
    PDF1=paste("Box_pval_",tip,"_points_",name,".pdf",sep="")
    PDF2=paste("Box_pval_",tip,"_notch_points_",name,".pdf",sep="")
    PDF3=paste("Box_pval_",tip,"_names_",name,".pdf",sep="")
    PDF4=paste("Box_pval_",tip,"_notch_names_",name,".pdf",sep="")
}else{
    PDF1=paste("Box_pval_",tip,"_stars_",name,".pdf",sep="")
    PDF2=paste("Box_pval_",tip,"_notch_stars_",name,".pdf",sep="")
    PDF3=paste("Box_pval_",tip,"_numbers_",name,".pdf",sep="")
    PDF4=paste("Box_pval_",tip,"_notch_numbers_",name,".pdf",sep="")
}

##################################################################################
if (count == 2){

pdf((PDF1),width=10,height=10)
par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

for(i in 1:ncol(dataM)){
    boxplot(dataM[,i] ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25)
    stripchart(dataM[,i] ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
    mtext("p=",3,at=2.2,cex=1.25)

    if (pvalue[i,1]>0.0001){
        a=round(pvalue[i,1],4)
    } else {
        a=signif(pvalue[i,1],3)
    }
    mtext(a,3,at=2.4,cex=1.25)
}
dev.off()


pdf((PDF2),width=10,height=10)
par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

for(i in 1:ncol(dataM)){
    boxplot(dataM[,i] ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,notch=TRUE)
    stripchart(dataM[,i] ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
    mtext("p=",3,at=2.2,cex=1.25)

    if (pvalue[i,1]>0.0001){
        a=round(pvalue[i,1],4)
    } else {
        a=signif(pvalue[i,1],3)
    }
    mtext(a,3,at=2.4,cex=1.25)
}
dev.off()


pdf((PDF3),width=10,height=10)
par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

labels=rownames(dataM)

for(i in 1:ncol(dataM)){
    box=boxplot(dataM[,i] ~ groups,outpch=NA,names=groupnames,main=names[i],cex.axis=1.25)
    par(new=TRUE)
    text(groups,dataM[,i],label=labels,col="red")
    mtext("p=",3,at=2.2,cex=1.25)

    if (pvalue[i,1]>0.0001){
        a=round(pvalue[i,1],4)
    } else {
        a=signif(pvalue[i,1],3)
    }
    mtext(a,3,at=2.4,cex=1.25)
    }
dev.off()

pdf((PDF4),width=10,height=10)
par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

labels=rownames(dataM)

for(i in 1:ncol(dataM)){
    box=boxplot(dataM[,i] ~ groups,outpch=NA,names=groupnames,main=names[i],cex.axis=1.25,notch=TRUE)
    par(new=TRUE)
    text(groups,dataM[,i],label=labels,col="red")
    mtext("p=",3,at=2.2,cex=1.25)

    if (pvalue[i,1]>0.0001){
        a=round(pvalue[i,1],4)
    } else {
        a=signif(pvalue[i,1],3)
    }
    mtext(a,3,at=2.4,cex=1.25)
}
dev.off()
}
##################################################################################
if (count == 3){

    pdf((PDF1),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=3.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=3.5,cex=1.25)
        ma2=max(dataM[,i])+konst/3
        ma3=max(dataM[,i])+konst*2/3
        ma4=max(dataM[,i])+konst/3
        ma22=max(dataM[,i])+konst/2
        ma32=max(dataM[,i])+konst
        ma42=max(dataM[,i])+konst/2


        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.52,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.02,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(2.1,ma4,2.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(2.1,ma4,2.9,ma4)
            text(2.52,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(2.1,ma4,2.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }


    }
    dev.off()

    pdf((PDF2),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=3.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=3.5,cex=1.25)
        ma2=max(dataM[,i])+konst/3
        ma3=max(dataM[,i])+konst*2/3
        ma4=max(dataM[,i])+konst/3
        ma22=max(dataM[,i])+konst/2
        ma32=max(dataM[,i])+konst
        ma42=max(dataM[,i])+konst/2


        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.52,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.02,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(2.1,ma4,2.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(2.1,ma4,2.9,ma4)
            text(2.52,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(2.1,ma4,2.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF3),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=3.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=3.5,cex=1.25)
        ma2=max(dataM[,i])+konst/3
        ma3=max(dataM[,i])+konst*2/3
        ma4=max(dataM[,i])+konst/3
        ma22=max(dataM[,i])+konst/2
        ma32=max(dataM[,i])+konst
        ma42=max(dataM[,i])+konst/2

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(2.1,ma4,2.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(2.1,ma4,2.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=3.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=3.5,cex=1.25)
        ma2=max(dataM[,i])+konst/3
        ma3=max(dataM[,i])+konst*2/3
        ma4=max(dataM[,i])+konst/3
        ma22=max(dataM[,i])+konst/2
        ma32=max(dataM[,i])+konst
        ma42=max(dataM[,i])+konst/2

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(2.1,ma4,2.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(2.1,ma4,2.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}
    }
    dev.off()

}
##################################################################################
if (count == 4){

    pdf((PDF1),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/5
        ma3=max(dataM[,i])+konst*2/5
        ma4=max(dataM[,i])+konst*3/5
        ma5=max(dataM[,i])+konst/5
        ma6=max(dataM[,i])+konst*4/5
        ma7=max(dataM[,i])+konst/5
        ma22=max(dataM[,i])+konst/4
        ma32=max(dataM[,i])+konst/2
        ma42=max(dataM[,i])+konst*3/4
        ma52=max(dataM[,i])+konst/4
        ma62=max(dataM[,i])+konst
        ma72=max(dataM[,i])+konst/4

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.53,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.03,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.53,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(2.1,ma5,2.9,ma5)
            text(2.5,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(2.1,ma5,2.9,ma5)
            text(2.53,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(2.1,ma5,2.9,ma5)
            text(2.5,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(2.1,ma6,3.9,ma6)
            text(3,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(2.1,ma6,3.9,ma6)
            text(3.03,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(2.1,ma6,3.9,ma6)
            text(3,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(3.1,ma7,3.9,ma7)
            text(3.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(3.1,ma7,3.9,ma7)
            text(3.53,ma72,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(3.1,ma7,3.9,ma7)
            text(3.5,ma72,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF2),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/5
        ma3=max(dataM[,i])+konst*2/5
        ma4=max(dataM[,i])+konst*3/5
        ma5=max(dataM[,i])+konst/5
        ma6=max(dataM[,i])+konst*4/5
        ma7=max(dataM[,i])+konst/5
        ma22=max(dataM[,i])+konst/4
        ma32=max(dataM[,i])+konst/2
        ma42=max(dataM[,i])+konst*3/4
        ma52=max(dataM[,i])+konst/4
        ma62=max(dataM[,i])+konst
        ma72=max(dataM[,i])+konst/4

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.53,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.03,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.53,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(2.1,ma5,2.9,ma5)
            text(2.5,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(2.1,ma5,2.9,ma5)
            text(2.53,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(2.1,ma5,2.9,ma5)
            text(2.5,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(2.1,ma6,3.9,ma6)
            text(3,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(2.1,ma6,3.9,ma6)
            text(3.03,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(2.1,ma6,3.9,ma6)
            text(3,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(3.1,ma7,3.9,ma7)
            text(3.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(3.1,ma7,3.9,ma7)
            text(3.53,ma72,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(3.1,ma7,3.9,ma7)
            text(3.5,ma72,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF3),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=4.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=4.5,cex=1.25)
        ma2=max(dataM[,i])+konst/5
        ma3=max(dataM[,i])+konst*2/5
        ma4=max(dataM[,i])+konst*3/5
        ma5=max(dataM[,i])+konst/5
        ma6=max(dataM[,i])+konst*4/5
        ma7=max(dataM[,i])+konst/5
        ma22=max(dataM[,i])+konst/4
        ma32=max(dataM[,i])+konst/2
        ma42=max(dataM[,i])+konst*3/4
        ma52=max(dataM[,i])+konst/4
        ma62=max(dataM[,i])+konst
        ma72=max(dataM[,i])+konst/4

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(2.1,ma5,2.9,ma5)
                text(2.25,ma52,"p=",cex=0.75)
                text(2.5,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(2.1,ma5,2.9,ma5)
                text(2.25,ma52,"p=",cex=0.75)
                text(2.5,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(2.1,ma6,3.9,ma6)
                text(2.75,ma62,"p=",cex=0.75)
                text(3,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(2.1,ma6,3.9,ma6)
                text(2.75,ma62,"p=",cex=0.75)
                text(3,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(3.1,ma7,3.9,ma7)
                text(3.25,ma72,"p=",cex=0.75)
                text(3.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(3.1,ma7,3.9,ma7)
                text(3.25,ma72,"p=",cex=0.75)
                text(3.5,ma72,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=10,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=4.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=4.5,cex=1.25)
        ma2=max(dataM[,i])+konst/5
        ma3=max(dataM[,i])+konst*2/5
        ma4=max(dataM[,i])+konst*3/5
        ma5=max(dataM[,i])+konst/5
        ma6=max(dataM[,i])+konst*4/5
        ma7=max(dataM[,i])+konst/5
        ma22=max(dataM[,i])+konst/4
        ma32=max(dataM[,i])+konst/2
        ma42=max(dataM[,i])+konst*3/4
        ma52=max(dataM[,i])+konst/4
        ma62=max(dataM[,i])+konst
        ma72=max(dataM[,i])+konst/4

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(2.1,ma5,2.9,ma5)
                text(2.25,ma52,"p=",cex=0.75)
                text(2.5,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(2.1,ma5,2.9,ma5)
                text(2.25,ma52,"p=",cex=0.75)
                text(2.5,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(2.1,ma6,3.9,ma6)
                text(2.75,ma62,"p=",cex=0.75)
                text(3,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(2.1,ma6,3.9,ma6)
                text(2.75,ma62,"p=",cex=0.75)
                text(3,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(3.1,ma7,3.9,ma7)
                text(3.25,ma72,"p=",cex=0.75)
                text(3.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(3.1,ma7,3.9,ma7)
                text(3.25,ma72,"p=",cex=0.75)
                text(3.5,ma72,a,cex=0.75)
            }}
    }
    dev.off()
}
##################################################################################
if (count == 5){

    pdf((PDF1),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/7
        ma3=max(dataM[,i])+konst*2/7
        ma4=max(dataM[,i])+konst*3/7
        ma5=max(dataM[,i])+konst*4/7
        ma6=max(dataM[,i])+konst/7
        ma7=max(dataM[,i])+konst*5/7
        ma8=max(dataM[,i])+konst*6/7
        ma9=max(dataM[,i])+konst/7
        ma10=max(dataM[,i])+konst*2/7
        ma11=max(dataM[,i])+konst/7
        ma22=max(dataM[,i])+konst/6
        ma32=max(dataM[,i])+konst*2/6
        ma42=max(dataM[,i])+konst*3/6
        ma52=max(dataM[,i])+konst*4/6
        ma62=max(dataM[,i])+konst/6
        ma72=max(dataM[,i])+konst*4.5/6
        ma82=max(dataM[,i])+konst*5.5/6
        ma92=max(dataM[,i])+konst/6
        ma102=max(dataM[,i])+konst*2/6
        ma112=max(dataM[,i])+konst/6

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.53,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.03,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.53,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.03,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(2.1,ma6,2.9,ma6)
            text(2.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(2.1,ma6,2.9,ma6)
            text(2.53,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(2.1,ma6,2.9,ma6)
            text(2.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.03,ma72,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.53,ma82,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(3.1,ma9,3.9,ma9)
            text(3.5,ma92,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(3.1,ma9,3.9,ma9)
            text(3.53,ma92,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(3.1,ma9,3.9,ma9)
            text(3.5,ma92,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(3.1,ma10,4.9,ma10)
            text(4,ma102,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(3.1,ma10,4.9,ma10)
            text(4.03,ma102,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(3.1,ma10,4.9,ma10)
            text(4,ma102,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(4.1,ma11,4.9,ma11)
            text(4.5,ma112,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(4.1,ma11,4.9,ma11)
            text(4.53,ma112,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(4.1,ma11,4.9,ma11)
            text(4.5,ma112,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF2),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/7
        ma3=max(dataM[,i])+konst*2/7
        ma4=max(dataM[,i])+konst*3/7
        ma5=max(dataM[,i])+konst*4/7
        ma6=max(dataM[,i])+konst/7
        ma7=max(dataM[,i])+konst*5/7
        ma8=max(dataM[,i])+konst*6/7
        ma9=max(dataM[,i])+konst/7
        ma10=max(dataM[,i])+konst*2/7
        ma11=max(dataM[,i])+konst/7
        ma22=max(dataM[,i])+konst/6
        ma32=max(dataM[,i])+konst*2/6
        ma42=max(dataM[,i])+konst*3/6
        ma52=max(dataM[,i])+konst*4/6
        ma62=max(dataM[,i])+konst/6
        ma72=max(dataM[,i])+konst*4.5/6
        ma82=max(dataM[,i])+konst*5.5/6
        ma92=max(dataM[,i])+konst/6
        ma102=max(dataM[,i])+konst*2/6
        ma112=max(dataM[,i])+konst/6

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.53,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.03,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.53,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.03,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(2.1,ma6,2.9,ma6)
            text(2.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(2.1,ma6,2.9,ma6)
            text(2.53,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(2.1,ma6,2.9,ma6)
            text(2.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.03,ma72,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.53,ma82,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(3.1,ma9,3.9,ma9)
            text(3.5,ma92,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(3.1,ma9,3.9,ma9)
            text(3.53,ma92,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(3.1,ma9,3.9,ma9)
            text(3.5,ma92,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(3.1,ma10,4.9,ma10)
            text(4,ma102,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(3.1,ma10,4.9,ma10)
            text(4.03,ma102,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(3.1,ma10,4.9,ma10)
            text(4,ma102,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(4.1,ma11,4.9,ma11)
            text(4.5,ma112,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(4.1,ma11,4.9,ma11)
            text(4.53,ma112,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(4.1,ma11,4.9,ma11)
            text(4.5,ma112,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF3),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=4.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=4.5,cex=1.25)
        ma2=max(dataM[,i])+konst/7
        ma3=max(dataM[,i])+konst*2/7
        ma4=max(dataM[,i])+konst*3/7
        ma5=max(dataM[,i])+konst*4/7
        ma6=max(dataM[,i])+konst/7
        ma7=max(dataM[,i])+konst*5/7
        ma8=max(dataM[,i])+konst*6/7
        ma9=max(dataM[,i])+konst/7
        ma10=max(dataM[,i])+konst*2/7
        ma11=max(dataM[,i])+konst/7
        ma22=max(dataM[,i])+konst/6
        ma32=max(dataM[,i])+konst*2/6
        ma42=max(dataM[,i])+konst*3/6
        ma52=max(dataM[,i])+konst*4/6
        ma62=max(dataM[,i])+konst/6
        ma72=max(dataM[,i])+konst*4.5/6
        ma82=max(dataM[,i])+konst*5.5/6
        ma92=max(dataM[,i])+konst/6
        ma102=max(dataM[,i])+konst*2/6
        ma112=max(dataM[,i])+konst/6

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(2.1,ma6,2.9,ma6)
                text(2.25,ma62,"p=",cex=0.75)
                text(2.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(2.1,ma6,2.9,ma6)
                text(2.25,ma62,"p=",cex=0.75)
                text(2.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(3.1,ma9,3.9,ma9)
                text(3.25,ma92,"p=",cex=0.75)
                text(3.5,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(3.1,ma9,3.9,ma9)
                text(3.25,ma92,"p=",cex=0.75)
                text(3.5,ma92,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(3.1,ma10,4.9,ma10)
                text(3.75,ma102,"p=",cex=0.75)
                text(4,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(3.1,ma10,4.9,ma10)
                text(3.75,ma102,"p=",cex=0.75)
                text(4,ma102,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(4.1,ma11,4.9,ma11)
                text(4.25,ma112,"p=",cex=0.75)
                text(4.5,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(4.1,ma11,4.9,ma11)
                text(4.25,ma112,"p=",cex=0.75)
                text(4.5,ma112,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/3
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=4.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=4.5,cex=1.25)
        ma2=max(dataM[,i])+konst/7
        ma3=max(dataM[,i])+konst*2/7
        ma4=max(dataM[,i])+konst*3/7
        ma5=max(dataM[,i])+konst*4/7
        ma6=max(dataM[,i])+konst/7
        ma7=max(dataM[,i])+konst*5/7
        ma8=max(dataM[,i])+konst*6/7
        ma9=max(dataM[,i])+konst/7
        ma10=max(dataM[,i])+konst*2/7
        ma11=max(dataM[,i])+konst/7
        ma22=max(dataM[,i])+konst/6
        ma32=max(dataM[,i])+konst*2/6
        ma42=max(dataM[,i])+konst*3/6
        ma52=max(dataM[,i])+konst*4/6
        ma62=max(dataM[,i])+konst/6
        ma72=max(dataM[,i])+konst*4.5/6
        ma82=max(dataM[,i])+konst*5.5/6
        ma92=max(dataM[,i])+konst/6
        ma102=max(dataM[,i])+konst*2/6
        ma112=max(dataM[,i])+konst/6

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(2.1,ma6,2.9,ma6)
                text(2.25,ma62,"p=",cex=0.75)
                text(2.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(2.1,ma6,2.9,ma6)
                text(2.25,ma62,"p=",cex=0.75)
                text(2.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(3.1,ma9,3.9,ma9)
                text(3.25,ma92,"p=",cex=0.75)
                text(3.5,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(3.1,ma9,3.9,ma9)
                text(3.25,ma92,"p=",cex=0.75)
                text(3.5,ma92,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(3.1,ma10,4.9,ma10)
                text(3.75,ma102,"p=",cex=0.75)
                text(4,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(3.1,ma10,4.9,ma10)
                text(3.75,ma102,"p=",cex=0.75)
                text(4,ma102,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(4.1,ma11,4.9,ma11)
                text(4.25,ma112,"p=",cex=0.75)
                text(4.5,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(4.1,ma11,4.9,ma11)
                text(4.25,ma112,"p=",cex=0.75)
                text(4.5,ma112,a,cex=0.75)
            }}
    }
dev.off()

}
##################################################################################
if (count == 6){

    pdf((PDF1),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=5.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=5.5,cex=1.25)
        ma2=max(dataM[,i])+konst/10
        ma3=max(dataM[,i])+konst*2/10
        ma4=max(dataM[,i])+konst*3/10
        ma5=max(dataM[,i])+konst*4/10
        ma6=max(dataM[,i])+konst*5/10
        ma7=max(dataM[,i])+konst*6/10
        ma8=max(dataM[,i])+konst*7/10
        ma9=max(dataM[,i])+konst*8/10
        ma10=max(dataM[,i])+konst*9/10
        ma22=max(dataM[,i])+konst*1.25/9
        ma32=max(dataM[,i])+konst*2/9
        ma42=max(dataM[,i])+konst*3/9
        ma52=max(dataM[,i])+konst*4/9
        ma62=max(dataM[,i])+konst*5/9
        ma72=max(dataM[,i])+konst*6/9
        ma82=max(dataM[,i])+konst*7/9
        ma92=max(dataM[,i])+konst*7.5/9
        ma102=max(dataM[,i])+konst*8.5/9

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }


        if (pvalue[i,8] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }


    }
    dev.off()

    pdf((PDF2),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=5.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=5.5,cex=1.25)
        ma2=max(dataM[,i])+konst/10
        ma3=max(dataM[,i])+konst*2/10
        ma4=max(dataM[,i])+konst*3/10
        ma5=max(dataM[,i])+konst*4/10
        ma6=max(dataM[,i])+konst*5/10
        ma7=max(dataM[,i])+konst*6/10
        ma8=max(dataM[,i])+konst*7/10
        ma9=max(dataM[,i])+konst*8/10
        ma10=max(dataM[,i])+konst*9/10
        ma22=max(dataM[,i])+konst*1.25/9
        ma32=max(dataM[,i])+konst*2/9
        ma42=max(dataM[,i])+konst*3/9
        ma52=max(dataM[,i])+konst*4/9
        ma62=max(dataM[,i])+konst*5/9
        ma72=max(dataM[,i])+konst*6/9
        ma82=max(dataM[,i])+konst*7/9
        ma92=max(dataM[,i])+konst*7.5/9
        ma102=max(dataM[,i])+konst*8.5/9

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }


        if (pvalue[i,8] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF3),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/10
        ma3=max(dataM[,i])+konst*2/10
        ma4=max(dataM[,i])+konst*3/10
        ma5=max(dataM[,i])+konst*4/10
        ma6=max(dataM[,i])+konst*5/10
        ma7=max(dataM[,i])+konst*6/10
        ma8=max(dataM[,i])+konst*7/10
        ma9=max(dataM[,i])+konst*8.25/10
        ma10=max(dataM[,i])+konst*9/10
        ma22=max(dataM[,i])+konst*1.25/9
        ma32=max(dataM[,i])+konst*2/9
        ma42=max(dataM[,i])+konst*3/9
        ma52=max(dataM[,i])+konst*4/9
        ma62=max(dataM[,i])+konst*5/9
        ma72=max(dataM[,i])+konst*6/9
        ma82=max(dataM[,i])+konst*7/9
        ma92=max(dataM[,i])+konst*7.75/9
        ma102=max(dataM[,i])+konst*8.5/9

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/10
        ma3=max(dataM[,i])+konst*2/10
        ma4=max(dataM[,i])+konst*3/10
        ma5=max(dataM[,i])+konst*4/10
        ma6=max(dataM[,i])+konst*5/10
        ma7=max(dataM[,i])+konst*6/10
        ma8=max(dataM[,i])+konst*7/10
        ma9=max(dataM[,i])+konst*8.25/10
        ma10=max(dataM[,i])+konst*9/10
        ma22=max(dataM[,i])+konst*1.25/9
        ma32=max(dataM[,i])+konst*2/9
        ma42=max(dataM[,i])+konst*3/9
        ma52=max(dataM[,i])+konst*4/9
        ma62=max(dataM[,i])+konst*5/9
        ma72=max(dataM[,i])+konst*6/9
        ma82=max(dataM[,i])+konst*7/9
        ma92=max(dataM[,i])+konst*7.75/9
        ma102=max(dataM[,i])+konst*8.5/9

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

}
##################################################################################
if (count == 7){

    pdf((PDF1),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=5.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=5.5,cex=1.25)
        ma2=max(dataM[,i])+konst/13
        ma3=max(dataM[,i])+konst*2/13
        ma4=max(dataM[,i])+konst*3/13
        ma5=max(dataM[,i])+konst*4/13
        ma6=max(dataM[,i])+konst*5/13
        ma7=max(dataM[,i])+konst*6/13
        ma8=max(dataM[,i])+konst*7/13
        ma9=max(dataM[,i])+konst*8/13
        ma10=max(dataM[,i])+konst*9/13
        ma11=max(dataM[,i])+konst*10/13
        ma12=max(dataM[,i])+konst*11/13
        ma13=max(dataM[,i])+konst*12/13
        ma22=max(dataM[,i])+konst*1.25/12
        ma32=max(dataM[,i])+konst*2/12
        ma42=max(dataM[,i])+konst*3/12
        ma52=max(dataM[,i])+konst*4/12
        ma62=max(dataM[,i])+konst*5/12
        ma72=max(dataM[,i])+konst*6/12
        ma82=max(dataM[,i])+konst*7/12
        ma92=max(dataM[,i])+konst*7.5/12
        ma102=max(dataM[,i])+konst*8.5/12
        ma112=max(dataM[,i])+konst*9.5/12
        ma122=max(dataM[,i])+konst*10.5/12
        ma132=max(dataM[,i])+konst*11.5/12

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(1.1,ma11,6.9,ma11)
            text(4.04,ma112,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(2.1,ma12,6.9,ma12)
            text(4.54,ma122,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(3.1,ma13,6.9,ma13)
            text(5.04,ma132,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"*",cex=1.5)
        }

        if (pvalue[i,16] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,16] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,16] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,17] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,17] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,17] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,18] < 0.001) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,18] < 0.01) {
            segments(4.1,ma7,6.9,ma7)
            text(5.54,ma72,"**",cex=1.5)
        }

        if (pvalue[i,18] < 0.05) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"*",cex=1.5)
        }

        if (pvalue[i,19] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,19] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,19] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,20] < 0.001) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"***",cex=1.5)
        }

        if (pvalue[i,20] < 0.01) {
            segments(5.1,ma3,6.9,ma3)
            text(6.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,20] < 0.05) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"*",cex=1.5)
        }

        if (pvalue[i,21] < 0.001) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,21] < 0.01) {
            segments(6.1,ma2,6.9,ma2)
            text(6.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,21] < 0.05) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF2),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=5.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=5.5,cex=1.25)
        ma2=max(dataM[,i])+konst/13
        ma3=max(dataM[,i])+konst*2/13
        ma4=max(dataM[,i])+konst*3/13
        ma5=max(dataM[,i])+konst*4/13
        ma6=max(dataM[,i])+konst*5/13
        ma7=max(dataM[,i])+konst*6/13
        ma8=max(dataM[,i])+konst*7/13
        ma9=max(dataM[,i])+konst*8/13
        ma10=max(dataM[,i])+konst*9/13
        ma11=max(dataM[,i])+konst*10/13
        ma12=max(dataM[,i])+konst*11/13
        ma13=max(dataM[,i])+konst*12/13
        ma22=max(dataM[,i])+konst*1.25/12
        ma32=max(dataM[,i])+konst*2/12
        ma42=max(dataM[,i])+konst*3/12
        ma52=max(dataM[,i])+konst*4/12
        ma62=max(dataM[,i])+konst*5/12
        ma72=max(dataM[,i])+konst*6/12
        ma82=max(dataM[,i])+konst*7/12
        ma92=max(dataM[,i])+konst*7.5/12
        ma102=max(dataM[,i])+konst*8.5/12
        ma112=max(dataM[,i])+konst*9.5/12
        ma122=max(dataM[,i])+konst*10.5/12
        ma132=max(dataM[,i])+konst*11.5/12

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(1.1,ma11,6.9,ma11)
            text(4.04,ma112,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(2.1,ma12,6.9,ma12)
            text(4.54,ma122,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(3.1,ma13,6.9,ma13)
            text(5.04,ma132,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"*",cex=1.5)
        }

        if (pvalue[i,16] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,16] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,16] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,17] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,17] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,17] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,18] < 0.001) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,18] < 0.01) {
            segments(4.1,ma7,6.9,ma7)
            text(5.54,ma72,"**",cex=1.5)
        }

        if (pvalue[i,18] < 0.05) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"*",cex=1.5)
        }

        if (pvalue[i,19] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,19] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,19] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,20] < 0.001) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"***",cex=1.5)
        }

        if (pvalue[i,20] < 0.01) {
            segments(5.1,ma3,6.9,ma3)
            text(6.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,20] < 0.05) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"*",cex=1.5)
        }

        if (pvalue[i,21] < 0.001) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,21] < 0.01) {
            segments(6.1,ma2,6.9,ma2)
            text(6.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,21] < 0.05) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF3),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/13
        ma3=max(dataM[,i])+konst*2/13
        ma4=max(dataM[,i])+konst*3/13
        ma5=max(dataM[,i])+konst*4/13
        ma6=max(dataM[,i])+konst*5/13
        ma7=max(dataM[,i])+konst*6/13
        ma8=max(dataM[,i])+konst*7/13
        ma9=max(dataM[,i])+konst*8/13
        ma10=max(dataM[,i])+konst*9/13
        ma11=max(dataM[,i])+konst*10/13
        ma12=max(dataM[,i])+konst*11/13
        ma13=max(dataM[,i])+konst*12/13
        ma22=max(dataM[,i])+konst*1.25/12
        ma32=max(dataM[,i])+konst*2/12
        ma42=max(dataM[,i])+konst*3/12
        ma52=max(dataM[,i])+konst*4/12
        ma62=max(dataM[,i])+konst*5/12
        ma72=max(dataM[,i])+konst*6/12
        ma82=max(dataM[,i])+konst*7/12
        ma92=max(dataM[,i])+konst*7.5/12
        ma102=max(dataM[,i])+konst*8.5/12
        ma112=max(dataM[,i])+konst*9.5/12
        ma122=max(dataM[,i])+konst*10.5/12
        ma132=max(dataM[,i])+konst*11.5/12

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }}

        if (pvalue[i,16] < 0.05) {
            if (pvalue[i,16]>0.0001){
                a=round(pvalue[i,16],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,16],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,17] < 0.05) {
            if (pvalue[i,17]>0.0001){
                a=round(pvalue[i,17],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,17],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,18] < 0.05) {
            if (pvalue[i,18]>0.0001){
                a=round(pvalue[i,18],4)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,18],3)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }}

        if (pvalue[i,19] < 0.05) {
            if (pvalue[i,19]>0.0001){
                a=round(pvalue[i,19],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,19],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,20] < 0.05) {
            if (pvalue[i,20]>0.0001){
                a=round(pvalue[i,20],4)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,20],3)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }}

        if (pvalue[i,21] < 0.05) {
            if (pvalue[i,21]>0.0001){
                a=round(pvalue[i,21],4)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,21],3)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=15,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)/2
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/13
        ma3=max(dataM[,i])+konst*2/13
        ma4=max(dataM[,i])+konst*3/13
        ma5=max(dataM[,i])+konst*4/13
        ma6=max(dataM[,i])+konst*5/13
        ma7=max(dataM[,i])+konst*6/13
        ma8=max(dataM[,i])+konst*7/13
        ma9=max(dataM[,i])+konst*8/13
        ma10=max(dataM[,i])+konst*9/13
        ma11=max(dataM[,i])+konst*10/13
        ma12=max(dataM[,i])+konst*11/13
        ma13=max(dataM[,i])+konst*12/13
        ma22=max(dataM[,i])+konst*1.25/12
        ma32=max(dataM[,i])+konst*2/12
        ma42=max(dataM[,i])+konst*3/12
        ma52=max(dataM[,i])+konst*4/12
        ma62=max(dataM[,i])+konst*5/12
        ma72=max(dataM[,i])+konst*6/12
        ma82=max(dataM[,i])+konst*7/12
        ma92=max(dataM[,i])+konst*7.5/12
        ma102=max(dataM[,i])+konst*8.5/12
        ma112=max(dataM[,i])+konst*9.5/12
        ma122=max(dataM[,i])+konst*10.5/12
        ma132=max(dataM[,i])+konst*11.5/12

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }}

        if (pvalue[i,16] < 0.05) {
            if (pvalue[i,16]>0.0001){
                a=round(pvalue[i,16],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,16],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,17] < 0.05) {
            if (pvalue[i,17]>0.0001){
                a=round(pvalue[i,17],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,17],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,18] < 0.05) {
            if (pvalue[i,18]>0.0001){
                a=round(pvalue[i,18],4)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,18],3)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }}

        if (pvalue[i,19] < 0.05) {
            if (pvalue[i,19]>0.0001){
                a=round(pvalue[i,19],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,19],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,20] < 0.05) {
            if (pvalue[i,20]>0.0001){
                a=round(pvalue[i,20],4)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,20],3)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }}

        if (pvalue[i,21] < 0.05) {
            if (pvalue[i,21]>0.0001){
                a=round(pvalue[i,21],4)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,21],3)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

}
##################################################################################
if (count == 8){

    pdf((PDF1),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=5.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=5.5,cex=1.25)
        ma2=max(dataM[,i])+konst/17
        ma3=max(dataM[,i])+konst*2/17
        ma4=max(dataM[,i])+konst*3/17
        ma5=max(dataM[,i])+konst*4/17
        ma6=max(dataM[,i])+konst*5/17
        ma7=max(dataM[,i])+konst*6/17
        ma8=max(dataM[,i])+konst*7/17
        ma9=max(dataM[,i])+konst*8/17
        ma10=max(dataM[,i])+konst*9/17
        ma11=max(dataM[,i])+konst*10/17
        ma12=max(dataM[,i])+konst*11/17
        ma13=max(dataM[,i])+konst*12/17
        ma14=max(dataM[,i])+konst*13/17
        ma15=max(dataM[,i])+konst*14/17
        ma16=max(dataM[,i])+konst*15/17
        ma17=max(dataM[,i])+konst*16/17
        ma22=max(dataM[,i])+konst*1.25/16
        ma32=max(dataM[,i])+konst*2/16
        ma42=max(dataM[,i])+konst*3/16
        ma52=max(dataM[,i])+konst*4/16
        ma62=max(dataM[,i])+konst*5/16
        ma72=max(dataM[,i])+konst*6/16
        ma82=max(dataM[,i])+konst*7/16
        ma92=max(dataM[,i])+konst*7.5/16
        ma102=max(dataM[,i])+konst*8.5/16
        ma112=max(dataM[,i])+konst*9.5/16
        ma122=max(dataM[,i])+konst*10.5/16
        ma132=max(dataM[,i])+konst*11.5/16
        ma142=max(dataM[,i])+konst*12.5/16
        ma152=max(dataM[,i])+konst*13.5/16
        ma162=max(dataM[,i])+konst*14.5/16
        ma172=max(dataM[,i])+konst*15.5/16

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(1.1,ma11,6.9,ma11)
            text(4.04,ma112,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(1.1,ma14,7.9,ma14)
            text(4.54,ma142,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(2.1,ma12,6.9,ma12)
            text(4.54,ma122,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(2.1,ma15,7.9,ma15)
            text(5.04,ma152,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,16] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,16] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,16] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,17] < 0.001) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"***",cex=1.5)
        }

        if (pvalue[i,17] < 0.01) {
            segments(3.1,ma13,6.9,ma13)
            text(5.04,ma132,"**",cex=1.5)
        }

        if (pvalue[i,17] < 0.05) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"*",cex=1.5)
        }

        if (pvalue[i,18] < 0.001) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"***",cex=1.5)
        }

        if (pvalue[i,18] < 0.01) {
            segments(3.1,ma16,7.9,ma16)
            text(5.54,ma162,"**",cex=1.5)
        }

        if (pvalue[i,18] < 0.05) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"*",cex=1.5)
        }

        if (pvalue[i,19] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,19] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,19] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,20] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,20] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,20] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,21] < 0.001) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,21] < 0.01) {
            segments(4.1,ma7,6.9,ma7)
            text(5.54,ma72,"**",cex=1.5)
        }

        if (pvalue[i,21] < 0.05) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"*",cex=1.5)
        }

        if (pvalue[i,22] < 0.001) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"***",cex=1.5)
        }

        if (pvalue[i,22] < 0.01) {
            segments(4.1,ma17,7.9,ma17)
            text(6.04,ma172,"**",cex=1.5)
        }

        if (pvalue[i,22] < 0.05) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"*",cex=1.5)
        }

        if (pvalue[i,23] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,23] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,23] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,24] < 0.001) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"***",cex=1.5)
        }

        if (pvalue[i,24] < 0.01) {
            segments(5.1,ma3,6.9,ma3)
            text(6.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,24] < 0.05) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"*",cex=1.5)
        }

        if (pvalue[i,25] < 0.001) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"***",cex=1.5)
        }

        if (pvalue[i,25] < 0.01) {
            segments(5.1,ma5,7.9,ma5)
            text(6.54,ma52,"**",cex=1.5)
        }

        if (pvalue[i,25] < 0.05) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"*",cex=1.5)
        }

        if (pvalue[i,26] < 0.001) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,26] < 0.01) {
            segments(6.1,ma2,6.9,ma2)
            text(6.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,26] < 0.05) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,27] < 0.001) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"***",cex=1.5)
        }

        if (pvalue[i,27] < 0.01) {
            segments(6.1,ma4,7.9,ma4)
            text(7.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,27] < 0.05) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"*",cex=1.5)
        }

        if (pvalue[i,28] < 0.001) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,28] < 0.01) {
            segments(7.1,ma2,7.9,ma2)
            text(7.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,28] < 0.05) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF2),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        #mtext("p=",3,at=5.2,cex=1.25)
        #mtext(round(c[i,1],3),3,at=5.5,cex=1.25)
        ma2=max(dataM[,i])+konst/17
        ma3=max(dataM[,i])+konst*2/17
        ma4=max(dataM[,i])+konst*3/17
        ma5=max(dataM[,i])+konst*4/17
        ma6=max(dataM[,i])+konst*5/17
        ma7=max(dataM[,i])+konst*6/17
        ma8=max(dataM[,i])+konst*7/17
        ma9=max(dataM[,i])+konst*8/17
        ma10=max(dataM[,i])+konst*9/17
        ma11=max(dataM[,i])+konst*10/17
        ma12=max(dataM[,i])+konst*11/17
        ma13=max(dataM[,i])+konst*12/17
        ma14=max(dataM[,i])+konst*13/17
        ma15=max(dataM[,i])+konst*14/17
        ma16=max(dataM[,i])+konst*15/17
        ma17=max(dataM[,i])+konst*16/17
        ma22=max(dataM[,i])+konst*1.25/16
        ma32=max(dataM[,i])+konst*2/16
        ma42=max(dataM[,i])+konst*3/16
        ma52=max(dataM[,i])+konst*4/16
        ma62=max(dataM[,i])+konst*5/16
        ma72=max(dataM[,i])+konst*6/16
        ma82=max(dataM[,i])+konst*7/16
        ma92=max(dataM[,i])+konst*7.5/16
        ma102=max(dataM[,i])+konst*8.5/16
        ma112=max(dataM[,i])+konst*9.5/16
        ma122=max(dataM[,i])+konst*10.5/16
        ma132=max(dataM[,i])+konst*11.5/16
        ma142=max(dataM[,i])+konst*12.5/16
        ma152=max(dataM[,i])+konst*13.5/16
        ma162=max(dataM[,i])+konst*14.5/16
        ma172=max(dataM[,i])+konst*15.5/16

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(1.1,ma11,6.9,ma11)
            text(4.04,ma112,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(1.1,ma14,7.9,ma14)
            text(4.54,ma142,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(2.1,ma12,6.9,ma12)
            text(4.54,ma122,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(2.1,ma15,7.9,ma15)
            text(5.04,ma152,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,16] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,16] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,16] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,17] < 0.001) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"***",cex=1.5)
        }

        if (pvalue[i,17] < 0.01) {
            segments(3.1,ma13,6.9,ma13)
            text(5.04,ma132,"**",cex=1.5)
        }

        if (pvalue[i,17] < 0.05) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"*",cex=1.5)
        }

        if (pvalue[i,18] < 0.001) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"***",cex=1.5)
        }

        if (pvalue[i,18] < 0.01) {
            segments(3.1,ma16,7.9,ma16)
            text(5.54,ma162,"**",cex=1.5)
        }

        if (pvalue[i,18] < 0.05) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"*",cex=1.5)
        }

        if (pvalue[i,19] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,19] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,19] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,20] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,20] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,20] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,21] < 0.001) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,21] < 0.01) {
            segments(4.1,ma7,6.9,ma7)
            text(5.54,ma72,"**",cex=1.5)
        }

        if (pvalue[i,21] < 0.05) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"*",cex=1.5)
        }

        if (pvalue[i,22] < 0.001) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"***",cex=1.5)
        }

        if (pvalue[i,22] < 0.01) {
            segments(4.1,ma17,7.9,ma17)
            text(6.04,ma172,"**",cex=1.5)
        }

        if (pvalue[i,22] < 0.05) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"*",cex=1.5)
        }

        if (pvalue[i,23] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,23] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,23] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,24] < 0.001) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"***",cex=1.5)
        }

        if (pvalue[i,24] < 0.01) {
            segments(5.1,ma3,6.9,ma3)
            text(6.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,24] < 0.05) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"*",cex=1.5)
        }

        if (pvalue[i,25] < 0.001) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"***",cex=1.5)
        }

        if (pvalue[i,25] < 0.01) {
            segments(5.1,ma5,7.9,ma5)
            text(6.54,ma52,"**",cex=1.5)
        }

        if (pvalue[i,25] < 0.05) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"*",cex=1.5)
        }

        if (pvalue[i,26] < 0.001) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,26] < 0.01) {
            segments(6.1,ma2,6.9,ma2)
            text(6.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,26] < 0.05) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,27] < 0.001) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"***",cex=1.5)
        }

        if (pvalue[i,27] < 0.01) {
            segments(6.1,ma4,7.9,ma4)
            text(7.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,27] < 0.05) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"*",cex=1.5)
        }

        if (pvalue[i,28] < 0.001) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,28] < 0.01) {
            segments(7.1,ma2,7.9,ma2)
            text(7.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,28] < 0.05) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF3),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/17
        ma3=max(dataM[,i])+konst*2/17
        ma4=max(dataM[,i])+konst*3/17
        ma5=max(dataM[,i])+konst*4/17
        ma6=max(dataM[,i])+konst*5/17
        ma7=max(dataM[,i])+konst*6/17
        ma8=max(dataM[,i])+konst*7/17
        ma9=max(dataM[,i])+konst*8/17
        ma10=max(dataM[,i])+konst*9/17
        ma11=max(dataM[,i])+konst*10/17
        ma12=max(dataM[,i])+konst*11/17
        ma13=max(dataM[,i])+konst*12/17
        ma14=max(dataM[,i])+konst*13/17
        ma15=max(dataM[,i])+konst*14/17
        ma16=max(dataM[,i])+konst*15/17
        ma17=max(dataM[,i])+konst*16/17
        ma22=max(dataM[,i])+konst*1.25/16
        ma32=max(dataM[,i])+konst*2/16
        ma42=max(dataM[,i])+konst*3/16
        ma52=max(dataM[,i])+konst*4/16
        ma62=max(dataM[,i])+konst*5/16
        ma72=max(dataM[,i])+konst*6/16
        ma82=max(dataM[,i])+konst*7/16
        ma92=max(dataM[,i])+konst*7.5/16
        ma102=max(dataM[,i])+konst*8.5/16
        ma112=max(dataM[,i])+konst*9.5/16
        ma122=max(dataM[,i])+konst*10.5/16
        ma132=max(dataM[,i])+konst*11.5/16
        ma142=max(dataM[,i])+konst*12.5/16
        ma152=max(dataM[,i])+konst*13.5/16
        ma162=max(dataM[,i])+konst*14.5/16
        ma172=max(dataM[,i])+konst*15.5/16

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,16] < 0.05) {
            if (pvalue[i,16]>0.0001){
                a=round(pvalue[i,16],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,16],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,17] < 0.05) {
            if (pvalue[i,17]>0.0001){
                a=round(pvalue[i,17],4)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }else{
                a=signif(pvalue[i,17],3)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }}

        if (pvalue[i,18] < 0.05) {
            if (pvalue[i,18]>0.0001){
                a=round(pvalue[i,18],4)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }else{
                a=signif(pvalue[i,18],3)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }}

        if (pvalue[i,19] < 0.05) {
            if (pvalue[i,19]>0.0001){
                a=round(pvalue[i,19],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,19],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,20] < 0.05) {
            if (pvalue[i,20]>0.0001){
                a=round(pvalue[i,20],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,20],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,21] < 0.05) {
            if (pvalue[i,21]>0.0001){
                a=round(pvalue[i,21],4)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,21],3)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }}

        if (pvalue[i,22] < 0.05) {
            if (pvalue[i,22]>0.0001){
                a=round(pvalue[i,22],4)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }else{
                a=signif(pvalue[i,22],3)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }}

        if (pvalue[i,23] < 0.05) {
            if (pvalue[i,23]>0.0001){
                a=round(pvalue[i,23],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,23],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,24] < 0.05) {
            if (pvalue[i,24]>0.0001){
                a=round(pvalue[i,24],4)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,24],3)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }}

        if (pvalue[i,25] < 0.05) {
            if (pvalue[i,25]>0.0001){
                a=round(pvalue[i,25],4)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,25],3)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }}

        if (pvalue[i,26] < 0.05) {
            if (pvalue[i,26]>0.0001){
                a=round(pvalue[i,26],4)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,26],3)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,27] < 0.05) {
            if (pvalue[i,27]>0.0001){
                a=round(pvalue[i,27],4)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,27],3)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }}

        if (pvalue[i,28] < 0.05) {
            if (pvalue[i,28]>0.0001){
                a=round(pvalue[i,28],4)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,28],3)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/17
        ma3=max(dataM[,i])+konst*2/17
        ma4=max(dataM[,i])+konst*3/17
        ma5=max(dataM[,i])+konst*4/17
        ma6=max(dataM[,i])+konst*5/17
        ma7=max(dataM[,i])+konst*6/17
        ma8=max(dataM[,i])+konst*7/17
        ma9=max(dataM[,i])+konst*8/17
        ma10=max(dataM[,i])+konst*9/17
        ma11=max(dataM[,i])+konst*10/17
        ma12=max(dataM[,i])+konst*11/17
        ma13=max(dataM[,i])+konst*12/17
        ma14=max(dataM[,i])+konst*13/17
        ma15=max(dataM[,i])+konst*14/17
        ma16=max(dataM[,i])+konst*15/17
        ma17=max(dataM[,i])+konst*16/17
        ma22=max(dataM[,i])+konst*1.25/16
        ma32=max(dataM[,i])+konst*2/16
        ma42=max(dataM[,i])+konst*3/16
        ma52=max(dataM[,i])+konst*4/16
        ma62=max(dataM[,i])+konst*5/16
        ma72=max(dataM[,i])+konst*6/16
        ma82=max(dataM[,i])+konst*7/16
        ma92=max(dataM[,i])+konst*7.5/16
        ma102=max(dataM[,i])+konst*8.5/16
        ma112=max(dataM[,i])+konst*9.5/16
        ma122=max(dataM[,i])+konst*10.5/16
        ma132=max(dataM[,i])+konst*11.5/16
        ma142=max(dataM[,i])+konst*12.5/16
        ma152=max(dataM[,i])+konst*13.5/16
        ma162=max(dataM[,i])+konst*14.5/16
        ma172=max(dataM[,i])+konst*15.5/16

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,16] < 0.05) {
            if (pvalue[i,16]>0.0001){
                a=round(pvalue[i,16],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,16],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,17] < 0.05) {
            if (pvalue[i,17]>0.0001){
                a=round(pvalue[i,17],4)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }else{
                a=signif(pvalue[i,17],3)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }}

        if (pvalue[i,18] < 0.05) {
            if (pvalue[i,18]>0.0001){
                a=round(pvalue[i,18],4)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }else{
                a=signif(pvalue[i,18],3)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }}

        if (pvalue[i,19] < 0.05) {
            if (pvalue[i,19]>0.0001){
                a=round(pvalue[i,19],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,19],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,20] < 0.05) {
            if (pvalue[i,20]>0.0001){
                a=round(pvalue[i,20],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,20],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,21] < 0.05) {
            if (pvalue[i,21]>0.0001){
                a=round(pvalue[i,21],4)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,21],3)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }}

        if (pvalue[i,22] < 0.05) {
            if (pvalue[i,22]>0.0001){
                a=round(pvalue[i,22],4)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }else{
                a=signif(pvalue[i,22],3)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }}

        if (pvalue[i,23] < 0.05) {
            if (pvalue[i,23]>0.0001){
                a=round(pvalue[i,23],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,23],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,24] < 0.05) {
            if (pvalue[i,24]>0.0001){
                a=round(pvalue[i,24],4)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,24],3)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }}

        if (pvalue[i,25] < 0.05) {
            if (pvalue[i,25]>0.0001){
                a=round(pvalue[i,25],4)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,25],3)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }}

        if (pvalue[i,26] < 0.05) {
            if (pvalue[i,26]>0.0001){
                a=round(pvalue[i,26],4)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,26],3)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,27] < 0.05) {
            if (pvalue[i,27]>0.0001){
                a=round(pvalue[i,27],4)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,27],3)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }}

        if (pvalue[i,28] < 0.05) {
            if (pvalue[i,28]>0.0001){
                a=round(pvalue[i,28],4)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,28],3)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

    }
##################################################################################
if (count == 9){

    pdf((PDF1),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/21
        ma3=max(dataM[,i])+konst*2/21
        ma4=max(dataM[,i])+konst*3/21
        ma5=max(dataM[,i])+konst*4/21
        ma6=max(dataM[,i])+konst*5/21
        ma7=max(dataM[,i])+konst*6/21
        ma8=max(dataM[,i])+konst*7/21
        ma9=max(dataM[,i])+konst*8/21
        ma10=max(dataM[,i])+konst*9/21
        ma11=max(dataM[,i])+konst*10/21
        ma12=max(dataM[,i])+konst*11/21
        ma13=max(dataM[,i])+konst*12/21
        ma14=max(dataM[,i])+konst*13/21
        ma15=max(dataM[,i])+konst*14/21
        ma16=max(dataM[,i])+konst*15/21
        ma17=max(dataM[,i])+konst*16/21
        ma18=max(dataM[,i])+konst*17/21
        ma19=max(dataM[,i])+konst*18/21
        ma20=max(dataM[,i])+konst*19/21
        ma21=max(dataM[,i])+konst*20/21
        ma22=max(dataM[,i])+konst*1.25/20
        ma32=max(dataM[,i])+konst*2/20
        ma42=max(dataM[,i])+konst*3/20
        ma52=max(dataM[,i])+konst*4/20
        ma62=max(dataM[,i])+konst*5/20
        ma72=max(dataM[,i])+konst*6/20
        ma82=max(dataM[,i])+konst*7/20
        ma92=max(dataM[,i])+konst*7.5/20
        ma102=max(dataM[,i])+konst*8.5/20
        ma112=max(dataM[,i])+konst*9.5/20
        ma122=max(dataM[,i])+konst*10.5/20
        ma132=max(dataM[,i])+konst*11.5/20
        ma142=max(dataM[,i])+konst*12.5/20
        ma152=max(dataM[,i])+konst*13.5/20
        ma162=max(dataM[,i])+konst*14.5/20
        ma172=max(dataM[,i])+konst*15.5/20
        ma182=max(dataM[,i])+konst*16.5/20
        ma192=max(dataM[,i])+konst*17.5/20
        ma202=max(dataM[,i])+konst*18.5/20
        ma212=max(dataM[,i])+konst*19.5/20

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(1.1,ma11,6.9,ma11)
            text(4.04,ma112,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(1.1,ma14,7.9,ma14)
            text(4.54,ma142,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(1.1,ma18,8.9,ma18)
            text(5,ma182,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(1.1,ma18,8.9,ma18)
            text(5.04,ma182,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(1.1,ma18,8.9,ma18)
            text(5,ma182,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(2.1,ma12,6.9,ma12)
            text(4.54,ma122,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(2.1,ma15,7.9,ma15)
            text(5.04,ma152,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(2.1,ma19,8.9,ma19)
            text(5.5,ma192,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(2.1,ma19,8.9,ma19)
            text(5.54,ma192,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(2.1,ma19,8.9,ma19)
            text(5.5,ma192,"*",cex=1.5)
        }

        if (pvalue[i,16] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,16] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,16] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,17] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,17] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,17] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,18] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,18] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,18] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,19] < 0.001) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"***",cex=1.5)
        }

        if (pvalue[i,19] < 0.01) {
            segments(3.1,ma13,6.9,ma13)
            text(5.04,ma132,"**",cex=1.5)
        }

        if (pvalue[i,19] < 0.05) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"*",cex=1.5)
        }

        if (pvalue[i,20] < 0.001) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"***",cex=1.5)
        }

        if (pvalue[i,20] < 0.01) {
            segments(3.1,ma16,7.9,ma16)
            text(5.54,ma162,"**",cex=1.5)
        }

        if (pvalue[i,20] < 0.05) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"*",cex=1.5)
        }

        if (pvalue[i,21] < 0.001) {
            segments(3.1,ma20,8.9,ma20)
            text(6,ma202,"***",cex=1.5)
        }

        if (pvalue[i,21] < 0.01) {
            segments(3.1,ma20,8.9,ma20)
            text(6.04,ma202,"**",cex=1.5)
        }

        if (pvalue[i,21] < 0.05) {
            segments(3.1,ma20,8.9,ma20)
            text(6,ma202,"*",cex=1.5)
        }

        if (pvalue[i,22] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,22] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,22] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,23] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,23] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,23] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,24] < 0.001) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,24] < 0.01) {
            segments(4.1,ma7,6.9,ma7)
            text(5.54,ma72,"**",cex=1.5)
        }

        if (pvalue[i,24] < 0.05) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"*",cex=1.5)
        }

        if (pvalue[i,25] < 0.001) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"***",cex=1.5)
        }

        if (pvalue[i,25] < 0.01) {
            segments(4.1,ma17,7.9,ma17)
            text(6.04,ma172,"**",cex=1.5)
        }

        if (pvalue[i,25] < 0.05) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"*",cex=1.5)
        }

        if (pvalue[i,26] < 0.001) {
            segments(4.1,ma21,8.9,ma21)
            text(6.5,ma212,"***",cex=1.5)
        }

        if (pvalue[i,26] < 0.01) {
            segments(4.1,ma21,8.9,ma21)
            text(6.54,ma212,"**",cex=1.5)
        }

        if (pvalue[i,26] < 0.05) {
            segments(4.1,ma21,8.9,ma21)
            text(6.5,ma212,"*",cex=1.5)
        }

        if (pvalue[i,27] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,27] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,27] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,28] < 0.001) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"***",cex=1.5)
        }

        if (pvalue[i,28] < 0.01) {
            segments(5.1,ma3,6.9,ma3)
            text(6.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,28] < 0.05) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"*",cex=1.5)
        }

        if (pvalue[i,29] < 0.001) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"***",cex=1.5)
        }

        if (pvalue[i,29] < 0.01) {
            segments(5.1,ma5,7.9,ma5)
            text(6.54,ma52,"**",cex=1.5)
        }

        if (pvalue[i,29] < 0.05) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"*",cex=1.5)
        }

        if (pvalue[i,30] < 0.001) {
            segments(5.1,ma8,8.9,ma8)
            text(7,ma82,"***",cex=1.5)
        }

        if (pvalue[i,30] < 0.01) {
            segments(5.1,ma8,8.9,ma8)
            text(7.04,ma82,"**",cex=1.5)
        }

        if (pvalue[i,30] < 0.05) {
            segments(5.1,ma8,8.9,ma8)
            text(7,ma82,"*",cex=1.5)
        }

        if (pvalue[i,31] < 0.001) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,31] < 0.01) {
            segments(6.1,ma2,6.9,ma2)
            text(6.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,31] < 0.05) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,32] < 0.001) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"***",cex=1.5)
        }

        if (pvalue[i,32] < 0.01) {
            segments(6.1,ma4,7.9,ma4)
            text(7.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,32] < 0.05) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"*",cex=1.5)
        }

        if (pvalue[i,33] < 0.001) {
            segments(6.1,ma10,8.9,ma10)
            text(7.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,33] < 0.01) {
            segments(6.1,ma10,8.9,ma10)
            text(7.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,33] < 0.05) {
            segments(6.1,ma10,8.9,ma10)
            text(7.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,34] < 0.001) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,34] < 0.01) {
            segments(7.1,ma2,7.9,ma2)
            text(7.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,34] < 0.05) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,35] < 0.001) {
            segments(7.1,ma7,8.9,ma7)
            text(8,ma72,"***",cex=1.5)
        }

        if (pvalue[i,35] < 0.01) {
            segments(7.1,ma7,8.9,ma7)
            text(8.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,35] < 0.05) {
            segments(7.1,ma7,8.9,ma7)
            text(8,ma72,"*",cex=1.5)
        }

        if (pvalue[i,36] < 0.001) {
            segments(8.1,ma2,8.9,ma2)
            text(8.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,36] < 0.01) {
            segments(8.1,ma2,8.9,ma2)
            text(8.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,36] < 0.05) {
            segments(8.1,ma2,8.9,ma2)
            text(8.5,ma22,"*",cex=1.5)
        }

    }
    dev.off()

    pdf((PDF2),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/21
        ma3=max(dataM[,i])+konst*2/21
        ma4=max(dataM[,i])+konst*3/21
        ma5=max(dataM[,i])+konst*4/21
        ma6=max(dataM[,i])+konst*5/21
        ma7=max(dataM[,i])+konst*6/21
        ma8=max(dataM[,i])+konst*7/21
        ma9=max(dataM[,i])+konst*8/21
        ma10=max(dataM[,i])+konst*9/21
        ma11=max(dataM[,i])+konst*10/21
        ma12=max(dataM[,i])+konst*11/21
        ma13=max(dataM[,i])+konst*12/21
        ma14=max(dataM[,i])+konst*13/21
        ma15=max(dataM[,i])+konst*14/21
        ma16=max(dataM[,i])+konst*15/21
        ma17=max(dataM[,i])+konst*16/21
        ma18=max(dataM[,i])+konst*17/21
        ma19=max(dataM[,i])+konst*18/21
        ma20=max(dataM[,i])+konst*19/21
        ma21=max(dataM[,i])+konst*20/21
        ma22=max(dataM[,i])+konst*1.25/20
        ma32=max(dataM[,i])+konst*2/20
        ma42=max(dataM[,i])+konst*3/20
        ma52=max(dataM[,i])+konst*4/20
        ma62=max(dataM[,i])+konst*5/20
        ma72=max(dataM[,i])+konst*6/20
        ma82=max(dataM[,i])+konst*7/20
        ma92=max(dataM[,i])+konst*7.5/20
        ma102=max(dataM[,i])+konst*8.5/20
        ma112=max(dataM[,i])+konst*9.5/20
        ma122=max(dataM[,i])+konst*10.5/20
        ma132=max(dataM[,i])+konst*11.5/20
        ma142=max(dataM[,i])+konst*12.5/20
        ma152=max(dataM[,i])+konst*13.5/20
        ma162=max(dataM[,i])+konst*14.5/20
        ma172=max(dataM[,i])+konst*15.5/20
        ma182=max(dataM[,i])+konst*16.5/20
        ma192=max(dataM[,i])+konst*17.5/20
        ma202=max(dataM[,i])+konst*18.5/20
        ma212=max(dataM[,i])+konst*19.5/20

        if (pvalue[i,1] < 0.001) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,1] < 0.01) {
            segments(1.1,ma2,1.9,ma2)
            text(1.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,1] < 0.05) {
            segments(1.1,ma2,1.9,ma2)
            text(1.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,2] < 0.001) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"***",cex=1.5)
        }

        if (pvalue[i,2] < 0.01) {
            segments(1.1,ma3,2.9,ma3)
            text(2.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,2] < 0.05) {
            segments(1.1,ma3,2.9,ma3)
            text(2,ma32,"*",cex=1.5)
        }

        if (pvalue[i,3] < 0.001) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,3] < 0.01) {
            segments(1.1,ma4,3.9,ma4)
            text(2.54,ma42,"**",cex=1.5)
        }

        if (pvalue[i,3] < 0.05) {
            segments(1.1,ma4,3.9,ma4)
            text(2.5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,4] < 0.001) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"***",cex=1.5)
        }

        if (pvalue[i,4] < 0.01) {
            segments(1.1,ma5,4.9,ma5)
            text(3.04,ma52,"**",cex=1.5)
        }

        if (pvalue[i,4] < 0.05) {
            segments(1.1,ma5,4.9,ma5)
            text(3,ma52,"*",cex=1.5)
        }

        if (pvalue[i,5] < 0.001) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"***",cex=1.5)
        }

        if (pvalue[i,5] < 0.01) {
            segments(1.1,ma6,5.9,ma6)
            text(3.54,ma62,"**",cex=1.5)
        }

        if (pvalue[i,5] < 0.05) {
            segments(1.1,ma6,5.9,ma6)
            text(3.5,ma62,"*",cex=1.5)
        }

        if (pvalue[i,6] < 0.001) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"***",cex=1.5)
        }

        if (pvalue[i,6] < 0.01) {
            segments(1.1,ma11,6.9,ma11)
            text(4.04,ma112,"**",cex=1.5)
        }

        if (pvalue[i,6] < 0.05) {
            segments(1.1,ma11,6.9,ma11)
            text(4,ma112,"*",cex=1.5)
        }

        if (pvalue[i,7] < 0.001) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"***",cex=1.5)
        }

        if (pvalue[i,7] < 0.01) {
            segments(1.1,ma14,7.9,ma14)
            text(4.54,ma142,"**",cex=1.5)
        }

        if (pvalue[i,7] < 0.05) {
            segments(1.1,ma14,7.9,ma14)
            text(4.5,ma142,"*",cex=1.5)
        }

        if (pvalue[i,8] < 0.001) {
            segments(1.1,ma18,8.9,ma18)
            text(5,ma182,"***",cex=1.5)
        }

        if (pvalue[i,8] < 0.01) {
            segments(1.1,ma18,8.9,ma18)
            text(5.04,ma182,"**",cex=1.5)
        }

        if (pvalue[i,8] < 0.05) {
            segments(1.1,ma18,8.9,ma18)
            text(5,ma182,"*",cex=1.5)
        }

        if (pvalue[i,9] < 0.001) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,9] < 0.01) {
            segments(2.1,ma2,2.9,ma2)
            text(2.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,9] < 0.05) {
            segments(2.1,ma2,2.9,ma2)
            text(2.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,10] < 0.001) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"***",cex=1.5)
        }

        if (pvalue[i,10] < 0.01) {
            segments(2.1,ma7,3.9,ma7)
            text(3.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,10] < 0.05) {
            segments(2.1,ma7,3.9,ma7)
            text(3,ma72,"*",cex=1.5)
        }

        if (pvalue[i,11] < 0.001) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"***",cex=1.5)
        }

        if (pvalue[i,11] < 0.01) {
            segments(2.1,ma8,4.9,ma8)
            text(3.54,ma82,"**",cex=1.5)
        }

        if (pvalue[i,11] < 0.05) {
            segments(2.1,ma8,4.9,ma8)
            text(3.5,ma82,"*",cex=1.5)
        }

        if (pvalue[i,12] < 0.001) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"***",cex=1.5)
        }

        if (pvalue[i,12] < 0.01) {
            segments(2.1,ma9,5.9,ma9)
            text(4.04,ma92,"**",cex=1.5)
        }

        if (pvalue[i,12] < 0.05) {
            segments(2.1,ma9,5.9,ma9)
            text(4,ma92,"*",cex=1.5)
        }

        if (pvalue[i,13] < 0.001) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"***",cex=1.5)
        }

        if (pvalue[i,13] < 0.01) {
            segments(2.1,ma12,6.9,ma12)
            text(4.54,ma122,"**",cex=1.5)
        }

        if (pvalue[i,13] < 0.05) {
            segments(2.1,ma12,6.9,ma12)
            text(4.5,ma122,"*",cex=1.5)
        }

        if (pvalue[i,14] < 0.001) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"***",cex=1.5)
        }

        if (pvalue[i,14] < 0.01) {
            segments(2.1,ma15,7.9,ma15)
            text(5.04,ma152,"**",cex=1.5)
        }

        if (pvalue[i,14] < 0.05) {
            segments(2.1,ma15,7.9,ma15)
            text(5,ma152,"*",cex=1.5)
        }

        if (pvalue[i,15] < 0.001) {
            segments(2.1,ma19,8.9,ma19)
            text(5.5,ma192,"***",cex=1.5)
        }

        if (pvalue[i,15] < 0.01) {
            segments(2.1,ma19,8.9,ma19)
            text(5.54,ma192,"**",cex=1.5)
        }

        if (pvalue[i,15] < 0.05) {
            segments(2.1,ma19,8.9,ma19)
            text(5.5,ma192,"*",cex=1.5)
        }

        if (pvalue[i,16] < 0.001) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,16] < 0.01) {
            segments(3.1,ma2,3.9,ma2)
            text(3.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,16] < 0.05) {
            segments(3.1,ma2,3.9,ma2)
            text(3.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,17] < 0.001) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"***",cex=1.5)
        }

        if (pvalue[i,17] < 0.01) {
            segments(3.1,ma3,4.9,ma3)
            text(4.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,17] < 0.05) {
            segments(3.1,ma3,4.9,ma3)
            text(4,ma32,"*",cex=1.5)
        }

        if (pvalue[i,18] < 0.001) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,18] < 0.01) {
            segments(3.1,ma10,5.9,ma10)
            text(4.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,18] < 0.05) {
            segments(3.1,ma10,5.9,ma10)
            text(4.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,19] < 0.001) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"***",cex=1.5)
        }

        if (pvalue[i,19] < 0.01) {
            segments(3.1,ma13,6.9,ma13)
            text(5.04,ma132,"**",cex=1.5)
        }

        if (pvalue[i,19] < 0.05) {
            segments(3.1,ma13,6.9,ma13)
            text(5,ma132,"*",cex=1.5)
        }

        if (pvalue[i,20] < 0.001) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"***",cex=1.5)
        }

        if (pvalue[i,20] < 0.01) {
            segments(3.1,ma16,7.9,ma16)
            text(5.54,ma162,"**",cex=1.5)
        }

        if (pvalue[i,20] < 0.05) {
            segments(3.1,ma16,7.9,ma16)
            text(5.5,ma162,"*",cex=1.5)
        }

        if (pvalue[i,21] < 0.001) {
            segments(3.1,ma20,8.9,ma20)
            text(6,ma202,"***",cex=1.5)
        }

        if (pvalue[i,21] < 0.01) {
            segments(3.1,ma20,8.9,ma20)
            text(6.04,ma202,"**",cex=1.5)
        }

        if (pvalue[i,21] < 0.05) {
            segments(3.1,ma20,8.9,ma20)
            text(6,ma202,"*",cex=1.5)
        }

        if (pvalue[i,22] < 0.001) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,22] < 0.01) {
            segments(4.1,ma2,4.9,ma2)
            text(4.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,22] < 0.05) {
            segments(4.1,ma2,4.9,ma2)
            text(4.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,23] < 0.001) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"***",cex=1.5)
        }

        if (pvalue[i,23] < 0.01) {
            segments(4.1,ma4,5.9,ma4)
            text(5.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,23] < 0.05) {
            segments(4.1,ma4,5.9,ma4)
            text(5,ma42,"*",cex=1.5)
        }

        if (pvalue[i,24] < 0.001) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"***",cex=1.5)
        }

        if (pvalue[i,24] < 0.01) {
            segments(4.1,ma7,6.9,ma7)
            text(5.54,ma72,"**",cex=1.5)
        }

        if (pvalue[i,24] < 0.05) {
            segments(4.1,ma7,6.9,ma7)
            text(5.5,ma72,"*",cex=1.5)
        }

        if (pvalue[i,25] < 0.001) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"***",cex=1.5)
        }

        if (pvalue[i,25] < 0.01) {
            segments(4.1,ma17,7.9,ma17)
            text(6.04,ma172,"**",cex=1.5)
        }

        if (pvalue[i,25] < 0.05) {
            segments(4.1,ma17,7.9,ma17)
            text(6,ma172,"*",cex=1.5)
        }

        if (pvalue[i,26] < 0.001) {
            segments(4.1,ma21,8.9,ma21)
            text(6.5,ma212,"***",cex=1.5)
        }

        if (pvalue[i,26] < 0.01) {
            segments(4.1,ma21,8.9,ma21)
            text(6.54,ma212,"**",cex=1.5)
        }

        if (pvalue[i,26] < 0.05) {
            segments(4.1,ma21,8.9,ma21)
            text(6.5,ma212,"*",cex=1.5)
        }

        if (pvalue[i,27] < 0.001) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,27] < 0.01) {
            segments(5.1,ma2,5.9,ma2)
            text(5.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,27] < 0.05) {
            segments(5.1,ma2,5.9,ma2)
            text(5.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,28] < 0.001) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"***",cex=1.5)
        }

        if (pvalue[i,28] < 0.01) {
            segments(5.1,ma3,6.9,ma3)
            text(6.04,ma32,"**",cex=1.5)
        }

        if (pvalue[i,28] < 0.05) {
            segments(5.1,ma3,6.9,ma3)
            text(6,ma32,"*",cex=1.5)
        }

        if (pvalue[i,29] < 0.001) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"***",cex=1.5)
        }

        if (pvalue[i,29] < 0.01) {
            segments(5.1,ma5,7.9,ma5)
            text(6.54,ma52,"**",cex=1.5)
        }

        if (pvalue[i,29] < 0.05) {
            segments(5.1,ma5,7.9,ma5)
            text(6.5,ma52,"*",cex=1.5)
        }

        if (pvalue[i,30] < 0.001) {
            segments(5.1,ma8,8.9,ma8)
            text(7,ma82,"***",cex=1.5)
        }

        if (pvalue[i,30] < 0.01) {
            segments(5.1,ma8,8.9,ma8)
            text(7.04,ma82,"**",cex=1.5)
        }

        if (pvalue[i,30] < 0.05) {
            segments(5.1,ma8,8.9,ma8)
            text(7,ma82,"*",cex=1.5)
        }

        if (pvalue[i,31] < 0.001) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,31] < 0.01) {
            segments(6.1,ma2,6.9,ma2)
            text(6.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,31] < 0.05) {
            segments(6.1,ma2,6.9,ma2)
            text(6.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,32] < 0.001) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"***",cex=1.5)
        }

        if (pvalue[i,32] < 0.01) {
            segments(6.1,ma4,7.9,ma4)
            text(7.04,ma42,"**",cex=1.5)
        }

        if (pvalue[i,32] < 0.05) {
            segments(6.1,ma4,7.9,ma4)
            text(7,ma42,"*",cex=1.5)
        }

        if (pvalue[i,33] < 0.001) {
            segments(6.1,ma10,8.9,ma10)
            text(7.5,ma102,"***",cex=1.5)
        }

        if (pvalue[i,33] < 0.01) {
            segments(6.1,ma10,8.9,ma10)
            text(7.54,ma102,"**",cex=1.5)
        }

        if (pvalue[i,33] < 0.05) {
            segments(6.1,ma10,8.9,ma10)
            text(7.5,ma102,"*",cex=1.5)
        }

        if (pvalue[i,34] < 0.001) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,34] < 0.01) {
            segments(7.1,ma2,7.9,ma2)
            text(7.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,34] < 0.05) {
            segments(7.1,ma2,7.9,ma2)
            text(7.5,ma22,"*",cex=1.5)
        }

        if (pvalue[i,35] < 0.001) {
            segments(7.1,ma7,8.9,ma7)
            text(8,ma72,"***",cex=1.5)
        }

        if (pvalue[i,35] < 0.01) {
            segments(7.1,ma7,8.9,ma7)
            text(8.04,ma72,"**",cex=1.5)
        }

        if (pvalue[i,35] < 0.05) {
            segments(7.1,ma7,8.9,ma7)
            text(8,ma72,"*",cex=1.5)
        }

        if (pvalue[i,36] < 0.001) {
            segments(8.1,ma2,8.9,ma2)
            text(8.5,ma22,"***",cex=1.5)
        }

        if (pvalue[i,36] < 0.01) {
            segments(8.1,ma2,8.9,ma2)
            text(8.54,ma22,"**",cex=1.5)
        }

        if (pvalue[i,36] < 0.05) {
            segments(8.1,ma2,8.9,ma2)
            text(8.5,ma22,"*",cex=1.5)
        }
    }
    dev.off()

    pdf((PDF3),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma))
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/21
        ma3=max(dataM[,i])+konst*2/21
        ma4=max(dataM[,i])+konst*3/21
        ma5=max(dataM[,i])+konst*4/21
        ma6=max(dataM[,i])+konst*5/21
        ma7=max(dataM[,i])+konst*6/21
        ma8=max(dataM[,i])+konst*7/21
        ma9=max(dataM[,i])+konst*8/21
        ma10=max(dataM[,i])+konst*9/21
        ma11=max(dataM[,i])+konst*10/21
        ma12=max(dataM[,i])+konst*11/21
        ma13=max(dataM[,i])+konst*12/21
        ma14=max(dataM[,i])+konst*13/21
        ma15=max(dataM[,i])+konst*14/21
        ma16=max(dataM[,i])+konst*15/21
        ma17=max(dataM[,i])+konst*16/21
        ma18=max(dataM[,i])+konst*17/21
        ma19=max(dataM[,i])+konst*18/21
        ma20=max(dataM[,i])+konst*19/21
        ma21=max(dataM[,i])+konst*20/21
        ma22=max(dataM[,i])+konst*1.25/20
        ma32=max(dataM[,i])+konst*2/20
        ma42=max(dataM[,i])+konst*3/20
        ma52=max(dataM[,i])+konst*4/20
        ma62=max(dataM[,i])+konst*5/20
        ma72=max(dataM[,i])+konst*6/20
        ma82=max(dataM[,i])+konst*7/20
        ma92=max(dataM[,i])+konst*7.5/20
        ma102=max(dataM[,i])+konst*8.5/20
        ma112=max(dataM[,i])+konst*9.5/20
        ma122=max(dataM[,i])+konst*10.5/20
        ma132=max(dataM[,i])+konst*11.5/20
        ma142=max(dataM[,i])+konst*12.5/20
        ma152=max(dataM[,i])+konst*13.5/20
        ma162=max(dataM[,i])+konst*14.5/20
        ma172=max(dataM[,i])+konst*15.5/20
        ma182=max(dataM[,i])+konst*16.5/20
        ma192=max(dataM[,i])+konst*17.5/20
        ma202=max(dataM[,i])+konst*18.5/20
        ma212=max(dataM[,i])+konst*19.5/20

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(1.1,ma18,8.9,ma18)
                text(4.75,ma182,"p=",cex=0.75)
                text(5,ma182,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(1.1,ma18,8.9,ma18)
                text(4.75,ma182,"p=",cex=0.75)
                text(5,ma182,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(2.1,ma19,8.9,ma19)
                text(5.25,ma192,"p=",cex=0.75)
                text(5.5,ma192,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(2.1,ma19,8.9,ma19)
                text(5.25,ma192,"p=",cex=0.75)
                text(5.5,ma192,a,cex=0.75)
            }}

        if (pvalue[i,16] < 0.05) {
            if (pvalue[i,16]>0.0001){
                a=round(pvalue[i,16],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,16],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,17] < 0.05) {
            if (pvalue[i,17]>0.0001){
                a=round(pvalue[i,17],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,17],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,18] < 0.05) {
            if (pvalue[i,18]>0.0001){
                a=round(pvalue[i,18],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,18],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,19] < 0.05) {
            if (pvalue[i,19]>0.0001){
                a=round(pvalue[i,19],4)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }else{
                a=signif(pvalue[i,19],3)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }}

        if (pvalue[i,20] < 0.05) {
            if (pvalue[i,20]>0.0001){
                a=round(pvalue[i,20],4)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }else{
                a=signif(pvalue[i,20],3)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }}

        if (pvalue[i,21] < 0.05) {
            if (pvalue[i,21]>0.0001){
                a=round(pvalue[i,21],4)
                segments(3.1,ma20,8.9,ma20)
                text(5.75,ma202,"p=",cex=0.75)
                text(6,ma202,a,cex=0.75)
            }else{
                a=signif(pvalue[i,21],3)
                segments(3.1,ma20,8.9,ma20)
                text(5.75,ma202,"p=",cex=0.75)
                text(6,ma202,a,cex=0.75)
            }}

        if (pvalue[i,22] < 0.05) {
            if (pvalue[i,22]>0.0001){
                a=round(pvalue[i,22],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,22],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,23] < 0.05) {
            if (pvalue[i,23]>0.0001){
                a=round(pvalue[i,23],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,23],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,24] < 0.05) {
            if (pvalue[i,24]>0.0001){
                a=round(pvalue[i,24],4)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,24],3)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }}

        if (pvalue[i,25] < 0.05) {
            if (pvalue[i,25]>0.0001){
                a=round(pvalue[i,25],4)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }else{
                a=signif(pvalue[i,25],3)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }}

        if (pvalue[i,26] < 0.05) {
            if (pvalue[i,26]>0.0001){
                a=round(pvalue[i,26],4)
                segments(4.1,ma21,8.9,ma21)
                text(6.25,ma212,"p=",cex=0.75)
                text(6.5,ma212,a,cex=0.75)
            }else{
                a=signif(pvalue[i,26],3)
                segments(4.1,ma21,8.9,ma21)
                text(6.25,ma212,"p=",cex=0.75)
                text(6.5,ma212,a,cex=0.75)
            }}

        if (pvalue[i,27] < 0.05) {
            if (pvalue[i,27]>0.0001){
                a=round(pvalue[i,27],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,27],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,28] < 0.05) {
            if (pvalue[i,28]>0.0001){
                a=round(pvalue[i,28],4)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,28],3)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }}

        if (pvalue[i,29] < 0.05) {
            if (pvalue[i,29]>0.0001){
                a=round(pvalue[i,29],4)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,29],3)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }}

        if (pvalue[i,30] < 0.05) {
            if (pvalue[i,30]>0.0001){
                a=round(pvalue[i,30],4)
                segments(5.1,ma8,8.9,ma8)
                text(6.75,ma82,"p=",cex=0.75)
                text(7,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,30],3)
                segments(5.1,ma8,8.9,ma8)
                text(6.75,ma82,"p=",cex=0.75)
                text(7,ma82,a,cex=0.75)
            }}

        if (pvalue[i,31] < 0.05) {
            if (pvalue[i,31]>0.0001){
                a=round(pvalue[i,31],4)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,31],3)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,32] < 0.05) {
            if (pvalue[i,32]>0.0001){
                a=round(pvalue[i,32],4)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,32],3)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }}

        if (pvalue[i,33] < 0.05) {
            if (pvalue[i,33]>0.0001){
                a=round(pvalue[i,33],4)
                segments(6.1,ma10,8.9,ma10)
                text(7.25,ma102,"p=",cex=0.75)
                text(7.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,33],3)
                segments(6.1,ma10,8.9,ma10)
                text(7.25,ma102,"p=",cex=0.75)
                text(7.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,34] < 0.05) {
            if (pvalue[i,34]>0.0001){
                a=round(pvalue[i,34],4)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,34],3)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,35] < 0.05) {
            if (pvalue[i,35]>0.0001){
                a=round(pvalue[i,35],4)
                segments(7.1,ma7,8.9,ma7)
                text(7.75,ma72,"p=",cex=0.75)
                text(8,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,35],3)
                segments(7.1,ma7,8.9,ma7)
                text(7.75,ma72,"p=",cex=0.75)
                text(8,ma72,a,cex=0.75)
            }}

        if (pvalue[i,36] < 0.05) {
            if (pvalue[i,36]>0.0001){
                a=round(pvalue[i,36],4)
                segments(8.1,ma2,8.9,ma2)
                text(8.25,ma22,"p=",cex=0.75)
                text(8.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,36],3)
                segments(8.1,ma2,8.9,ma2)
                text(8.25,ma22,"p=",cex=0.75)
                text(8.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

    pdf((PDF4),width=20,height=10)
    par(mar = c(4, 4, 6, 4) + 0.1,oma=c(1,1,1,1))

    for(i in 1:ncol(dataM)){
        data=dataM[,i]
        groups=groups
        ma=max(dataM[,i])
        mi=min(dataM[,i])
        konst=abs(ma-mi)
        ma=max(dataM[,i])+konst
        mi=min(dataM[,i])
        boxplot(data ~ groups, names=groupnames,main=names[i],outpch = NA,cex.axis=1.25,ylim=c(mi,ma),notch=TRUE)
        stripchart(data ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        ma2=max(dataM[,i])+konst/21
        ma3=max(dataM[,i])+konst*2/21
        ma4=max(dataM[,i])+konst*3/21
        ma5=max(dataM[,i])+konst*4/21
        ma6=max(dataM[,i])+konst*5/21
        ma7=max(dataM[,i])+konst*6/21
        ma8=max(dataM[,i])+konst*7/21
        ma9=max(dataM[,i])+konst*8/21
        ma10=max(dataM[,i])+konst*9/21
        ma11=max(dataM[,i])+konst*10/21
        ma12=max(dataM[,i])+konst*11/21
        ma13=max(dataM[,i])+konst*12/21
        ma14=max(dataM[,i])+konst*13/21
        ma15=max(dataM[,i])+konst*14/21
        ma16=max(dataM[,i])+konst*15/21
        ma17=max(dataM[,i])+konst*16/21
        ma18=max(dataM[,i])+konst*17/21
        ma19=max(dataM[,i])+konst*18/21
        ma20=max(dataM[,i])+konst*19/21
        ma21=max(dataM[,i])+konst*20/21
        ma22=max(dataM[,i])+konst*1.25/20
        ma32=max(dataM[,i])+konst*2/20
        ma42=max(dataM[,i])+konst*3/20
        ma52=max(dataM[,i])+konst*4/20
        ma62=max(dataM[,i])+konst*5/20
        ma72=max(dataM[,i])+konst*6/20
        ma82=max(dataM[,i])+konst*7/20
        ma92=max(dataM[,i])+konst*7.5/20
        ma102=max(dataM[,i])+konst*8.5/20
        ma112=max(dataM[,i])+konst*9.5/20
        ma122=max(dataM[,i])+konst*10.5/20
        ma132=max(dataM[,i])+konst*11.5/20
        ma142=max(dataM[,i])+konst*12.5/20
        ma152=max(dataM[,i])+konst*13.5/20
        ma162=max(dataM[,i])+konst*14.5/20
        ma172=max(dataM[,i])+konst*15.5/20
        ma182=max(dataM[,i])+konst*16.5/20
        ma192=max(dataM[,i])+konst*17.5/20
        ma202=max(dataM[,i])+konst*18.5/20
        ma212=max(dataM[,i])+konst*19.5/20

        if (pvalue[i,1] < 0.05) {
            if (pvalue[i,1]>0.0001){
                a=round(pvalue[i,1],4)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,1],3)
                segments(1.1,ma2,1.9,ma2)
                text(1.25,ma22,"p=",cex=0.75)
                text(1.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,2] < 0.05) {
            if (pvalue[i,2]>0.0001){
                a=round(pvalue[i,2],4)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,2],3)
                segments(1.1,ma3,2.9,ma3)
                text(1.75,ma32,"p=",cex=0.75)
                text(2,ma32,a,cex=0.75)
            }}

        if (pvalue[i,3] < 0.05) {
            if (pvalue[i,3]>0.0001){
                a=round(pvalue[i,3],4)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,3],3)
                segments(1.1,ma4,3.9,ma4)
                text(2.25,ma42,"p=",cex=0.75)
                text(2.5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,4] < 0.05) {
            if (pvalue[i,4]>0.0001){
                a=round(pvalue[i,4],4)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,4],3)
                segments(1.1,ma5,4.9,ma5)
                text(2.75,ma52,"p=",cex=0.75)
                text(3,ma52,a,cex=0.75)
            }}

        if (pvalue[i,5] < 0.05) {
            if (pvalue[i,5]>0.0001){
                a=round(pvalue[i,5],4)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }else{
                a=signif(pvalue[i,5],3)
                segments(1.1,ma6,5.9,ma6)
                text(3.25,ma62,"p=",cex=0.75)
                text(3.5,ma62,a,cex=0.75)
            }}

        if (pvalue[i,6] < 0.05) {
            if (pvalue[i,6]>0.0001){
                a=round(pvalue[i,6],4)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }else{
                a=signif(pvalue[i,6],3)
                segments(1.1,ma11,6.9,ma11)
                text(3.75,ma112,"p=",cex=0.75)
                text(4,ma112,a,cex=0.75)
            }}

        if (pvalue[i,7] < 0.05) {
            if (pvalue[i,7]>0.0001){
                a=round(pvalue[i,7],4)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }else{
                a=signif(pvalue[i,7],3)
                segments(1.1,ma14,7.9,ma14)
                text(4.25,ma142,"p=",cex=0.75)
                text(4.5,ma142,a,cex=0.75)
            }}

        if (pvalue[i,8] < 0.05) {
            if (pvalue[i,8]>0.0001){
                a=round(pvalue[i,8],4)
                segments(1.1,ma18,8.9,ma18)
                text(4.75,ma182,"p=",cex=0.75)
                text(5,ma182,a,cex=0.75)
            }else{
                a=signif(pvalue[i,8],3)
                segments(1.1,ma18,8.9,ma18)
                text(4.75,ma182,"p=",cex=0.75)
                text(5,ma182,a,cex=0.75)
            }}

        if (pvalue[i,9] < 0.05) {
            if (pvalue[i,9]>0.0001){
                a=round(pvalue[i,9],4)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,9],3)
                segments(2.1,ma2,2.9,ma2)
                text(2.25,ma22,"p=",cex=0.75)
                text(2.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,10] < 0.05) {
            if (pvalue[i,10]>0.0001){
                a=round(pvalue[i,10],4)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,10],3)
                segments(2.1,ma7,3.9,ma7)
                text(2.75,ma72,"p=",cex=0.75)
                text(3,ma72,a,cex=0.75)
            }}

        if (pvalue[i,11] < 0.05) {
            if (pvalue[i,11]>0.0001){
                a=round(pvalue[i,11],4)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,11],3)
                segments(2.1,ma8,4.9,ma8)
                text(3.25,ma82,"p=",cex=0.75)
                text(3.5,ma82,a,cex=0.75)
            }}

        if (pvalue[i,12] < 0.05) {
            if (pvalue[i,12]>0.0001){
                a=round(pvalue[i,12],4)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }else{
                a=signif(pvalue[i,12],3)
                segments(2.1,ma9,5.9,ma9)
                text(3.75,ma92,"p=",cex=0.75)
                text(4,ma92,a,cex=0.75)
            }}

        if (pvalue[i,13] < 0.05) {
            if (pvalue[i,13]>0.0001){
                a=round(pvalue[i,13],4)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }else{
                a=signif(pvalue[i,13],3)
                segments(2.1,ma12,6.9,ma12)
                text(4.25,ma122,"p=",cex=0.75)
                text(4.5,ma122,a,cex=0.75)
            }}

        if (pvalue[i,14] < 0.05) {
            if (pvalue[i,14]>0.0001){
                a=round(pvalue[i,14],4)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }else{
                a=signif(pvalue[i,14],3)
                segments(2.1,ma15,7.9,ma15)
                text(4.75,ma152,"p=",cex=0.75)
                text(5,ma152,a,cex=0.75)
            }}

        if (pvalue[i,15] < 0.05) {
            if (pvalue[i,15]>0.0001){
                a=round(pvalue[i,15],4)
                segments(2.1,ma19,8.9,ma19)
                text(5.25,ma192,"p=",cex=0.75)
                text(5.5,ma192,a,cex=0.75)
            }else{
                a=signif(pvalue[i,15],3)
                segments(2.1,ma19,8.9,ma19)
                text(5.25,ma192,"p=",cex=0.75)
                text(5.5,ma192,a,cex=0.75)
            }}

        if (pvalue[i,16] < 0.05) {
            if (pvalue[i,16]>0.0001){
                a=round(pvalue[i,16],4)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,16],3)
                segments(3.1,ma2,3.9,ma2)
                text(3.25,ma22,"p=",cex=0.75)
                text(3.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,17] < 0.05) {
            if (pvalue[i,17]>0.0001){
                a=round(pvalue[i,17],4)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,17],3)
                segments(3.1,ma3,4.9,ma3)
                text(3.75,ma32,"p=",cex=0.75)
                text(4,ma32,a,cex=0.75)
            }}

        if (pvalue[i,18] < 0.05) {
            if (pvalue[i,18]>0.0001){
                a=round(pvalue[i,18],4)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,18],3)
                segments(3.1,ma10,5.9,ma10)
                text(4.25,ma102,"p=",cex=0.75)
                text(4.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,19] < 0.05) {
            if (pvalue[i,19]>0.0001){
                a=round(pvalue[i,19],4)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }else{
                a=signif(pvalue[i,19],3)
                segments(3.1,ma13,6.9,ma13)
                text(4.75,ma132,"p=",cex=0.75)
                text(5,ma132,a,cex=0.75)
            }}

        if (pvalue[i,20] < 0.05) {
            if (pvalue[i,20]>0.0001){
                a=round(pvalue[i,20],4)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }else{
                a=signif(pvalue[i,20],3)
                segments(3.1,ma16,7.9,ma16)
                text(5.25,ma162,"p=",cex=0.75)
                text(5.5,ma162,a,cex=0.75)
            }}

        if (pvalue[i,21] < 0.05) {
            if (pvalue[i,21]>0.0001){
                a=round(pvalue[i,21],4)
                segments(3.1,ma20,8.9,ma20)
                text(5.75,ma202,"p=",cex=0.75)
                text(6,ma202,a,cex=0.75)
            }else{
                a=signif(pvalue[i,21],3)
                segments(3.1,ma20,8.9,ma20)
                text(5.75,ma202,"p=",cex=0.75)
                text(6,ma202,a,cex=0.75)
            }}

        if (pvalue[i,22] < 0.05) {
            if (pvalue[i,22]>0.0001){
                a=round(pvalue[i,22],4)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,22],3)
                segments(4.1,ma2,4.9,ma2)
                text(4.25,ma22,"p=",cex=0.75)
                text(4.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,23] < 0.05) {
            if (pvalue[i,23]>0.0001){
                a=round(pvalue[i,23],4)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,23],3)
                segments(4.1,ma4,5.9,ma4)
                text(4.75,ma42,"p=",cex=0.75)
                text(5,ma42,a,cex=0.75)
            }}

        if (pvalue[i,24] < 0.05) {
            if (pvalue[i,24]>0.0001){
                a=round(pvalue[i,24],4)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,24],3)
                segments(4.1,ma7,6.9,ma7)
                text(5.25,ma72,"p=",cex=0.75)
                text(5.5,ma72,a,cex=0.75)
            }}

        if (pvalue[i,25] < 0.05) {
            if (pvalue[i,25]>0.0001){
                a=round(pvalue[i,25],4)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }else{
                a=signif(pvalue[i,25],3)
                segments(4.1,ma17,7.9,ma17)
                text(5.75,ma172,"p=",cex=0.75)
                text(6,ma172,a,cex=0.75)
            }}

        if (pvalue[i,26] < 0.05) {
            if (pvalue[i,26]>0.0001){
                a=round(pvalue[i,26],4)
                segments(4.1,ma21,8.9,ma21)
                text(6.25,ma212,"p=",cex=0.75)
                text(6.5,ma212,a,cex=0.75)
            }else{
                a=signif(pvalue[i,26],3)
                segments(4.1,ma21,8.9,ma21)
                text(6.25,ma212,"p=",cex=0.75)
                text(6.5,ma212,a,cex=0.75)
            }}

        if (pvalue[i,27] < 0.05) {
            if (pvalue[i,27]>0.0001){
                a=round(pvalue[i,27],4)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,27],3)
                segments(5.1,ma2,5.9,ma2)
                text(5.25,ma22,"p=",cex=0.75)
                text(5.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,28] < 0.05) {
            if (pvalue[i,28]>0.0001){
                a=round(pvalue[i,28],4)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }else{
                a=signif(pvalue[i,28],3)
                segments(5.1,ma3,6.9,ma3)
                text(5.75,ma32,"p=",cex=0.75)
                text(6,ma32,a,cex=0.75)
            }}

        if (pvalue[i,29] < 0.05) {
            if (pvalue[i,29]>0.0001){
                a=round(pvalue[i,29],4)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }else{
                a=signif(pvalue[i,29],3)
                segments(5.1,ma5,7.9,ma5)
                text(6.25,ma52,"p=",cex=0.75)
                text(6.5,ma52,a,cex=0.75)
            }}

        if (pvalue[i,30] < 0.05) {
            if (pvalue[i,30]>0.0001){
                a=round(pvalue[i,30],4)
                segments(5.1,ma8,8.9,ma8)
                text(6.75,ma82,"p=",cex=0.75)
                text(7,ma82,a,cex=0.75)
            }else{
                a=signif(pvalue[i,30],3)
                segments(5.1,ma8,8.9,ma8)
                text(6.75,ma82,"p=",cex=0.75)
                text(7,ma82,a,cex=0.75)
            }}

        if (pvalue[i,31] < 0.05) {
            if (pvalue[i,31]>0.0001){
                a=round(pvalue[i,31],4)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,31],3)
                segments(6.1,ma2,6.9,ma2)
                text(6.25,ma22,"p=",cex=0.75)
                text(6.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,32] < 0.05) {
            if (pvalue[i,32]>0.0001){
                a=round(pvalue[i,32],4)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }else{
                a=signif(pvalue[i,32],3)
                segments(6.1,ma4,7.9,ma4)
                text(6.75,ma42,"p=",cex=0.75)
                text(7,ma42,a,cex=0.75)
            }}

        if (pvalue[i,33] < 0.05) {
            if (pvalue[i,33]>0.0001){
                a=round(pvalue[i,33],4)
                segments(6.1,ma10,8.9,ma10)
                text(7.25,ma102,"p=",cex=0.75)
                text(7.5,ma102,a,cex=0.75)
            }else{
                a=signif(pvalue[i,33],3)
                segments(6.1,ma10,8.9,ma10)
                text(7.25,ma102,"p=",cex=0.75)
                text(7.5,ma102,a,cex=0.75)
            }}

        if (pvalue[i,34] < 0.05) {
            if (pvalue[i,34]>0.0001){
                a=round(pvalue[i,34],4)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,34],3)
                segments(7.1,ma2,7.9,ma2)
                text(7.25,ma22,"p=",cex=0.75)
                text(7.5,ma22,a,cex=0.75)
            }}

        if (pvalue[i,35] < 0.05) {
            if (pvalue[i,35]>0.0001){
                a=round(pvalue[i,35],4)
                segments(7.1,ma7,8.9,ma7)
                text(7.75,ma72,"p=",cex=0.75)
                text(8,ma72,a,cex=0.75)
            }else{
                a=signif(pvalue[i,35],3)
                segments(7.1,ma7,8.9,ma7)
                text(7.75,ma72,"p=",cex=0.75)
                text(8,ma72,a,cex=0.75)
            }}

        if (pvalue[i,36] < 0.05) {
            if (pvalue[i,36]>0.0001){
                a=round(pvalue[i,36],4)
                segments(8.1,ma2,8.9,ma2)
                text(8.25,ma22,"p=",cex=0.75)
                text(8.5,ma22,a,cex=0.75)
            }else{
                a=signif(pvalue[i,36],3)
                segments(8.1,ma2,8.9,ma2)
                text(8.25,ma22,"p=",cex=0.75)
                text(8.5,ma22,a,cex=0.75)
            }}
    }
    dev.off()

}

setwd(dirout)
}
