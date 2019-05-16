#' Quality control samples (QCs) checking
#'
#' Quality control samples (QCs) are checked to data irregularities. It is used for data from untargeted metabolomic analysis.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @details Values of QCs are evaluated and questionable values for particular variables are denoted. There are two steps of evaluation: 1. QCs with completely higher values than the maximum of data, 2. QCs higher than majority of data.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return Boxplots of QCs and the other data groups.
#' @return Excel file with the list of questionable variables from two steps of evaluation.
#' @import openxlsx
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' bigQC(data,name,groupnames)
#' @export
bigQC=function(data,name,groupnames){

    ################################################################################################################################
    #data=as.matrix(data)

    ##########################################################################################################################
    basecolor=c("blue","magenta","forestgreen","darkorange","deepskyblue","mediumaquamarine","lightslateblue","saddlebrown",
                "gray40","darkslateblue","firebrick","darkcyan","darkmagenta", "deeppink1","limegreen","gold2","bisque2",
                "lightcyan3","red","darkolivegreen3")           # Basic colours from: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    basemarks=c(15,17,18,8,11,2,0,16,5,6,4,10,3,7,9,12)

    groupnames=groupnames
    #groupnames=unique(gsub("[[:digit:]]","",rownames(data)))
    count=length(groupnames)
    groups=NULL
    marks=NULL
    color=NULL
    for (i in 1:count){
        Gr=grep(groupnames[i],rownames(data))
        gr=rep(i,length(Gr))
        groups=c(groups,gr)
        zn=rep(basemarks[i],length(Gr))
        marks=c(marks,zn)
        cl=rep(basecolor[i],length(Gr))
        color=c(color,cl)
    }
################################################################################################################################
# denoting of QCs

QCi=grep("QC",rownames(data))
dataQC=data[QCi,]

################################################################################################################################
# rule 1 - comparison of maximum of samples and minimum of QCs

rule1=matrix(rep(NA,ncol(data)),ncol=1)
for(i in 1:ncol(data)){
    maxs=max(data[-QCi,i])
    b=boxplot(data[QCi,i] ~ groups[QCi], names=groupnames[1],main=colnames(data)[i],notch=FALSE,plot=FALSE)
    minQC=b$stats[1,1]
    if (maxs<minQC){
        rule1[i,1]=1
    } else {
        rule1[i,1]=0
    }
}

rownames(rule1)=colnames(data)

#head(rule1)

idxrule1 = which(rule1 == 1)

if (length(idxrule1)!=0){
    data2=data[,-idxrule1]
    dataout=matrix(rep(0,nrow(data)*length(idxrule1)),nrow=nrow(data))
    rownames(dataout)=rownames(data)
    for (k in 1:length(idxrule1)){
        dataout[,k]=data[,idxrule1[k]]
        colnames(dataout)=colnames(data)[idxrule1]
    }

    write.xlsx(dataout,file = paste("Box_out_rule_1_",name,".xlsx",sep=""),sheetName="Out",
           col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)

    labels=rownames(dataout)

    pdf(paste("Box_out_rule_1_",name,".pdf",sep=""))
    for(i in 1:ncol(dataout)){
        boxplot(dataout[,i] ~ groups, names=groupnames,main=colnames(dataout)[i],notch=TRUE,outpch = NA)
        text(groups,dataout[,i],label=labels,col="red",cex=0.5)
        }
    dev.off()
    }else{
        data2 = data
        print("No questionable QCs in rule 1.")
        }
#unique(gsub("[[:digit:]]","",rownames(dataSet)))
################################################################################################################################
# rule 2 - QCs higher than majority of data (some samples are higher than QCs)

rule2=matrix(rep(NA,ncol(data2)*1),ncol=1)
for(i in 1:ncol(data2)){
    b=boxplot(data2[,i] ~ groups, names=groupnames,main=colnames(data2)[i],notch=FALSE,plot=FALSE)
    qc=grep("QC",groupnames)
    cAQC=b$conf[1,qc]
    cBs=max(b$stats[4,-qc])

    if (cAQC>cBs){
        rule2[i,1]=1
        } else {
        rule2[i,1]=0
        }
}

rownames(rule2)=colnames(data2)

#head(rule2)

idxrule2 = which(apply(rule2,1,sum) == 1)

if (length(idxrule2)!=0){
    data3=data2[,-idxrule2]
    dataout2=matrix(rep(0,nrow(data2)*length(idxrule2)),nrow=nrow(data2))
    rownames(dataout2)=rownames(data)
    for (k in 1:length(idxrule2)){
        dataout2[,k]=data2[,idxrule2[k]]
        colnames(dataout2)=colnames(data2)[idxrule2]
    }

    write.xlsx(dataout2,file = paste("Box_out_rule_2_",name,".xlsx",sep=""),sheetName="Out",
           col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)

    labels=rownames(dataout2)

    pdf(paste("Box_out_rule_2_",name,".pdf",sep=""))
    for(i in 1:ncol(dataout2)){
        b=boxplot(dataout2[,i] ~ groups, names=groupnames,main=colnames(dataout2)[i],notch=TRUE,outpch = NA)
        text(groups,dataout2[,i],label=labels,col="red",cex=0.5)
    }
    dev.off()

    pdf(paste("Box_rest_",name,".pdf",sep=""))
    for(i in 1:ncol(data3)){
        b=boxplot(data3[,i] ~ groups, names=groupnames,main=colnames(data3)[i],notch=TRUE,outpch = NA)
        stripchart(data3[,i] ~ groups, vertical = TRUE, method = "jitter",pch = unique(marks), col = unique(color), add = TRUE)
        }
    dev.off()
    }else{
        data3 = data2
        print("No questionable QCs in rule 2.")
        }
}
