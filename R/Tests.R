#' Standard tests
#'
#' Makes Shapiro-Wilk test, ANOVA, Kruskal-Wallis test, t-test for every variable and Wilcoxon test for every variable.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param groupnames A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @param pair logical. If TRUE then the paired tests are done. Default is FALSE.
#' @param QCs logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @return Excel file with p-values from Shapiro-Wilk test, ANOVA, Kruskal-Wallis test, t-test for every variable and Wilcoxon test for every variable and their summary.
#' @importFrom stats kruskal.test oneway.test shapiro.test t.test wilcox.test
#' @importFrom robCompositions cenLR
#' @import openxlsx
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' Tests(data,name,groupnames)
#' @export
Tests=function(data,name,groupnames,tsf="clr",pair=FALSE,QCs=FALSE){

  ######################################################################################
  dirout = paste(getwd(),"/",sep = "")
  #dir.create(dirout, recursive=TRUE, showWarnings = FALSE)

  wb <- createWorkbook()
  sheet1  <- addWorksheet(wb, sheetName = "Summary")
  sheet2  <- addWorksheet(wb, sheetName = "Normality")
  sheet3  <- addWorksheet(wb, sheetName = "ANOVA")
  sheet4  <- addWorksheet(wb, sheetName = "Kruskal")

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
  ######################################################################################
  count=length(groupnames)
  groupss=matrix(rep(NA,count*2),ncol=2)
  colnames(groupss)=c("min","max")
  rownames(groupss)=groupnames
  groups=NULL
  for (i in 1:count){
      Gr=grep(groupnames[i],rownames(dataM))
      groupss[i,1]=min(Gr)
      groupss[i,2]=max(Gr)
      gr=rep(i,length(Gr))
      groups=c(groups,gr)
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

######################################################################################
force(groupss)

names<-colnames(dataM)

######################################################################################
#Normality test:

shap=matrix(rep(0,ncol(dataM)),nrow=ncol(dataM))
rownames(shap)=names

for(i in 1:ncol(dataM)){
    shap[i]=shapiro.test(dataM[,i])$p.value
}

colnames(shap)="P_value"

writeData(wb,sheet2,shap,colNames = TRUE, rowNames = TRUE)

pocet=function(x){length(x[x==TRUE])}
sumnorm=pocet(shap>0.05)

#################################################################################################################
# ANOVA:

data=as.matrix(cbind(dataM,groups))

anova=rep(0,ncol(data)-1)
for(i in 1:(ncol(data)-1)){
    anova[i]=oneway.test(data[,i] ~ groups,data,var.equal=TRUE)$p.value
}

anova=as.matrix(anova)
rownames(anova)=names

colnames(anova)="P_value"

writeData(wb,sheet3,anova,colNames = TRUE, rowNames = TRUE)

sumanova=pocet(anova<0.05)

#################################################################################################################
#Kruskal-Wallis:
data=as.matrix(cbind(dataM,groups))

KW=rep(0,ncol(data)-1)
for(i in 1:(ncol(data)-1)){
    KW[i]=kruskal.test(data[,i] ~ groups,data)$p.value
}

KW=as.matrix(KW)
rownames(KW)=names

colnames(KW)="P_value"

writeData(wb,sheet4,KW,colNames = TRUE, rowNames = TRUE)

sumKW=pocet(KW<0.05)

####################################################################################################

if (pair == TRUE) {

    sheet5  <- addWorksheet(wb, sheetName = "Wilcoxon")
    sheet6  <- addWorksheet(wb, sheetName = "T-test")

    file00=paste(dirout, "Tests_pair_",name,".xlsx",sep="")

    #Wilcoxon test:

    pvaluew=NULL
    for(k in 1:(ncol(dataM))){
        pvaluew2=NULL
        for(i in 1:(nrow(groupss)-1)){
            A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
            for(j in (i+1):nrow(groupss)){
                B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                c=wilcox.test(A,B,alternative = "two.sided",paired =TRUE)$p.value
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

    writeData(wb,sheet5,pvaluew,colNames = TRUE, rowNames = TRUE)

    sumpw=apply(pvaluew<0.05,2,pocet)

    #t-test:

    pvaluet=NULL
    for(k in 1:(ncol(dataM))){
        pvaluet2=NULL
        for(i in 1:(nrow(groupss)-1)){
            A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
            for(j in (i+1):nrow(groupss)){
                B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                c=t.test(A,B,alternative = "two.sided",paired =TRUE)$p.value
                pvaluet2=cbind(pvaluet2,c)
            }
        }
        pvaluet=rbind(pvaluet,pvaluet2)
    }

    nazvy=NULL
    for(i in 1:(nrow(groupss)-1)){
        for(j in (i+1):nrow(groupss)){
            col=paste(groupnames[i],"_",groupnames[j],sep="")
            nazvy=c(nazvy,col)
        }
    }

    colnames(pvaluet)=nazvy
    row.names(pvaluet)=names

    writeData(wb,sheet6,pvaluet,colNames = TRUE, rowNames = TRUE)

    sumpv=apply(pvaluet<0.05,2,pocet)

    ################################################################################################################################
    #summaries:

    sumtabA2=cbind(rep(NA,2),rbind(sumpw,sumpv))

    sumtabA1=cbind(rbind(sumnorm,sumanova,sumKW),rbind(rbind(rep(NA,(ncol(sumtabA2)-1))),
                                                       rbind(rep(NA,(ncol(sumtabA2)-1))),rbind(rep(NA,(ncol(sumtabA2)-1)))))

    sumtabA=rbind(sumtabA1,sumtabA2)

    rownames(sumtabA)=c("Normality","ANOVA","Kruskal","Wilcoxon","T-test")

    colnames(sumtabA)[1]="All"

    sumtabB2=cbind(rep(NA,2),rbind(round(100*sumpw/ncol(dataM),2),round(100*sumpv/ncol(dataM),2)))

    sumtabB1=cbind(rbind(round(100*sumnorm/ncol(dataM),2),round(100*sumanova/ncol(dataM),2),round(100*sumKW/ncol(dataM),2)),
                   rbind(rbind(rep(NA,(ncol(sumtabB2)-1))),rbind(rep(NA,(ncol(sumtabB2)-1))),rbind(rep(NA,(ncol(sumtabB2)-1)))))

    sumtabB=rbind(sumtabB1,sumtabB2)

    rownames(sumtabB)=c("Normality (%)","ANOVA (%)","Kruskal (%)","Wilcoxon (%)","T-test (%)")

    sumtab=rbind(sumtabA,sumtabB)

    saveWorkbook(wb, file = file00,overwrite = TRUE)

}

####################################################################################################
else {
    sheet5  <- addWorksheet(wb, sheetName = "Wilcoxon")
    sheet6  <- addWorksheet(wb, sheetName = "T-test")

    file00=paste(dirout, "Tests_",name,".xlsx",sep="")

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

    writeData(wb,sheet5,pvaluew,colNames = TRUE, rowNames = TRUE)

    sumpw=apply(pvaluew<0.05,2,pocet)

    #t-test:

    pvaluet=NULL
    for(k in 1:(ncol(dataM))){
        pvaluet2=NULL
        for(i in 1:(nrow(groupss)-1)){
            A=as.matrix(dataM[groupss[i,1]:groupss[i,2],k])
            for(j in (i+1):nrow(groupss)){
                B=as.matrix(dataM[groupss[j,1]:groupss[j,2],k])
                c=t.test(A,B,alternative = "two.sided")$p.value
                pvaluet2=cbind(pvaluet2,c)
            }
        }
        pvaluet=rbind(pvaluet,pvaluet2)
    }

    nazvy=NULL
    for(i in 1:(nrow(groupss)-1)){
        for(j in (i+1):nrow(groupss)){
            col=paste(groupnames[i],"_",groupnames[j],sep="")
            nazvy=c(nazvy,col)
        }
    }

    colnames(pvaluet)=nazvy
    row.names(pvaluet)=names

    writeData(wb,sheet6,pvaluet,colNames = TRUE, rowNames = TRUE)

    sumpv=apply(pvaluet<0.05,2,pocet)

    ################################################################################################################################
    #summaries:

    sumtabA2=cbind(rep(NA,2),rbind(sumpw,sumpv))

    sumtabA1=cbind(rbind(sumnorm,sumanova,sumKW),rbind(rbind(rep(NA,(ncol(sumtabA2)-1))),
                                                       rbind(rep(NA,(ncol(sumtabA2)-1))),rbind(rep(NA,(ncol(sumtabA2)-1)))))

    sumtabA=rbind(sumtabA1,sumtabA2)

    rownames(sumtabA)=c("Normality","ANOVA","Kruskal","Wilcoxon","T-test")

    colnames(sumtabA)[1]="All"

    sumtabB2=cbind(rep(NA,2),rbind(round(100*sumpw/ncol(dataM),2),round(100*sumpv/ncol(dataM),2)))

    sumtabB1=cbind(rbind(round(100*sumnorm/ncol(dataM),2),round(100*sumanova/ncol(dataM),2),round(100*sumKW/ncol(dataM),2)),
                   rbind(rbind(rep(NA,(ncol(sumtabB2)-1))),rbind(rep(NA,(ncol(sumtabB2)-1))),rbind(rep(NA,(ncol(sumtabB2)-1)))))

    sumtabB=rbind(sumtabB1,sumtabB2)

    rownames(sumtabB)=c("Normality (%)","ANOVA (%)","Kruskal (%)","Wilcoxon (%)","T-test (%)")

    sumtab=rbind(sumtabA,sumtabB)

    writeData(wb,sheet1,sumtab,colNames = TRUE, rowNames = TRUE)

    saveWorkbook(wb, file = file00,overwrite = TRUE)
}
}

