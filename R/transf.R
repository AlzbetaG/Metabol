#' Transformation
#'
#' Transform data according to chosen function.
#' @param data Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.
#' @param name A character string or expression indicating a name of data set. It occurs in names of every output.
#' @param tsf Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.
#' @details Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.
#' @details Up to twenty different groups can be distinguished in data (including QCs).
#' @return Excel file with data transformed according to chosen function.
#' @importFrom robCompositions cenLR
#' @references Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.
#' @examples data=metabol
#' name="Metabolomics"    #name of the project
#' groupnames=c("Con","Pat","QC")
#' transf(data,name)
#' @export
transf=function(data,name,tsf="clr"){

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

    write.table(dataM, file = paste("12_",tsf,"_",name,".csv",sep=""), sep=";",col.names =NA)     # Final table

}


