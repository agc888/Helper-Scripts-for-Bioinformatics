library(ctc)
library(heatmap.plus)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(utils)
library(stats)
library(grDevices)

#' Estimate PAM50 Subtype and Risk Scores
#'
#' This function performs PAM50 subtype prediction and calculates genomic and risk scores 
#' based on a given Seurat object. It supports multiple input configurations and provides 
#' options for visualizations and data export.
#'
#' @param seurat_obj A Seurat object containing the expression data for the analysis.
#' @param group.by A character string specifying the metadata column to group cells by (e.g., cluster labels).
#' @param n An integer specifying the number of groups to pool for the analysis. Default is 3.
#' @param assay The assay to use for retrieving expression data. Defaults to the active assay of the Seurat object.
#' @param slot The slot from the assay to retrieve data from. Default is "counts".
#' @param save_plots_to A character string specifying the directory to save generated plots. If NULL, no plots are saved.
#' @param collapseMethod A character string specifying the method to collapse data. Options are "mean" (default) or "iqr".
#' @param tumour_mass_column A character string specifying the column name in the Seurat object's metadata 
#'        containing tumour mass information (if available). If NULL, tumour size data is not included.
#' @param calibrationParameters A numeric or NA value specifying the column of the calibration file to use. 
#'        NA will center within the test set; -1 will apply no calibration.
#' @param stdArray A logical value indicating whether to standardize the array data for visualization. Default is TRUE.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#'
#' @return A list containing the PAM50 subtype predictions, risk scores, and related data.
#'
#' # REQUIRED PACKAGES: ctc, heatmap.plus, Seurat, ggplot2, dplyr, stringr, SingleCellExperiment, scater, utils, stats, grDevices
#'
#'@export
EstimatePAM50 <- function(seurat_obj, group.by, n = 3, assay = DefaultAssay(seurat_obj),slot = "counts", save_plots_to = NULL,collapseMethod = "mean", tumour_mass_column = NULL, calibrationParameters = NA,stdArray = TRUE, verbose = TRUE){
    
    pooled_data <- run_pooling(data.filt = seurat_obj,
                               idents = group.by, 
                               n = n, 
                               assay = assay, 
                               slot = slot, 
                               verbose = verbose)

    
    pam50 <- read.table("https://raw.githubusercontent.com/agc888/Helper-Scripts-for-Bioinformatics/refs/heads/main/R%20Scripts/PAM50_R/pam50_annotation.txt", header = TRUE, sep = "\t", fill = TRUE)

    mtx <- pooled_data$mtx
    mtx <- mtx[rownames(mtx) %in% pam50$pcrID,]


    short<-group.by # short name that will be used for output files
    
    calibrationParameters<- calibrationParameters 	#the column of the "mediansPerDataset.txt" file to use for calibration; 
    														#NA will force centering within the test set & -1 will not do any 
    														#adjustment (when adjustment performed by used)
    if (is.null(tumour_mass_column)){
        hasClinical<-FALSE 	#may include tumor size as second row, with 'T' as the gene name, 
    										#and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
    										#set this variable to FALSE if tumor size is not available
    } else {

        tumour_mass <- unlist(lapply(colnames(mtx), function(x){
            cells <- names(pooled_data$ids[pooled_data$ids == x])
            tumour_size_df <- seurat_obj@meta.data[rownames(seurat_obj@meta.data) %in% cells,]
            mean(tumour_size_df[[tumour_mass_column]], na.rm = TRUE)
        }))
        
        new_row <- matrix(tumour_mass, nrow = 1, ncol = ncol(mtx))  # Replace 0 with desired values
        # Set row name for the new row
        rownames(new_row) <- "T"
        
        # Combine the new row with the existing matrix
        mtx <- rbind(new_row, mtx)
        
        hasClinical<-TRUE
    }
    
    collapseMethod<- collapseMethod # can be mean or iqr (probe with max iqr is selected)
    											# typically, mean is preferred for long oligo and
    											# iqr is preferred for short oligo platforms


    ####
    # Get Additional Functions
    ####
    
    source("https://raw.githubusercontent.com/agc888/Helper-Scripts-for-Bioinformatics/refs/heads/main/R%20Scripts/PAM50_R/subtypePrediction_functions.R")
    
    ## NOTE: ALL CODE BELOW HAS BEEN ADAPTED FROM 'https://genome-publications.bioinf.unc.edu/PAM50/' specifically from the 'PAM50_R/subtypePrediction_distributed.R' file. 


    # this code & protocol has been validated externally by MDACC
    #1. Download normalized expression data
    #2. Log transform expression estimates
    #3. Optionally median center or appropriately adjust each probeset, which carries the population assumption
    #4. Map probesets to Entrez gene names
    #5. Format as tab delimited text, with first row of sample names and first column of gene names
    #---- This software then provides the following steps
    #5. For probesets that map to identical Entrez gene names, select the one with highest IQR (for Affy, select mean for Agilent)
    #6. Extract the 50 genes of interest (in pam50_centroids.txt)
    #7. Calculate Spearman's rank correlation between each sample and each subtype centroid (in pam50_centroids.txt)
    #8. Assign the class of the most highly correlated centroid to each sample
    
    calibrationFile<- "https://raw.githubusercontent.com/agc888/Helper-Scripts-for-Bioinformatics/refs/heads/main/R%20Scripts/PAM50_R/mediansPerDataset_v2.txt"
    

    trainCentroids<- "https://raw.githubusercontent.com/agc888/Helper-Scripts-for-Bioinformatics/refs/heads/main/R%20Scripts/PAM50_R/pam50_centroids.txt"
    trainFile<- "https://raw.githubusercontent.com/agc888/Helper-Scripts-for-Bioinformatics/refs/heads/main/R%20Scripts/PAM50_R/220arrays_nonUBCcommon+12normal_50g.txt"
    proliferationGenes<-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
    stdArray<-stdArray # just for visualization, and only set to F if many missing genes
    
    #NOTE: predFiles has now been changed to read directly from the mtx file created by pooling, rather then writing and reading a temp .txt file
    #predFiles<- paste(inputDir,inputFile,sep="/")
    
    verbose_message(message_text = "Starting PAM50 estimation analysis ... ", verbose = verbose)
    
    ###
    # some constants
    ###
    
    # for subtype only model
    glthreshold<- -0.15
    ghthreshold<-  0.1
    
    # for subtype + proliferation model
    gplthreshold<- -0.25
    gphthreshold<-  0.1
    
    # for combined model
    clthreshold<- -0.1
    chthreshold<-  0.2
    
    # for combined + proliferation model
    cplthreshold<- -0.2
    cphthreshold<-  0.2
    
    # begin analyses

    verbose_message(message_text = "Loading in training data and other pre-established constants ... ", verbose = verbose)

    
    # only need train data for visualizations
    x<-readarray(trainFile,hr=2)
    x$xd<-standardize(medianCtr(x$xd))

    
    # load the published centroids for classifcation
    pamout.centroids<-read.table(trainCentroids,sep="\t",header=T,row.names=1)
    
    if (is.null(save_plots_to)){
        warning("`save_plot_to = NULL`! general plots generated by base function will not be generated or saved. To save plots, please specify a filepath to `save_plots_to = ...`")
    } else {
        verbose_message(message_text = paste0("Output files and plots will be saved here: ", save_plots_to), verbose = verbose)

    	pdfname1<-paste(save_plots_to,paste("predictionScores_pam50RankCorrelation_1_",short,".pdf",sep=""),sep="/")
    	pdfname2<-paste(save_plots_to,paste("predictionScores_pam50RankCorrelation_2_",short,".pdf",sep=""),sep="/")
    	clustername<-paste(save_plots_to,paste(short,"_PAM50_normalized_heatmap",sep=""),sep="/")
    	outFile<- paste(save_plots_to,paste(short,"_pam50scores.txt",sep=""),sep="/")
    }
    
    
   
    
    # read in the data file
    if(hasClinical){
        xhr=2
    }else{
        xhr=1
    }

    mtx <- as.data.frame(mtx)
    mtx$Gene <- rownames(mtx)
    mtx <- rbind(colnames(mtx), mtx)  # Insert column names as the first row
    mtx <- mtx[c("Gene", colnames(mtx)[1:(length(colnames(mtx))-1)])]
    rownames(mtx) <- NULL
    mtx$Gene[1] <- ""
    colnames(mtx) <- paste0("V", 1:ncol(mtx))
    
    y<-mtxToArray(mtx,hr=xhr,method=collapseMethod,impute=F)
    
    # normalization
    if(is.na(calibrationParameters)){
        y$xd<-medianCtr(y$xd)
    }else{
        if(calibrationParameters != -1){
            medians<-readarray(calibrationFile,hr=1)
            verbose_message(message_text = paste("calibration to:",dimnames(medians$xd)[[2]][calibrationParameters]), verbose = verbose)
            tm<-overlapSets(medians$xd,y$xd)
            y$xd<-(tm$y-tm$x[,calibrationParameters])
            #y$xd<-(tm$y-tm$x[,calibrationParameters])/tm$x[,15]
        }
    }
    	
    num.missing<- NA

    if(stdArray){
        y$xd<-standardize(y$xd)
    }


    erScore<-as.vector(t(y$xd["ESR1",]))


    her2Score<-as.vector(t(y$xd["ERBB2",]))



    # assign the subtype scores and calculate the proliferation score
    this.proliferationGenes<-dimnames(y$xd)[[1]] %in% proliferationGenes

    prolifScore<-apply(y$xd[this.proliferationGenes,],2,mean,na.rm=T)

    out<-sspPredict(pamout.centroids,classes="",y$xd,std=F,distm="spearman",centroids=T)
    out$distances<- -1*out$distances
    
    call.conf<-c()
    for(j in 1:length(out$predictions)){
        call.conf[j]<- 1-cor.test(out$testData[,j],out$centroids[,which(colnames(pamout.centroids)==out$predictions[j])],method="spearman")$p.value
    }
    call.conf<-round(call.conf,2)
    
    # calculate the risk scores
    genomic <- 0.04210193*out$distances[,1] + 0.12466938*out$distances[,2] + -0.35235561*out$distances[,3] + 0.14213283*out$distances[,4]
    genomicWprolif <- -0.0009299747*out$distances[,1] + 0.0692289192*out$distances[,2] + -0.0951505484*out$distances[,3] +  0.0493487685*out$distances[,4] + 0.3385116381*prolifScore
    if(hasClinical){
        xT<-as.numeric(as.vector(y$classes$T))
        combined <- 0.0442770*out$distances[,1] + 0.1170297*out$distances[,2] + -0.2608388*out$distances[,3] + 0.1055908*out$distances[,4] + 0.1813751*xT
        combinedWprolif <- -0.009383416*out$distances[,1] +  0.073725503*out$distances[,2] + -0.090436516*out$distances[,3] + 0.053013865*out$distances[,4] + 0.131605960*xT + 0.327259375*prolifScore
    }
    
    # threshold the risk score
    griskgroups<-genomic
    griskgroups[genomic>ghthreshold]<-"high"
    griskgroups[genomic>glthreshold & genomic<ghthreshold]<-"med"
    griskgroups[genomic<glthreshold]<-"low"
    gpriskgroups<-genomicWprolif
    gpriskgroups[genomicWprolif>gphthreshold]<-"high"
    gpriskgroups[genomicWprolif>gplthreshold & genomicWprolif<gphthreshold]<-"med"
    gpriskgroups[genomicWprolif<gplthreshold]<-"low"
    
    genomic<- 100* (genomic + 0.35 ) / 0.85
    genomicWprolif<- 100* (genomicWprolif + 0.35 ) / 0.85
    
    # write output files
    if(hasClinical){
        criskgroups<-combined
        criskgroups[combined>chthreshold]<-"high"
        criskgroups[combined>clthreshold & combined<chthreshold]<-"med"
        criskgroups[combined<clthreshold]<-"low"
        cpriskgroups<-combinedWprolif
        cpriskgroups[combinedWprolif>cphthreshold]<-"high"
        cpriskgroups[combinedWprolif>cplthreshold & combinedWprolif<cphthreshold]<-"med"
        cpriskgroups[combinedWprolif<cplthreshold]<-"low"
    
        combined<- 100* (combined + 0.35 ) / 0.85
        combinedWprolif<- 100* (combinedWprolif + 0.35 ) / 0.85
        
        outtable<-cbind(out$distances, out$predictions, call.conf, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, combined, criskgroups, combinedWprolif, cpriskgroups, erScore, her2Score)
        dimnames(outtable)[[2]]<-c("Basal","Her2","LumA","LumB","Normal","Call","Confidence",
                                                                "ROR-S (Subtype Only)","ROR-S Group (Subtype Only)","Proliferation Score", 
                                                                "ROR-P (Subtype + Proliferation)","ROR-P Group (Subtype + Proliferation)",
                                                                "ROR-C (Subtype + Clinical)","ROR-C Group (Subtype + Clinical)",
                                                                "ROR-PC (Subtype + Clinical + Proliferation)","ROR-PC Group (Subtype + Clinical + Proliferation)",
                                                                "ER","Her2")
    }else{
        outtable<-cbind(out$distances, out$predictions, call.conf, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, erScore, her2Score)
        dimnames(outtable)[[2]]<-c("Basal","Her2","LumA","LumB","Normal","Call","Confidence",
                                                                "ROR-S (Subtype Only)","ROR-S Group (Subtype Only)","Proliferation Score", 
                                                                "ROR-P (Subtype + Proliferation)","ROR-P Group (Subtype + Proliferation)",
                                                                "ER","Her2")
    }

    verbose_message(message_text = paste("ER range:",quantile(erScore,.9,na.rm=T)-quantile(erScore,.1,na.rm=T)), verbose = verbose)
    
    if (!is.null(save_plots_to)){
        write.table(outtable,outFile,sep="\t",col.names=NA)
        
        # make some plots for evaluation
        verbose_message(message_text = paste("ER range:",quantile(erScore,.9,na.rm=T)-quantile(erScore,.1,na.rm=T)), verbose = verbose)
        
        subtypeColors<-out$predictions
        subtypeColors[subtypeColors=="Basal"]<-"red"
        subtypeColors[subtypeColors=="Her2"]<-"hotpink"
        subtypeColors[subtypeColors=="LumA"]<-"darkblue"
        subtypeColors[subtypeColors=="LumB"]<-"skyblue"
        subtypeColors[subtypeColors=="Normal"]<-"green"
        conf.colors<-call.conf
        conf.colors[call.conf>=0.95]<-"black"
        conf.colors[call.conf<0.95]<-"red"
    
 
        pdf(paste(clustername,".pdf",sep=""))
        myHeatmap(out$testData,cbind(subtypeColors,conf.colors),file=paste(clustername,".cdt",sep=""),rowNames=rownames(out$testData))
        dev.off()
        
        pdf(pdfname1,height=10,width=12)
        pars<-par(no.readonly=T)
        myplot(out,short,prolifScore)
        dev.off()
    
        tm<-overlapSets(x$xd,y$xd)
        tm$x<-tm$x[,!is.na(x$classes$subtype)]
        tm<-cbind(tm$x,	impute.knn(as.matrix(tm$y))$data)
        classes<-matrix(nrow=4,ncol=dim(tm)[2])
        nTrainSamples<-length(x$classes$subtype[!is.na(x$classes$subtype)])
    
        classes[1,]<-c(rep("train",nTrainSamples),rep(short,dim(y$xd)[2]))
        classes[2,]<-c(x$classes$subtype[!is.na(x$classes$subtype)],rep(NA,dim(tm)[2]-dim(x$xd[,!is.na(x$classes$subtype)])[2]))
        classes[3,]<-c(rep(NA,dim(tm)[2]-length(out$predictions)),out$predictions)
        
        pdf(pdfname2,height=6,width=12)
        tm<-scale(tm,center=F)
        par(mfrow=c(1,3))
        pcaEA(tm,classes[1,],mainStr="Traing and Test sets",showNames=F,showClasses=F)
        pcaEA(tm[,!is.na(as.vector(t(classes[2,])))],classes[2,!is.na(as.vector(t(classes[2,])))],mainStr="Training cases",showNames=F,showClasses=F,groupColors=c("red","hotpink","darkblue","skyblue","green"))
        pcaEA(tm[,!is.na(as.vector(t(classes[3,])))],classes[3,!is.na(as.vector(t(classes[3,])))],mainStr="Test cases",showNames=F,showClasses=F,groupColors=c("red","hotpink","darkblue","skyblue","green"))
        par(pars)
        dev.off()
    }

    verbose_message(message_text = "Analysis Finished!", verbose = verbose)
    
    return(outtable)

}

verbose_message <- function(message_text, verbose) {
  if (verbose) {
    message(message_text)
  }
}

#' Runs pooling of a merged Seurat Dataset to generate pseudo-replicates for each sample
#'       - This function is used by run_edgeR_annotations()
#'
#' @param data.filt A Seurat Object containing count values for pooling.
#' @param idents A character string defining the idents column to pool the data against.
#' @param n An integer defining the amount of pseudo-replicates to generate for each sample (default = 3).
#' @param assay Character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param slot Character string defining the assay storage slot to pull the relative mz intensity values from (default = "counts").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A SinglCellExpereiment object which contains pooled (n)-pseudo-replicate counts data based on the Seurat Object input
#' @export
#'
#' @examples
#' # run_pooling <- list(seuratObj, idents = "sample", n = 3, assay = "Spatial", slot = "counts")
run_pooling <- function(data.filt, idents, n, assay, slot, verbose = TRUE) {

  cell_metadata <- data.filt@meta.data
  samples <- unique(cell_metadata[[idents]])

  verbose_message(message_text = paste0("Pooling one sample into ", n ," replicates..."), verbose = verbose)

  nrg <- n
  
  for(i in c(1:length(samples))){
    set.seed(i)
    wo<-which(cell_metadata[[idents]]== samples[i])
    cell_metadata[wo,'orig.ident2']<-paste(samples[i],sample(c(1:n),length(wo)
                                                             ,replace=T,prob=rep(1/nrg,nrg)),sep='_')
  }
  gene_data <- row.names(data.filt)
  filtered.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.filt[[assay]][slot]),
                                       colData = cell_metadata)


  tempf=strsplit(filtered.sce@colData[["orig.ident2"]],'_')
  pid=NULL
  for(i in 1:length(tempf)){
    pidone=tempf[[i]]
    if(length(pidone)!=3){
      pidone=c(pidone[1],'yes',pidone[2])
    }
    pid=rbind(pid,pidone)
  }

  filtered.sce@colData$type=pid[,2]

  summed <- scater::aggregateAcrossCells(filtered.sce,
                                 id=SingleCellExperiment::colData(filtered.sce)[,'orig.ident2'])


  ids <- SingleCellExperiment::colData(filtered.sce)[,'orig.ident2']
  names(ids) <- rownames(cell_metadata)
    
  return(list("mtx" = SingleCellExperiment::counts(summed),
             "ids" = ids))
}


#' @title Overlap two datasets
#'
#' @description
#' Formatting function to read arrays and format for use in the claudinLow classifier.
#'
#' @usage
#' readArray(dataFile,designFile=NA,hr=1,impute=TRUE,method="mean")
#'
#' @param x matrix to convert
#' @param designFile Design of file.
#' @param hr Header rows as Present (2) or Absent (1).
#' @param impute whether data will be imputed or not.
#' @param method Default method is "mean".
#'
#' @return
#' A list
#'
#' @references
#' citation("claudinLow")
#'
#' @seealso
#' [genefu::claudinLow]
#'
#' @md
#' @importFrom impute impute.knn
#' @export
mtxToArray <- function(x, designFile=NA, hr=1, impute=FALSE, method="mean")
{

  headerRows <- hr

  if(headerRows==1){
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]
    classes<-NULL
    ids<-x[,1]
    xd<-x[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)
  }else{
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]

    classes<-x[1:(headerRows-1),]
    dimnames(classes)[[1]]<-classes[,1]
    classes<-classes[,-1]
    classes[classes==""]<-NA
    classes<-t(classes)
    rownames(classes)<-sampleNames
    classes<-as.data.frame(classes)

    xd<-x[(-1:-(headerRows-1)),]
    ids<-as.vector(t(xd[,1]))
    xd<-xd[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)
  }

  features<- dim(xd)[1]
  samples<- dim(xd)[2]
  geneNames<-rownames(xd)
  xd<-apply(xd,2,as.numeric)
  rownames(xd)<-geneNames
  colnames(xd)<-sampleNames

  if(!is.na(designFile)){
    x<-read.table(designFile,sep="\t", header=TRUE, row.names=1, fill=TRUE,
                  stringsAsFactors=FALSE)
    xd<-xd[,sort.list(colnames(xd))]
    xd<-xd[,colnames(xd) %in% rownames(x)]
    x<-x[rownames(x) %in% colnames(xd),]
    x<-x[sort.list(rownames(x)),]
    classes<-as.data.frame(x)
  }

  if(sum(apply(xd,2,is.na))>0 & impute){
    #library(impute)
    allAnn<-dimnames(xd)
    data.imputed<-impute.knn(as.matrix(xd))$data
    xd<-data.imputed[1:features,]
    dimnames(xd)<-allAnn
  }

  return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}



