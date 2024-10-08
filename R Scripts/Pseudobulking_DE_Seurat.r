library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(pheatmap)
library(limma)
library(utils)
library(stats)
library(grDevices)

#### SpaMTP Differential Peaks Analysis Functions ########################################################################################################################################################################################

#' Helper function for suppressing function progress messages
#'
#' @param message_text Character string containing the message being shown
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the messsage will be suppressed (default = TRUE).
#'
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

  return(summed)
}




#' Runs EdgeR analysis for pooled data
#'       - This function is used by run_edgeR_annotations()
#'
#' @param pooled_data A SingleCellExperiment object which contains the pooled pseudo-replicate data.
#' @param seurat_data A Seurat object containing the merged Xenium data being analysed (this is subset).
#' @param ident A character string defining the ident column to perform differential expression analysis against.
#' @param output_dir A character string defining the ident column to perform differential expression analysis against.
#' @param run_name A character string defining the title of this DE analysis (will be used when saving DEGs to .csv file).
#' @param n An integer that defines the number of pseudo-replicates per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param return.individual Boolean value defining whether to return a list of individual edgeR objects for each designated ident. If FALSE, one merged edgeR object will be returned (default = FALSE).
#'
#' @returns A modified edgeR object which contains the relative pseudo-bulking analysis outputs, including a DEGs data.frame with a list of differential expressed m/z peaks
#' @export
#'
#' @examples
#' # pooled_obj <- run_pooling(SeuratObj, "sample", n = 3)
#' # run_DE(pooled_obj, SeuratObj, "sample", "~/Documents/DE_output/", "run_1", n = 3, logFC_threshold = 1.2, annotation.column = "all_IsomerNames", assay = "Spatial")
run_DE <- function(pooled_data, seurat_data, ident, output_dir, run_name, n, logFC_threshold, assay, return.individual = FALSE, verbose = TRUE){

  verbose_message(message_text = paste("Running edgeR DE Analysis for ", run_name, " -> with samples [", paste(unique(unlist(seurat_data@meta.data[[ident]])), collapse = ", "), "]"), verbose = verbose)

  annotation_result <- list()
  for (condition in unique(seurat_data@meta.data[[ident]])){
    verbose_message(message_text = paste0("Starting condition: ",condition), verbose = verbose)

    groups <- SingleCellExperiment::colData(pooled_data)[[ident]]
    groups <- gsub(condition, "Comp_A", groups)
    groups <- ifelse(groups != "Comp_A", "Comp_B", groups)



    y <- edgeR::DGEList(SingleCellExperiment::counts(pooled_data), samples=SingleCellExperiment::colData(pooled_data)$orig.ident2, group = groups)


    y$samples$condition <- groups
    y$samples$ident <- sub("_(.*)", "", y$samples$samples)

    keep <- edgeR::filterByExpr(y, group = groups, min.count = 2, min.total.count = 10)
    y <- y[keep,]
    y <- edgeR::calcNormFactors(y)

    design <- stats::model.matrix(~groups)

    design[,2] <- 1-design[,2]

    y <- edgeR::estimateDisp(y, design, robust=TRUE)
    fit <- edgeR::glmQLFit(y, design, robust=TRUE)
    res <- edgeR::glmTreat(fit, coef=ncol(fit$design), lfc=log2(logFC_threshold))
    summary(limma::decideTests(res))
    res$table$regulate <- dplyr::recode(as.character(limma::decideTests(res)),"0"="Normal","1"="Up","-1"="Down")
    de_group_edgeR <- res$table[order(res$table$PValue),]
    table(limma::decideTests(res))
    res <- edgeR::topTags(res,n = nrow(y))
    res$table$regulate <- "Normal"
    res$table$regulate[res$table$logFC>0 & res$table$FDR<0.05] <- "Up"
    res$table$regulate[res$table$logFC<0 & res$table$FDR<0.05] <- "Down"
    de_group_edgeR <- res$table[order(res$table$FDR),]
    #table(de_group_edgeR$regulate)


    if (!(is.null(output_dir))){
      utils::write.csv(de_group_edgeR, paste0(output_dir,condition,"_",run_name, ".csv"))
    }

    de_group_edgeR$gene <- rownames(de_group_edgeR)

    y$DEGs <- de_group_edgeR
    annotation_result[[condition]] <- y

  }

  if (return.individual){
    annotation_result
    return(annotation_result)
  } else {



    edger <- edgeR::DGEList(
      counts = annotation_result[[1]]$counts,
      samples = annotation_result[[1]]$samples
    )
    edger$samples$group <- edger$samples$ident
    edger$samples$condition <- NULL

    degs <- lapply(names(annotation_result), function(x){
      annotation_result[[x]]$DEGs$cluster <- x
      rownames(annotation_result[[x]]$DEGs) <- NULL
      annotation_result[[x]]$DEGs
    })

    combined_degs <- do.call(rbind, degs)
    rownames(combined_degs) <- 1:length(combined_degs$cluster)

    edger$DEGs <- combined_degs
    return(edger)
  }

}


#' Finds all differentially expressed genes between comparison groups specified
#'       - This function uses run_pooling() and run_DE() to pool and run EdgeR analysis
#'
#' @param data A Seurat object
#' @param ident A character string defining the metadata column or groups to compare genes values between.
#' @param n An integer that defines the number of pseudo-replicates (pools) per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param DE_output_dir A character string defining the directory path for all output files to be stored. This path must a new directory. Else, set to NULL as default.
#' @param run_name A character string defining the title of this DE analysis that will be used when saving DEDs to .csv file (default = 'FindAllDEGs').
#' @param assay A character string defining the assay where the count data and annotations are stored (default = "Spatial").
#' @param slot Character string defining the assay storage slot to pull the relative expression values from. Note: EdgeR requires raw counts, all values must be positive (default = "counts").
#' @param return.individual Boolean value defining whether to return a list of individual edgeR objects for each designated ident. If FALSE, one merged edgeR object will be returned (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns Returns an list() contains the EdgeR DE results. Pseudo-bulk counts are stored in $counts and DEGs are in $DEGs.
#' @export
#'
#' @examples
#' # FindAllDEGs(SeuratObj, "sample",DE_output_dir = "~/Documents/DE_output/")
FindAllDEGs <- function(data, ident, n = 3, logFC_threshold = 1.2, DE_output_dir = NULL, run_name = "FindAllDEGs", assay = "Spatial", slot = "counts", return.individual = FALSE, verbose = TRUE){

  if (!(is.null(DE_output_dir))){
    if (dir.exists(DE_output_dir)){
      warning("Please supply a directory path that doesn't already exist")
      stop("dir.exists(DE_output_dir) = TRUE")
    } else{
      dir.create(DE_output_dir)
    }
  }

  #Step 1: Run Pooling to split each unique ident into 'n' number of pseudo-replicate pools
  pooled_data <- run_pooling(data,ident, n = n, assay = assay, slot = slot, verbose = verbose)

  #Step 2: Run EdgeR to calculate differentially expressed m/z peaks
  DEG_results <- run_DE(pooled_data, data, ident = ident, output_dir = DE_output_dir, run_name = run_name, n=n, logFC_threshold=logFC_threshold, assay = assay, verbose = verbose, return.individual = return.individual)

  # Returns an EDGEr object which contains the pseudo-bulk counts in $counts and DEGs in $DEGs
  return(DEG_results)

}




#' Generates a Heatmap of DEGs generated from edgeR analysis run using FindAllDEGs().
#'       - this function uses pheatmap() to plot data
#'
#' @param edgeR_output A list containing outputs from edgeR analysis (from FindAllDEGs()). This includes pseudo-bulked counts and DEGs.
#' @param n A numeric integer that defines the number of UP and DOWN regulated peaks to plot (default = 25).
#' @param only.pos Boolean indicating if only positive markers should be returned (default = FALSE).
#' @param FDR.threshold Numeric value that defines the FDR threshold to use for defining most significant results (default = 0.05).
#' @param logfc.threshold Numeric value that defines the logFC threshold to use for filtering significant results (default = 0.5).
#' @param order.by Character string defining which parameter to order markers by, options are either 'FDR' or 'logFC' (default = "FDR").
#' @param scale A character string indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param color A vector of colors used in heatmap (default = grDevices::colorRampPalette(c("navy", "white", "red"))(50)).
#' @param cluster_cols Boolean value determining if columns should be clustered or hclust object (default = F).
#' @param cluster_rows Boolean value determining if rows should be clustered or hclust object (default = T).
#' @param fontsize_row A numeric value defining the fontsize of rownames (default = 15).
#' @param fontsize_col A numeric value defining the fontsize of colnames (default = 15).
#' @param cutree_cols A numeric value defining the number of clusters the columns are divided into, based on the hierarchical clustering(using cutree), if cols are not clustered, the argument is ignored (default = 9).
#' @param silent Boolean value indicating if the plot should not be draw (default = TRUE).
#' @param plot_annotations_column Character string indicating the column name that contains the metabolite annotations to plot. Annotations = TRUE must be used in FindAllDEGs() for edgeR output to include annotations. If plot_annotations_column = NULL, m/z vaues will be plotted (default = NULL).
#' @param save_to_path Character string defining the full filepath and name of the plot to be saved as.
#' @param plot.save.width Integer value representing the width of the saved pdf plot (default = 20).
#' @param plot.save.height Integer value representing the height of the saved pdf plot (default = 20).
#' @param nlabels.to.show Numeric value defining the number of annotations to show per m/z (default = NULL).
#' @param annotation_colors List of colour palette used to colour idents. If NULL default colour scheme will be used (default = NULL).
#'
#' @returns A heatmap plot of significantly differentially expressed peaks defined in the edgeR ouput object.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' # DEGs <- FindAllDEGs(SeuratObj, "sample")
#'
#' # DEGsHeatmap(DEGs)
DEGsHeatmap <- function(edgeR_output,
                         n = 5,
                         only.pos = FALSE,
                         FDR.threshold = 0.05,
                         logfc.threshold = 0.5,
                         order.by = "FDR",
                         scale ="row",
                         color = grDevices::colorRampPalette(c("navy", "white", "red"))(50),
                         cluster_cols = F,
                         cluster_rows = T,
                         fontsize_row = 15,
                         fontsize_col = 15,
                         cutree_cols = 9,
                         silent = TRUE,
                         save_to_path = NULL,
                         plot.save.width = 20,
                         plot.save.height = 20,
                         annotation_colors = NULL){


  degs <- edgeR_output$DEGs
  degs <- subset(degs, FDR < FDR.threshold)

  if (order.by == "FDR"){

    grouped_pos<- degs %>%
      group_by(cluster) %>%
      filter( logFC > logfc.threshold) %>%
      arrange(desc(regulate)) %>%
      slice_head(n = n)


    if (only.pos) {
      grouped_neg <- NULL

    } else {
      grouped_neg <- degs %>%
        group_by(cluster) %>%
        filter(logFC < - logfc.threshold) %>%
        arrange(regulate) %>%
        slice_head(n = n)
    }
    df <- do.call(rbind, list(grouped_pos,grouped_neg))
    df <- df[order(df$cluster, dplyr::desc(df$regulate)), ]

  } else {
    if ( order.by != "logFC"){
      warning("order.by has invalid argument. Must be either 'FDR' or 'logFC'. Heatmap defaulting to order by logFC")
    }

    grouped_pos<- degs %>%
      group_by(cluster) %>%
      filter(logFC > logfc.threshold) %>%
      arrange(-logFC) %>%
      slice_head(n = n)


    if (only.pos) {
      grouped_neg <- NULL
    } else {
      grouped_neg <- degs %>%
        group_by(cluster) %>%
        filter(logFC < - logfc.threshold) %>%
        arrange(logFC) %>%
        slice_head(n = n)
    }
    df <- do.call(rbind, list(grouped_pos,grouped_neg))
    df <- df[order(df$cluster, -df$logFC), ]
  }



  col_annot <- data.frame(sample = edgeR_output$samples$ident)
  row.names(col_annot) <- colnames(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE)))

  col_annot <- data.frame(sample = edgeR_output$samples$ident)
  row.names(col_annot) <- colnames(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE)))
    if (!is.null(annotation_colors)){
       annotation_colors <- list(sample = unlist(annotation_colors))
  } else {
      annotation_colors <- NA
  }
  
  mtx <- as.matrix(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE))[unique(df$gene),])

  p <- pheatmap::pheatmap(mtx,scale=scale,color=color,cluster_cols = cluster_cols, annotation_col=col_annot, cluster_rows = cluster_rows,
                          fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_cols = cutree_cols, silent = silent, annotation_colors = annotation_colors)

   if (!(is.null(save_to_path))){
     save_pheatmap_as_pdf(pheatmap = p, filename = save_to_path, width = plot.save.width, height = plot.save.height)
   }

  return(p)
}


#' Saves the DEGs generated pheatmap as a PDF
#'
#' @param pheatmap A pheatmap plot object that is being saved.
#' @param filename Character string defining the full filepath and name of the plot to be saved as.
#' @param width Integer value representing the width of the saved pdf plot (default = 20).
#' @param height Integer value representing the height of the saved pdf plot (default = 20).
#'
#' @export
#'
#' @examples
#' # save_pheatmap_as_pdf(pheatmap, filename = "/Documents/plots/pheatmap1")
save_pheatmap_as_pdf <- function(pheatmap, filename, width=20, height=20){

  pdf(paste0(filename,".pdf"), width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(pheatmap$gtable)
  dev.off()
}

########################################################################################################################################################################################################################
