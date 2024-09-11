# Altered Function to work for merged seurat objects with multiple images/fov's

## This builds off the Gene set co-regulation analysis tutorial found here: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html#analysis-of-spatial-rna-seq 


## Libraries required 
library(msigdbr)
library(fgsea)

### Functions ###

#' Plots the spatial expression profile of a gene set
#'
#' @param pathway Gene set to plot or a list of gene sets (see details below)
#' @param object Seurat object
#' @param images Character string defining which image to use for plotting. 
#' @param title plot title
#' @param assay assay to use for obtaining scaled data, preferably with
#' the same universe of genes in the scaled data
#' @param colors vector of three colors to use in the color scheme
#' @param guide option for `ggplot2::scale_color_gradientn` to control for presence of the color legend
#' the same universe of genes in the scaled data
#'
#' @return ggplot object (or a list of objects) with the coregulation profile plot
#'
#' When the input is a list of pathways, pathway names are used for titles.
#' A list of ggplot objects a returned in that case.
#
#' @import ggplot2
#' @export
PlotHallmarksSpatially <- function(pathway,
                                           object,
                                            images, 
                                           title=NULL,
                                           assay=DefaultAssay(object),
                                           colors=c("darkblue", "lightgrey", "darkred"),
                                           guide="colourbar") {
    stopifnot(requireNamespace("Seurat"))
    # TODO duplicated code with plotCoregulationProfileReduction
    if (is.list(pathway)) {
        if (is.null(title)) {
            titles <- names(pathway)
        } else {
            if (length(title) != length(pathway)) {
                stop("Length of the specified titles does not match count of pathways")
            }
            titles <- title
        }
        ps <- lapply(seq_along(pathway), function(i)
            plotCoregulationProfileSpatialX(pathway[[i]],
                                           object=object,
                                           title=titles[i],
                                           images = images,
                                           assay=assay,
                                           colors=colors))
        names(ps) <- names(pathway)
        return(ps)
    }



    obj2 <- addGesecaScores(list(pathway=pathway), object, assay=assay,
                            scale=TRUE)

    p <- Seurat::SpatialFeaturePlot(obj2, features = "pathway", images = images,
                                    combine = FALSE, image.alpha = 0)[[1]]

    p$scales$scales[p$scales$find("fill")] <- NULL

    # suppress message of replacing existing color palette
    suppressMessages({
        p2 <- p +
            scale_fill_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
                                 oob=scales::squish,
                                 colors=colors,
                                 guide = guide,
                                 name = "z-score"
            ) + theme(legend.position = theme_get()$legend.position)
    })

    if (!is.null(title)) {
        p2 <- p2 + ggtitle(title)
    }
    p2
}






#' Additional Helper function for 'PlotHallmarksSpatially()'
#'
#' @param pathways Gene set to plot or a list of gene sets (see details below)
#' @param object Seurat object
#' @param assay assay to use for obtaining scaled data, preferably with the same universe of genes in the scaled data
#' @param prefix Character string defining what to put before the metadata_column title (default = "").
#' @param scale Boolean value indicating whether to scale the GESECA score (default = FALSE)
#'
#' @return seurat object with the GESECA score added to metadata
addGesecaScores <- function(pathways,
                            object,
                            assay=DefaultAssay(object),
                            prefix="",
                            scale=FALSE) {
    x <- GetAssay(object, assay)
    E <- x@scale.data

    res <- object


    for (i in seq_along(pathways)) {
        pathway <- pathways[[i]]
        pathway <- intersect(unique(pathway), rownames(E))

        score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
        score <- scale(score, center=TRUE, scale=scale)
        res@meta.data[[paste0(prefix, names(pathways)[i])]] <- score
    }


    return(res)
}





### Example ###

## General steps before plotting spatially
#set.seed(1)
#pathwaysDF <- msigdbr("human", category="H")
#pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
#E <- vis.combined@reductions$pca@feature.loadings #where vis.combined is your merged seurat obejct
#gesecaRes <- geseca(pathways, E, minSize = 15, maxSize = 500, center = FALSE)

## Use case of PlotHallmarksSpatially()
# topPathways <- gesecaRes[, pathway] |> head(2)
# titles <- sub("HALLMARK_", "", topPathways)


# options(repr.plot.width = 30, repr.plot.height = 8) #used to change the plot size only need in jupyter notebook

# for (i in names(vis.combined@images)){
    
#     ps <- PlotHallmarksSpatially(pathways[topPathways], vis.combined, image = i,
#                                        title=titles)
#     print(cowplot::plot_grid(plotlist=ps, ncol=2))
# }
