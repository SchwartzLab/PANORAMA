#' Boxplots of $RES calculated metrics grouping by modified nucleotides
#'
#' @param PANORAMA list. PANORAMA object as output by \code{\link{panorama_PANORAMA}}
#' @param metrics character. Selected metrics to plot. By default plots all metrics
#' @param title character. Title to be shown on top of plot
#'
#' @return
#' @export
#'
#' @examples
panorama_metricsBoxPlot_byNuc <- function(PANORAMA, metrics = NULL, title = ""){
    if(!"nuc" %in% names(PANORAMA$RES)){stop("'nuc' column with RNAmod identity is not present in PANORAMA$RES")}
    if(!is.null(metrics)){
        tmpDT <- dplyr::filter(PANORAMA$RES, metric %in% metrics) %>% stats::na.omit()
    }else{
        tmpDT <- PANORAMA$RES %>% stats::na.omit()
    }
    ggplot2::ggplot(tmpDT, ggplot2::aes(x = nuc, y = score)) + ggplot2::geom_boxplot(outlier.colour = NA) +
        ggplot2::geom_point(ggplot2::aes(colour = metric), alpha = 0.2) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::ggtitle(title)
}

#' PANORAMA scatterplot of scores by txcoor
#'
#' Plot scatterplots of scores in PANORAMA$RES by txcoor and colored by gene
#'
#' @param PANORAMA list. PANORAMA object as output by \code{\link{panorama_PANORAMA}}
#' @param metrics character. Selected metrics to plot. By default plots all metrics
#' @param title character. Title to be shown on top of plot
#'
#' @return
#' @export
#'
#' @examples
panorama_metricsScatterPlot_byPos <- function(PANORAMA, metrics = NULL, title = "", na.rm = FALSE){
    if(!is.null(metrics)){
        tmpDT <- dplyr::filter(PANORAMA$RES, metric %in% metrics)
    }else{
        tmpDT <- PANORAMA$RES
    }
    if(na.rm){tmpDT <- stats::na.omit(tmpDT)}
    ggplot2::ggplot(tmpDT) +
        ggplot2::geom_point(ggplot2::aes(x = .data$pos, y = .data$score, colour = .data$gene), alpha = 0.2) +
        ggplot2::facet_grid(.data$set ~ .data$metric, scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::ggtitle(title) + ggplot2::xlab("txCoor") +
        ggplot2::theme_minimal()
}

#' Confussion matrix plot
#'
#' @param PANORAMA
#' @param title
#' @param pointAlpha
#' @param removeOtherRNAmods
#' @param rmNotModGenes
#' @param blacklist
#'
#' @return
#' @export
#'
panorama_plot_confMat <- function(PANORAMA, title = "", pointAlpha = 0.5,
                               removeOtherRNAmods = TRUE, rmNotModGenes = TRUE,
                               blacklist = NULL){
    df <- panorama_confMatrix4plot(PANORAMA, removeOtherRNAmods = removeOtherRNAmods,
                                rmNotModGenes = rmNotModGenes,
                                blacklist = blacklist)
    df$RNAmod <- factor(df$RNAmod, levels = rev(RNAmods_vec))
    ggplot(df, aes(x = outcome, y = RNAmod)) +
        geom_jitter(alpha = 0.5, width = 0.2, height = 0.2) +
        theme_bw() + ggtitle(title)
}
