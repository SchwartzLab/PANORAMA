#' Plot metrics NA rate
#'
#' @param PANORAMA PANORAMA object
#'
#' @return ggplot
#' @export
qc_plotMetricsNArate <- function(PANORAMA){
    SETS <- unique(PANORAMA$META$set)
    NUCS <- c("A", "C", "G", "T")
    METRICS <- unique(PANORAMA$RES$metric)
    tmpL <- lapply(SETS, function(set_i){
        lapply(METRICS, function(metric_i){
            lapply(NUCS, function(nuc_i){
                tmpO <- mean(is.na(dplyr::filter(PANORAMA$RES, refSeq == nuc_i,
                                          metric == metric_i, set == set_i)$score))
                data.table::data.table(set = set_i, refSeq = nuc_i, metric = metric_i, NA_rate = tmpO)
            }) %>% do.call(what = rbind)
        }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)
    ggplot(tmpL, aes(x = metric, y = set, fill = NA_rate)) + ggplot2::geom_tile() + ggplot2::facet_wrap(.~refSeq) +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' Plot predictions NA rate
#'
#' @param PANORAMA PANORAMA object
#'
#' @return ggplot
#' @export
qc_plotPredsNApct <- function(PANORAMA){
    tmpO <- qc_predsNApct(PANORAMA)
    ggplot(tmpO, aes(x = RNAmod, y = NApred_pct, fill = set)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        theme_bw() + theme(legend.position = "bottom")
}

#' QC Predictions NA percentage
#'
#' @param PANORAMA PANORAMA object
#'
#' @return data.table
#' @export
qc_predsNApct <- function(PANORAMA){
    if("CALLS" %in% names(PANORAMA)){
        sets <- names(PANORAMA$CALLS)
        tmpO <- lapply(sets, function(set_i){
            RNAmods <- names(PANORAMA$CALLS[[set_i]])
            lapply(RNAmods, function(RNAmod_i){
                tmpDT <- PANORAMA$CALLS[[set_i]][[RNAmod_i]]
                naRate <- tmpDT[tmpDT$refSeq %in% RNAMod_baseNuc[[RNAmod_i]],]$pred %>%
                    is.na() %>% mean() %>% multiply_by(100) %>% round(2)
                data.table(set = set_i, RNAmod = RNAmod_i, NApred_pct = naRate)
            }) %>% do.call(what = "rbind")
        }) %>% do.call(what = "rbind")
        tmpO
    }else{stop("CALLS is not an element of PANORAMA object. Use ",
               "panorama_assignScores() and call_RNAmods_logRegMods() to compute ",
               "RNAmod predictions")}
}
