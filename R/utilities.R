#' Results summary
#'
#' @param PANORAMA
#' @param NCORES
#'
#' @return
#' @export
resSummary <- function(PANORAMA, NCORES = 1){
    minMETA <- miniMETA(PANORAMA)
    SETS <- names(PANORAMA$CALLS)
    RNAmod_InPANORAMA <- lapply(SETS, function(set){
        names(PANORAMA$CALLS[[set]])
    }) %>% unlist() %>% unique()
    RNAmod_withPreds <- lapply(RNAmod_InPANORAMA, function(RNAmod){
        lapply(SETS, function(set){
            data.frame(set = set, RNAmod = RNAmod,
                       anyPred = any(!is.na(PANORAMA$CALLS[[set]][[RNAmod]]$pred)))
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind") %>% data.table::data.table()
    RNAmod_withPreds <- unique(RNAmod_withPreds[RNAmod_withPreds$anyPred == TRUE, ][["RNAmod"]])
    summRES <- parallel::mclapply(mc.cores = NCORES, SETS, function(set){
        tmpO <- lapply(RNAmod_withPreds, function(RNAmod){
            tmp <- PANORAMA$CALLS[[set]][[RNAmod]][, c(panorama_baseCoorCols,
                                                    "linear_Score", "pred"), with = FALSE]
            tmp$pred <- factor(ifelse(tmp$pred, yes = "putModified", no = "notModified"))
            tmp <- tmp[!is.na(tmp$linear_Score),]
            tmp$RNAmod <- RNAmod
            tmp$set <- set
            data.table::data.table(tmp)
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind") %>% dplyr::left_join(minMETA, by ="set")
    summRES$bioTreat <- factor_numOrd(summRES$bioTreat)
    summRES$replicate <- factor_numOrd(summRES$replicate)
    summRES$RNAmod <- factor(summRES$RNAmod, levels = PANORAMA::RNAmods_vec)
    summRES
}

#' Calls summary
#'
#' @param PANORAMA list. PANORAMA object
#'
#' @return
#' @export
callsSummary <- function(PANORAMA){
    tmpM <- miniMETA(PANORAMA)
    SETS <- names(PANORAMA$CALLS)
    RNAmod_InPANORAMA <- lapply(SETS, function(set){
        names(PANORAMA$CALLS[[set]])
    }) %>% unlist() %>% unique()
    RNAmod_withPreds <- lapply(RNAmod_InPANORAMA, function(RNAmod){
        lapply(SETS, function(set){
            data.frame(set = set, RNAmod = RNAmod,
                       anyPred = any(!is.na(PANORAMA$CALLS[[set]][[RNAmod]]$pred)))
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind") %>% data.table::data.table()
    RNAmod_withPreds <- unique(RNAmod_withPreds[RNAmod_withPreds$anyPred == TRUE, ][["RNAmod"]])
    callSumm <- lapply(SETS, function(set){
        lapply(RNAmod_withPreds, function(RNAmod){
            calls <- subset(PANORAMA$CALLS[[set]][[RNAmod]], refSeq %in% RNAMod_baseNuc[[RNAmod]])
            data.table(set = set,
                       RNAmod = RNAmod,
                       putSites = sum(calls$pred, na.rm = TRUE),
                       queried = sum(!is.na(calls$pred)),
                       queriedRate = round(mean(!is.na(calls$pred)), digits = 4),
                       putGenes = length(unique(base::subset(calls, pred == TRUE)$gene)),
                       queriedGenes = length(unique(calls[!is.na(calls$pred),]$gene)))
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind")
    callSumm$modRate <- round(callSumm$putSites / callSumm$queried, 4)
    callSumm$modRate[is.nan(callSumm$modRate)] <- NA
    callSumm <- dplyr::left_join(callSumm, tmpM, by ="set")
    callSumm$bioTreat <- factor_numOrd(callSumm$bioTreat)
    callSumm$replicate <- factor_numOrd(callSumm$replicate)
    callSumm$RNAmod <- factor(callSumm$RNAmod, levels = PANORAMA::RNAmods_vec)
    callSumm
}


#' Order factor's levels numerically
#'
#' Extracts the first numeric portion of factors' levels to order them numerically
#'
#' @param fct factor. Factor to re-leveled
#' @param decreasing logical. Order decreasingly, default = FALSE
#'
#' @return
#' @export
#'
#' @examples
#' factor_numOrd(c("100deg", "25deg", "50deg"))
factor_numOrd <- function(fct, decreasing = FALSE){
    lvls <- unique(fct)
    fctNumeric <- as.numeric(stringr::str_extract(lvls, "[[:digit:]]+"))
    factor(fct, levels = lvls[order(fctNumeric, decreasing = decreasing)])
}

#' Mini-Meta (with ordered treatments)
#'
#' Smaller version of experimental design for specific purposes
#'
#' @param PANORAMA list. PANORAMA object
#'
#' @return
#' @export
miniMETA <- function(PANORAMA){
    META <- PANORAMA$META
    BIOTREATS  <- sort(unique(META$bioTreat))
    REPLICATES <- sort(unique(META$replicate))
    SETS <- names(PANORAMA$CALLS)
    miniMETA <- data.table::data.table(set = SETS)
    if(!is.null(REPLICATES)){ #Order factors by numeric
        miniMETA$replicate <- unlist(lapply(miniMETA$set, function(set_i){
            unique(META[META$set == set_i,]$replicate)
        })) %>% factor_numOrd()
    }else{
        miniMETA$replicate <- "rep1"
    }
    if(!is.null(BIOTREATS)){ #Order factors by numeric
        miniMETA$bioTreat <- unlist(lapply(miniMETA$set, function(set_i){
            unique(META[META$set == set_i,]$bioTreat)
        })) %>% factor_numOrd()
    }else{
        miniMETA$bioTreat <- "none"
    }
    miniMETA[order(miniMETA$bioTreat, miniMETA$replicate),]
}
