# Magrittr Pipe Operator
`%>%` <- magrittr::`%>%`

# Make Temporary directory for Pan-Mod-seq processing
mkTmpDir <- function(){
    if(!dir.exists("./PANORAMA_tmpDir")){dir.create("./PANORAMA_tmpDir")}
}

#' Remove PANORAMA temporary directory
#'
#' Removes temporary directory "./PANORAMA_tmpDir" created for alignment steps.
#'
#' @return
#' @export
#'
#' @examples
rmTmpDir <- function(){
    if(dir.exists("./PANORAMA_tmpDir")){unlink("./PANORAMA_tmpDir", recursive = TRUE)}
}

# Extract gene sequences from genome and geneAnnotation
getGeneSeqsfromGenome <- function(geneAnnot, genome, nCores = 1){
    txtools:::check_mc_windows(nCores)
    txtools:::check_GA_genome_chrCompat(geneAnnot = geneAnnot, genome = genome)
    parallel::mclapply(mc.cores = nCores, seq_along(geneAnnot), function(i){
        selGene <- geneAnnot[i]
        iGene <- as.character(selGene$name)
        iChr <- as.character(GenomicRanges::seqnames(selGene))
        iStr <- as.character(selGene@strand)
        iGA <- selGene
        iBlocks <- S4Vectors::mcols(iGA)$blocks %>% txtools:::if_IRangesList_Unlist() %>%
            IRanges::shift(IRanges::start(iGA) - 1)
        SEQ <- stringr::str_sub(as.character(genome[[iChr]]), start = IRanges::start(iBlocks),
                                end = IRanges::end(iBlocks)) %>% paste(collapse = "") %>%
            Biostrings::DNAString()
        if(iStr == "-") {
            SEQ <- Biostrings::reverseComplement(SEQ)
        }
        SEQ
    }) %>% Biostrings::DNAStringSet()
}

# Add results to a PANORAMA object. Remove scores if metric is already present.
hlpr_add_REScols <- function(PANORAMA_RES, REScols){
    iMetric <- unique(REScols[,"metric"]) %>% as.character()
    # remove results for identical metric
    if("metric" %in% names(PANORAMA_RES)){
        PANORAMA_RES <- PANORAMA_RES[PANORAMA_RES$metric != iMetric,]
    }
    PANORAMA_RES <- rbind(PANORAMA_RES, REScols)
    PANORAMA_RES
}

#' Check RNAmods
#'
#' Check which RNAmods are possible to predict based on RNAmod metrics present
#' in PANORAMA object
#'
#' @param PANORAMA
#'
#' @return character
#' @export
check_whichRNAmods <- function(PANORAMA){
    if("RES" %in% names(PANORAMA)){
        tmp <- lapply(PANORAMA::metricsList, function(x){
            any(unique(PANORAMA$RES$metric) %in% x)
        }) %>% unlist()
        return(names(tmp[tmp]))
    }else{stop("PANORAMA object does not include 'RES' element")}
}

# Confusion matrix helper
hlpr_confusMatOutcome <- function(x){
    x$outcome <- paste0(ifelse(x$truth == x$pred, "T", "F"),
                        ifelse(as.logical(x$pred), "P", "N"))
    x$outcome[is.na(x$truth) | is.na(x$pred)] <- "NA"
    y <- x[x$outcome == "NA",]
    NAsum <- tapply(y$Freq, y$refSeq, "sum")
    out <- rbind(x[x$outcome %in% c("TN", "TP", "FP", "FN") & !is.na(x$refSeq),],
                 data.table::data.table(truth = NA, pred = NA,
                                        refSeq = names(NAsum),
                                        Freq = as.numeric(NAsum),
                                        outcome = "NA"))
    out$Freq[is.na(out$Freq)] <- 0
    return(out)
}

hlp_removeColumnIfPresent <- function(DT, colName){
    if(colName %in% colnames(DT)){
        DT[, colnames(DT) != colName, with = FALSE]
    }else{DT}
}

check_DT <- function(DT){
    if(!data.table::is.data.table(DT)){
        if(!is.data.frame(DT)){
            stop("DT must be either a data.table or a data.frame")
        }else{
            if(is.data.frame(DT)){DT <- data.table::data.table(DT)}
        }
    }
    return(DT)
}

myROC <- function(values, success, na.rm = TRUE){
    if(na.rm){
        tmpL <- !(is.na(values) | is.na(success))
        values <- values[tmpL]
        success <- success[tmpL]
    }
    AUC::roc(values, factor(as.numeric(success)))
}

# Calculate mid-point between threshold and median linear score of modified sites
hlp_calcMidpoint <- function(thr, medLs){
    C <- medLs - thr
    D <- C / 2
    thr + D
}

# Remove empty strings from character vector
rm_empty <- function(x){
    if(!is.character(x)){stop("x must be of class character")}
    x[x != ""]
}

mcc <- function(tp, tn, fp, fn){
    (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
}


#' Add biotype (rRNA or tRNA)
#'
#' Adds RNA gene biotype based on gene name patterns
#'
#' @param txDT data.table
#' @param geneInfo data.frame Table with columns gene and biotype
#'
#' @return txDT
#' @export
add_biotype <- function(txDT, geneInfo){
    if(is.null(geneInfo)){
        txDT
    }else{
        geneInfo <- dplyr::rename(geneInfo, gene = name)
        geneInfo <- dplyr::rename(geneInfo, gene_biotype = biotype)
        geneInfo$gene_biotype <- factor(geneInfo$gene_biotype)
        hlp_removeColumnIfPresent(txDT, colName = "gene_biotype") %>%
            dplyr::left_join(dplyr::select(geneInfo, c("gene", "gene_biotype")), by = "gene")
    }
}

# For multicall RNAmod detected decision
remove_Nm <- function(x){
    stringr::str_remove(x, "Am") %>%
        stringr::str_remove("Cm") %>%
        stringr::str_remove("Gm") %>%
        stringr::str_remove("Um")
}

remove_trailCommas <- function(x){
    stringr::str_remove(x, ",+$")
}

# Return all max instances
which.max_all <- function(x){
    tmpM <- max(x, na.rm = TRUE)
    which(x == tmpM)
}

limit_0_100 <- function(x){
    x[x < 0] <- 0
    x[x > 100] <- 100
    x
}

# Return all max instances
which.max_all <- function(x){
    tmpM <- max(x, na.rm = TRUE)
    which(x == tmpM)
}

# NArate metric by nuc in PANORAMA training
NArate_metricNuc <- function(PANORAMA, metrics = NULL, nucs = PANORAMA::RNAmods_vec){
    if(is.null(metrics)){metrics <- unique(PANORAMA$RES$metric)}
    lapply(metrics, function(metric_i){
        lapply(nucs, function(RNAmod_i){
            tmpNAr <- dplyr::filter(PANORAMA$RES, metric == metric_i & nuc == RNAmod_i)[["score"]] %>%
                is.na() %>% mean()
            tmpQueried <- dplyr::filter(PANORAMA$RES, metric == metric_i & nuc == RNAmod_i)[["score"]] %>%
                is.na() %>% magrittr::not() %>% sum()
            data.table::data.table(metric = metric_i, RNAmod = RNAmod_i, NArate = tmpNAr, measuredSites = tmpQueried)
        }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)
}

# Output RNAmods with predictions in calls per input set name
check_whichRNAmods_pred <- function(PANORAMA, set_i){
    RNAMODS <- names(PANORAMA$CALLS[[set_i]])
    tmpL <- purrr::map_lgl(RNAMODS, function(RNAmod_i){
        "pred" %in% colnames(PANORAMA$CALLS[[set_i]][[RNAmod_i]])
    })
    RNAMODS[tmpL]
}
