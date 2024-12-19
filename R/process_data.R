panorama_summary <- function(PANORAMA){
    RNAmods <- names(PANORAMA$CALLS[[1]])
    OUT <- lapply(RNAmods, function(RNAmod_i){
        tmp_out <- PANORAMA$DATA[[1]][,panorama_baseCoorCols, with = FALSE]
        SETS <- levels(PANORAMA$META$set)
        if(is.null(PANORAMA$CALLS[[SETS[1]]][[RNAmod_i]]$linear_Score)){
            warning("linear_Score was not found for ", RNAmod_i, ". Results will",
                    " not be summarized in PANORAMA$SUMMARY. Check if metrics where ",
                    "placed in PANORAMA$CALLS[[1]]$", RNAmod_i)
            return(NULL)
        }else{
            tmp_lin <- lapply(SETS, function(x){
                PANORAMA$CALLS[[x]][[RNAmod_i]]$linear_Score
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("linScore",
                                                                          SETS, sep = "_"))
            tmp_prd <- lapply(SETS, function(x){
                PANORAMA$CALLS[[x]][[RNAmod_i]]$pred
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("pred",
                                                                          SETS, sep = "_"))
            tmp_out <- cbind(tmp_out, tmp_lin, tmp_prd)
            return(tmp_out)
        }
    }) %>% magrittr::set_names(RNAmods)
    PANORAMA$SUMMARY <- OUT
    PANORAMA
}

# Tables of files from FASTQ to expected BAM and RDS targets
files_table <- function(META){
    fastq <- META$FASTQ
    bam <- META$BAM
    lce <- gsub(pattern = "Aligned.out.sorted.bam",
                replacement = "Aligned.out.sorted.lce.txt", META$BAM)
    rds <- gsub(pattern = "Aligned.out.sorted.bam",
                replacement = "Aligned.out.sorted.txDT.rds", META$BAM)
    # Table
    tmpDT <- data.table::data.table(FASTQ = c(fastq),
                                    BAM = c(bam),
                                    BAM_ok = file.exists(c(bam)),
                                    lce = c(lce),
                                    lce_ok = file.exists(c(lce)),
                                    rds = c(rds),
                                    rds_ok = file.exists(c(rds)))
    return(tmpDT)
}

# FASTQ Duplication rate (library complexity)
fastq_dupRate <- function(FASTQs_pahts, nCores){
    parallel::mclapply(mc.cores = nCores, FASTQs_pahts, function(file){
        tmp <- ShortRead::readFastq(file)
        tmp2 <- ShortRead::readFastq(gsub(file, pattern = "R1", replacement = "R2"))
        dupR <- paste(tmp@sread, tmp2@sread, sep = "") %>% duplicated() %>% mean
        return(dupR)
    }) %>% unlist
}

#' Nucleotide frequency report
#'
#' @param META
#' @param nCores
#' @param firstN
#' @param pairedEnd
#'
#' @return
#' @export
#'
#' @examples
fastq_nucFreq <- function(META, nCores, firstN = 1e4, pairedEnd = TRUE){
    parallel::mclapply(mc.cores = nCores, META$FASTQ, function(file){
        tmp <- readLines(file, firstN * 4)[seq(2, firstN*4, 4)] %>% Biostrings::DNAStringSet()
        r2File <- gsub(file, pattern = "R1", replacement = "R2")
        if(pairedEnd){
            tmp2 <- readLines(r2File, firstN * 4)[seq(2, firstN*4, 4)] %>%
                Biostrings::DNAStringSet() %>% Biostrings::complement()
            mR1R2 <- paste(tmp, tmp2, sep = "")
            nucFreq <- mR1R2 %>% as.character() %>% stringr::str_split(pattern = "") %>% unlist %>% table
        }else if(!pairedEnd){
            nucFreq <- tmp %>% as.character() %>%  stringr::str_split(pattern = "") %>% unlist %>% table
        }
        return(nucFreq[c("A", "C", "G", "T")])
    }) %>% do.call(what = cbind) %>% magrittr::set_colnames(META$id)
}

#' nucleotide frequency plot
#'
#' @param nucF_x
#' @param subtitle
#'
#' @return
#' @export
#'
#' @examples
gg_nucFreq <- function(nucF_x, subtitle){
    tmp <- prop.table(nucF_x, margin = 2) %>% data.frame %>%
        tibble::rownames_to_column(var = "nuc") %>%
        tidyr::pivot_longer(cols = -nuc, names_to = "Sample", values_to = "Ratio")
    ggplot2::ggplot(tmp, ggplot2::aes(x = Sample, y = Ratio, fill = nuc)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_brewer(palette="Set1") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle("Nucleotide Frequency per library", subtitle = subtitle)
}

#' Alignment report table
#'
#' @param META
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
reads_report <- function(META, nCores = 4){
    DT <- files_table(META)
    if(all(DT$BAM_ok) & all(DT$rds_ok)){
        res1 <- parallel::mclapply(mc.cores = nCores, DT$FASTQ, function(x){
            tmp <- ShortRead::readFastq(x)
            length(tmp)
        }) %>% unlist
        res2 <- parallel::mclapply(mc.cores = nCores, DT$BAM[DT$BAM_ok], function(x){
            tmp2 <- Rsamtools::scanBam(x)
            tmp2[[1]]$qname %>% unique %>% length
        }) %>% unlist
        res3 <- parallel::mclapply(mc.cores = nCores, DT$rds[DT$rds_ok], function(x){
            tmplog <- data.table::fread(gsub(pattern = "rds", "log", x), header = F)
            tmplog[grep(tmplog$V1, pattern =  "unique reads"),]$V2 %>% as.numeric()
        }) %>% unlist
        tmpDT <- data.table::data.table(sample = META$id,
                                        FASTQ_reads = res1,
                                        BAM_aligns = NA,
                                        tx_starts = NA)
        tmpDT$BAM_aligns[DT$BAM_ok] <- res2
        tmpDT$tx_starts[DT$rds_ok] <- res3
        tmpDT$pC_BAM <- round(tmpDT$BAM_aligns / tmpDT$FASTQ_reads * 100, 2)
        tmpDT$pC_tx <- round(tmpDT$tx_starts / tmpDT$BAM_aligns * 100, 2)
        return(tmpDT)
    }else{
        stop("Some files are missing.")
    }
}

#' Alignment report plots
#'
#' @param rReport
#' @param species
#'
#' @return
#' @export
#'
#' @examples
gg_readStats <- function(rReport, species){
    tmpDT <- rReport[, -c(5, 6)]
    tmpDT$FASTQ <- tmpDT$FASTQ_reads - tmpDT$BAM_aligns
    tmpDT$BAM <- tmpDT$BAM_aligns - tmpDT$tx_starts
    tmpDT$txDT <- tmpDT$tx_starts
    tmpDT$BAM[is.na(tmpDT$BAM)] <- 0; tmpDT$txDT[is.na(tmpDT$txDT)] <- 0
    tmpDT <- tidyr::pivot_longer(data = tmpDT[,c("sample", "FASTQ", "BAM", "txDT")],
                                 cols = c("FASTQ", "BAM", "txDT"), names_to = "Reads")
    tmpDT$Reads <- factor(tmpDT$Reads, levels = c("FASTQ", "BAM", "txDT"))
    tmpDT$value <- magrittr::divide_by(tmpDT$value, 1e6)
    t_GG1 <- ggplot2::ggplot(tmpDT, ggplot2::aes(x = sample, y = value, fill = Reads)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Reads in", species, "libraries by processing step"),) +
        ggplot2::ylab("Million reads") +
        ggplot2::xlab("Samples")
    # Proportion
    tmpDT <- rReport[, -c(5, 6)]
    tmpDT$FASTQ <- tmpDT$FASTQ_reads - tmpDT$BAM_aligns
    tmpDT$BAM <- tmpDT$BAM_aligns - tmpDT$tx_starts
    tmpDT$txDT <- tmpDT$tx_starts
    tmpDT$BAM[is.na(tmpDT$BAM)] <- 0; tmpDT$txDT[is.na(tmpDT$txDT)] <- 0
    tmpDT <- data.frame(sample = tmpDT$sample,
                        apply(tmpDT[,c("FASTQ", "BAM", "txDT")], 1, prop.table) %>% t)
    tmpDT <- tidyr::pivot_longer(data = tmpDT[,c("sample", "FASTQ", "BAM", "txDT")],
                                 cols = c("FASTQ", "BAM", "txDT"), names_to = "Reads")
    tmpDT$Reads <- factor(tmpDT$Reads, levels = c("FASTQ", "BAM", "txDT"))

    t_GG2 <- ggplot2::ggplot(tmpDT, ggplot2::aes(x = sample, y = value, fill = Reads)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Proportion of reads in", species, "libraries by processing step"),) +
        ggplot2::ylab("Proportion") +
        ggplot2::xlab("Samples")
    return(list(t_GG1, t_GG2))
}

# Alignment and transcript data processing efficiency
ggAlignEffPlot <- function(META, rReport){
    tmp <- cbind(META, rReport[,-1])
    tmpGG1 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$libTreat, y = tmp$pC_BAM, colour = tmp$libTreat)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    tmpGG2 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$bioTreat, y = tmp$pC_BAM, colour = tmp$bioTreat)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    tmpGG3 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$RTase, y = tmp$pC_BAM, colour = tmp$RTase)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    list(tmpGG1, tmpGG2, tmpGG3)
}

#' Library complexity extrapolation plot
#'
#' @param META
#' @param tab_name
#' @param speciesName
#'
#' @return
#' @export
#'
#' @examples
gg_lce <- function(META, tab_name, speciesName = ""){
    lceFiles <- gsub(META$BAM, pattern = ".bam$", replacement = ".lce.txt") %>%
        magrittr::set_names(META$id)
    if(!all(file.exists(lceFiles))){
        stop("Report files missing:\n", paste(lceFiles[!file.exists(lceFiles)], collapse = " \n"))
    }
    tmp <- lapply(lceFiles, function(x) data.table::fread(x)) %>%
        magrittr::set_names(META$id)
    tmp <- lapply(names(tmp), function(x){
        cbind(tmp[[x]], id = x)
    }) %>% do.call(what = rbind) %>% data.table::data.table()
    data.table::fwrite(x = tmp, file = tab_name, sep = "\t")
    cat("Library complexity table output:", tab_name)
    tmpT <- table(tmp$TOTAL_READS)
    e_reads <- names(tmpT)[tmpT == max(tmpT)] %>% as.numeric %>% max
    tmp <- tmp[tmp$TOTAL_READS == e_reads,]
    ggOUT <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$id, y = tmp$EXPECTED_DISTINCT)) +
        ggplot2::geom_bar(stat="identity", color="black", position= ggplot2::position_dodge()) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin= tmp$LOWER_0.95CI, ymax= tmp$UPPER_0.95CI), width=.2,
                               position= ggplot2::position_dodge(1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(label = paste("Library complexity extrapolation -", speciesName),
                         sub = paste("Expected distinct reads at", e_reads, "depth. CI = 95%")) +
        ggplot2::ylab("Expected distinct reads") +
        ggplot2::xlab("Samples")
    return(ggOUT)
}

# is.Modified functions ##########################################################

# RNAmod logical vectors from character or factor
is.m1A <- function(nucs){nucs %in% c("m1A")}
is.m66A <- function(nucs){nucs %in% c("m66A")}
is.m3C <- function(nucs){nucs %in% c("m3C")}
is.m5C <- function(nucs){nucs %in% c("m5C")}
is.ac4C <- function(nucs){nucs %in% c("ac4C")}
is.m1G <- function(nucs){nucs %in% c("m1G")}
is.m22G <- function(nucs){nucs %in% c("m22G")}
is.m7G <- function(nucs){nucs %in% c("m7G")}
is.D <- function(nucs){nucs %in% c("D")}
is.pseudoU <- function(nucs){nucs %in% c("Y", "Ym", "m3Y")}
is.m3U <- function(nucs){nucs %in% c("m3U", "m3Y")}
is.m1acp3Y <- function(nucs){nucs %in% c("m1acp3Y")}
is.Am <- function(nucs){nucs %in% c("Am", "m1Am", "m66Am")}
is.Cm <- function(nucs){nucs %in% c("Cm", "ac4Cm", "m3Cm", "m5Cm")}
is.Gm <- function(nucs){nucs %in% c("Gm", "m1Gm", "m22Gm", "m7Gm")}
is.Um <- function(nucs){nucs %in% c("Um", "Um", "Dm", "Ym", "m3Um", "m1acp3Ym")}

# others
is.Nm <- function(nucs){nucs %in% c("Am", "Gm", "Um", "Cm", "Um", paste0(utils::head(RNAmods_vec, -4), "m"))} # Any 2-O-methylation
is.m5U <- function(nucs){nucs %in% c("m5U")} # not detectable
is.m6A <- function(nucs){nucs %in% c("m6A")} # not detectable
is.m2A <- function(nucs){nucs %in% c("m2A")} # not detectable
is.ho5C <- function(nucs){nucs %in% c("ho5C")}
is.m2G <- function(nucs){nucs %in% c("m2G")}
is.m3Y <- function(nucs){nucs %in% c("m3Y")}
is.m4Cm <- function(nucs){nucs %in% c("m4Cm")}

# Modeling functions #####

#  Balanced groups assignment of two level vectors for K-fold cross validation
CV_balancedGroups <- function(x, k){
    x <- factor(x)
    blocks <- seq(1, k)
    groups <- rep(NA, length(x))
    tmpLog1 <- x == levels(x)[1]
    tmpLog2 <- x == levels(x)[2]
    groups[tmpLog1] <- sample(rep(blocks, ceiling(sum(tmpLog1)/length(blocks))))[seq(groups[tmpLog1])]
    groups[tmpLog2] <- sample(rep(blocks, ceiling(sum(tmpLog2)/length(blocks))))[seq(groups[tmpLog2])]
    groups
}

# Extract metrics from k-CFV res
extract_AUC_sen_spe <- function(res, RNAmod, strategy, organism){
    data.frame(organism = organism,
               strategy = strategy,
               RNAmod = RNAmod,
               AUC = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "AUC"])),
               sensitivity = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "Sensitivity"])),
               specificity = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "Specificity"])))
}

# Matthew Correlation Coefficient
matthewCorrCoeff <- function(scores, successes){
    MCC <- NULL
    known <- rep(0, length(scores)); known[successes] <- 1
    NAs <- which(is.na(scores))
    scores[is.na(scores)] <- 0
    for (i in 1:length(scores)){
        thr <- scores[i]
        calls <- as.numeric(scores >= thr)
        TP <-  sum(calls & known)
        FP <-  sum(calls) - TP
        TN <- sum(!calls & !known)
        FN <- sum(!calls) - TN
        MCC[i] <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    }
    MCC[is.nan(MCC)] <- 0
    MCC[is.infinite(MCC)] <- 0
    MCC[NAs] <- NA
    return(MCC)
}

#' Matthews Correlation COefficient
#'
#' @param scores numeric
#' @param successes logical
#' @param thr numeric
#'
#' @return numeric
MCC <- function(scores, successes, thr){
    successes <- successes[!is.na(scores)]
    scores <- scores[!is.na(scores)]
    calls <- as.numeric(scores >= thr)
    known <- rep(0, length(scores)); known[successes] <- 1
    TP <-  sum(calls & known)
    FP <-  sum(calls) - TP
    TN <- sum(!calls & !known)
    FN <- sum(!calls) - TN
    MCC <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if(is.nan(MCC)){MCC <- 0}
    MCC
}

#' Maximum MCC
#'
#' @param scores numeric
#' @param successes logical
#' @param ncores integer
#'
#' @return numeric
#' @export
maxMCC <- function(scores, successes, ncores = 1L){
    THRS <- sort(unique(scores))
    MCCs <- parallel::mclapply(mc.cores = ncores, THRS, function(thr){
        MCC(scores, successes, thr)
    }) %>% do.call(what = "c")
    max(MCCs) %>% magrittr::set_names(THRS[which.max(MCCs)])
}

# All combinations of strings in character vector
all_comb_allNelements <- function(x){
    lapply(seq_along(x), function(y) utils::combn(x, y, simplify = FALSE)) %>%
        unlist(recursive = FALSE)
}

#' Logistic Model training using PANORAMA$CALLS object.
#'
#' @param PANORAMA
#' @param isModFun
#' @param RNAmod
#' @param accNArate logical. Account for NA rate for model selection
#' @param thresholds
#' @param minMCC
#'
#' @return
#' @export
trainLogModel_RNAmods <- function(PANORAMA, isModFun, RNAmod, thresholds,
                                  accNArate = TRUE, minMCC = 0.7){
    CALLS <- lapply(PANORAMA$CALLS, function(x) data.frame(x[[RNAmod]])) %>% do.call(what = "rbind")
    CALLS$isMod <- isModFun(CALLS$nuc)
    varNames <- names(CALLS)[names(CALLS) %in% allMetrics]
    allCombVar <- all_comb_allNelements(varNames)
    combNames <- paste0("Comb_", seq_along(allCombVar))
    OUT <- lapply(allCombVar, function(selVars){
        tData <- CALLS[CALLS$refSeq %in% RNAMod_baseNuc[[RNAmod]], c(selVars, "isMod")]
        tData$isMod <- as.numeric(tData$isMod)
        logMod <- stats::glm(isMod ~ ., family = stats::binomial(), data = tData) %>% suppressWarnings()
        CALLS$logPred <- stats::predict(logMod, newdata = CALLS, type = "response")
        CALLS <- CALLS[CALLS$refSeq %in% RNAMod_baseNuc[[RNAmod]],]
        return(list(MCCs = sapply(thresholds, function(x) MCC(CALLS$logPred, CALLS$isMod, x)),
                    NArate = mean(is.na(CALLS$logPred))))
    })
    MCCmat <- lapply(OUT, function(x) x$MCCs) %>% do.call(what = cbind) %>%
        magrittr::set_rownames(thresholds) %>% magrittr::set_colnames(combNames)
    NArate <- lapply(OUT, function(x) x$NArate) %>% unlist()
    medianMCC <- apply(MCCmat, 2, "median")
    maxMCC <-  apply(MCCmat, 2, "max")
    if(accNArate){
        medianMCC <- medianMCC
        maxMCC <- maxMCC
        modScore <- ((1 - NArate) * medianMCC)
        modScore2 <- ((1 - NArate) * maxMCC)
        selComb <- which.max(modScore * modScore2)
        selThr <- which.max(MCCmat[,selComb])
    }else{
        selComb <- which.max(medianMCC)
        selThr <- which.max(MCCmat[,selComb])
    }
    if(MCCmat[selThr, selComb] < minMCC){
        warning("Model for RNAmod:", RNAmod, " was set to NULL as obtained MCC:",
                round(MCCmat[selThr, selComb], 4), " is less than minimum MCC:",
                minMCC)
        return(list(logiMod = NULL, thr = thresholds[selThr], vars = NULL,
                    MCCmat = MCCmat, MCC = NULL, NArate = NULL,
                    onNucs = RNAMod_baseNuc[[RNAmod]], medianMCC = medianMCC,
                    maxMCC = maxMCC, NArates = NArate, selComb = NULL))
    }
    # FinalModel
    tData <- CALLS[CALLS$refSeq %in% RNAMod_baseNuc[[RNAmod]], c(allCombVar[[selComb]], "isMod")]
    tData$isMod <- as.numeric(tData$isMod)
    logMod <- stats::glm(isMod ~ ., family = stats::binomial(), data = tData) %>% suppressWarnings()
    return(list(logiMod = logMod, thr = thresholds[selThr], vars = allCombVar[[selComb]],
                MCCmat = MCCmat, MCC = MCCmat[selThr, selComb], NArate = NArate[selComb],
                onNucs = RNAMod_baseNuc[[RNAmod]], medianMCC = medianMCC, maxMCC = maxMCC,
                NArates = NArate, selComb = selComb))
}

#' Train all logistic models
#'
#' @param PANORAMA list. PANORAMA object
#' @param thrRange numeric vector. set of thresholds used to search for a maximum
#' MatthewsCorrelationCoefficient
#' @param accNArate logical. Account for NA rate for model selection
#' @param minMCC numeric. Minimum MCC needed for RNAmod model to be output.
#' Set to -1 to force all models to be output.
#'
#' @return
#' @export
trainLogModel_allRNAmods <- function(PANORAMA, thrRange = seq(-0.1, 0.5, 0.01),
                                     accNArate = TRUE, minMCC = 0.7){
    RNAmods <- check_whichRNAmods(PANORAMA)
    lapply(RNAmods, function(RNAmod_i){
        trainLogModel_RNAmods(PANORAMA = PANORAMA, isModFun = RNAModFunList[[RNAmod_i]],
                              RNAmod = RNAmod_i, thresholds = thrRange,
                              accNArate = accNArate, minMCC = minMCC)
    }) %>% magrittr::set_names(RNAmods)
}


#' Logistic regressions scores
#'
#' Update PANORAMA$CALLS with logistic and linear scores and prediction
#'
#' @param PANORAMA
#' @param logMods
#'
#' @return
#' @export
call_RNAmods_logRegMods <- function(PANORAMA, logMods, nucBias = TRUE){
    PANORAMA$CALLS <- lapply(names(PANORAMA$CALLS), function(set){
        RNAmods <- intersect(names(PANORAMA$CALLS[[set]]), names(logMods))
        tmp <- lapply(RNAmods, function(RNAmod){
            if(is.null(logMods[[RNAmod]]$logiMod)){
                warning("There is no model for RNAmod:", RNAmod,
                        " predictions and logistic",
                        " scores, were not calculated.")
                return(PANORAMA$CALLS[[set]][[RNAmod]])
            }
            if(all(logMods[[RNAmod]]$vars %in% colnames(PANORAMA$CALLS[[set]][[RNAmod]]))){
                nData <- data.frame(PANORAMA$CALLS[[set]][[RNAmod]])
                PANORAMA$CALLS[[set]][[RNAmod]]$logist_Score <-
                    stats::predict.glm(logMods[[RNAmod]]$logiMod, newdata = nData, type = "response")
                PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score <-
                    stats::predict.glm(logMods[[RNAmod]]$logiMod, newdata = nData, type = "link")
                PANORAMA$CALLS[[set]][[RNAmod]]$pred <-
                    PANORAMA$CALLS[[set]][[RNAmod]]$logist_Score >= logMods[[RNAmod]]$thr
                if(nucBias){
                    notBaseNuc <- !(PANORAMA$CALLS[[set]][[RNAmod]]$refSeq %in% logMods[[RNAmod]]$onNucs)
                    PANORAMA$CALLS[[set]][[RNAmod]]$logist_Score[notBaseNuc] <- NA
                    PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score[notBaseNuc] <- NA
                    PANORAMA$CALLS[[set]][[RNAmod]]$pred[notBaseNuc] <- NA
                }
                PANORAMA$CALLS[[set]][[RNAmod]]$pred[is.na(PANORAMA$CALLS[[set]][[RNAmod]]$logist_Score)] <- NA
                return(PANORAMA$CALLS[[set]][[RNAmod]])
            }else{
                warning("Some model metrics are not found in PANORAMA$CALLS ",
                        "for RNAmod:", RNAmod, " predictions and logistic",
                        " scores, were not calculated.")
                return(PANORAMA$CALLS[[set]][[RNAmod]])
            }
            return(tmp)
        }) %>% magrittr::set_names(RNAmods)
    }) %>% magrittr::set_names(names(PANORAMA$CALLS))
    PANORAMA
}

# Save as RDS if object is not already RDS
saveAsRDSIfNotAlready <- function(object, fileName){
    if(!file.exists(fileName)){saveRDS(object, fileName)}
}

# Has duplicates?
hasDups <- function(x){
    sum(duplicated(x)) > 0
}

# Limit PANORAMA object to genes in geneAnnot argument
panorama_reduceToGA <- function(PANORAMA, geneAnnot){
    PANORAMA$DATA <- lapply(PANORAMA$DATA, function(DT) DT[gene %in% geneAnnot$name,])
    return(PANORAMA)
}

# panorama_calls: Makes predictions based on RF and CutPointR for each modification
panorama_calls <- function(PANORAMA, RNAmodList, RF_list, CP_models){
    RNAMODS <- names(RNAmodList)
    coorSys <- PANORAMA$RES[,names(PANORAMA$RES) %in% c("chr", "gencoor", "strand",
                                                  "gene", "txcoor", "pos",
                                                  "refSeq", "nuc"), with = FALSE]
    iSets <- levels(PANORAMA$META$set)
    PANORAMA$CALLS <- lapply(seq_along(RNAMODS), function(i){
        pred_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(PANORAMA$RES)) &
                    grepl(pattern = iSets[j], names(PANORAMA$RES))
                tmpDat <- data.frame(x = PANORAMA$RES[[which(selVar)]])
                CPmod$predictor <- "x"
                as.logical(stats::predict(CPmod, tmpDat))
            }) %>% do.call(what = data.frame) %>% unname %>% apply(1, all)
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>%
            magrittr::set_names(paste("pred", RNAMODS[i], iSets, sep = "_"))

        scor_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(PANORAMA$RES)) &
                    grepl(pattern = iSets[j], names(PANORAMA$RES))
                PANORAMA$RES[,which(selVar), with = FALSE]
            }) %>% do.call(what = cbind)
        }) %>% do.call(what = cbind) %>% data.table::data.table()

        tree_RES <- lapply(seq_along(iSets), function(j){
            RF_mod <- RF_list[[i]]
            RF_vars <- names(RF_mod$forest$xlevels)
            newDat <- lapply(RF_vars, function(RF_v){
                tmpP <- gsub(pattern = "_TGIRT", replacement = "", x = RF_v)
                selVar <- grepl(pattern = tmpP, x = names(PANORAMA$RES))
                selSmp <- grepl(pattern = iSets[j], x = names(PANORAMA$RES)) & selVar
                newDat <- data.frame(x = PANORAMA$RES[, selSmp, with = FALSE]) %>% magrittr::set_names(RF_v)
            }) %>% do.call(what = cbind)
            newDat[is.na(newDat)] <- 0
            out <- stats::predict(RF_mod, newDat, norm.votes = TRUE, type = "vote")
            unname(out[,"TRUE"])
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>%
            magrittr::set_names(paste("votes", RNAMODS[i], iSets, sep = "_"))

        data.table::data.table(coorSys, pred_RES, tree_RES, scor_RES)
    }) %>% magrittr::set_names(RNAMODS)
    return(PANORAMA)
}

#' Restore txDT Genomic Coordinate System
#'
#' @param txDT data.table
#' @param geneAnnot GenomicRanges Gene annotation loaded via tx_load_bed()
#' @param nCores integer
#'
#' @return data.table
#' @export
restoreGenCoors <- function (txDT, geneAnnot, nCores = 1){
    txtools:::check_mc_windows(nCores)
    if(all(txDT$gene %in% geneAnnot$name)){
        txLengths <- txtools::tx_get_geneLengths(txDT)
        genCoorSys <- parallel::mclapply(
            mc.cores = nCores,
            as.character(unique(txDT$gene)),
            function(iGene){
                tmp2 <- geneAnnot[which(geneAnnot$name == iGene)]
                tmp3 <- c(GenomicAlignments::seqnames(tmp2),
                          GenomicRanges::strand(tmp2)) %>% as.character() %>%
                    c(iGene)
                tmpDT <- rep(tmp3, txLengths[iGene]) %>%
                    matrix(ncol = 3, byrow = T) %>%
                    cbind(txtools:::exonBlockGen(iGene, geneAnnot)) %>%
                    cbind(seq(1, txLengths[iGene]))
                tmpDT <- tmpDT[, c(1, 4, 2, 3, 5)] %>% data.table::data.table() %>%
                    magrittr::set_colnames(c("chr", "gencoor",
                                             "strand", "gene", "txcoor"))
                tmpDT$chr <- as.factor(tmpDT$chr)
                tmpDT$gencoor <- as.integer(tmpDT$gencoor)
                tmpDT$strand <- as.factor(tmpDT$strand)
                tmpDT$gene <- as.factor(tmpDT$gene)
                tmpDT$txcoor <- as.integer(tmpDT$txcoor)
                return(tmpDT)
            }) %>% do.call(what = "rbind")
        txDT[, 1:5] <- genCoorSys
        txDT
    }
    else {
        stop("Genes of txDT are not contained in geneAnnot$name .\n")
    }
}

# Linear models

#' Train all linear models
#'
#' @param PANORAMA list. PANORAMA object
#' @param nCores integer. Number of cores to use
#' @param n_thresholds integer. Number of thresholds the range from the median value
#' of non-modified nucleotides to the median value of known modified nucleotides
#' is divided to test for ideal linear score threshold.
#' @param min_AUC numeric. Minimum AUC needed for a metric to be considered into
#' the linear integrated score
#' @param blacklist data.frame. Table that includes the genes and txcoor to be
#' ignored from the known RNAmod nucleotides reference
#' @param rmGeneNoRNAmod logical. Remove genes with not RNAmods
#' @param removeOtherRNAmods logical. Remove other RNAmods from training step
#' @param nucBias logical. Only consider the values measured in the base-nucleotides
#' relevant to the RNAmod for training the model, e.g. only adenines for m1A
#' @param min_MCC numeric Minimum MCC for model to be output. Default = 0.7
#' @param rmOtherRNAmodKeepConf logical. Remove other RNAmods keeping confounding
#' RNAmods. Defined by `confRNAmods_list`.
#' @param rmRNAmod_thrSel
#' @param rmUnmodTRNA
#' @param max_NArate
#' @param geneInfo
#'
#' @return list. Linear models and diagnostic plots for each trained RNAmod
#' @export
trainLinModel_allRNAmods <- function(PANORAMA,
                                     nCores = 6,
                                     n_thresholds = 100,
                                     min_AUC = 0.9,
                                     min_MCC = 0.5,
                                     blacklist = NULL,
                                     rmGeneNoRNAmod = FALSE,
                                     rmRNAmod_thrSel = TRUE,
                                     removeOtherRNAmods = FALSE,
                                     rmUnmodTRNA = TRUE,
                                     rmOtherRNAmodKeepConf = TRUE,
                                     nucBias = TRUE,
                                     max_NArate = 0.4,
                                     geneInfo = NULL){
    RNAmods <- check_whichRNAmods(PANORAMA)
    parallel::mclapply(mc.cores = nCores, RNAmods, function(RNAmod_i){
        trainLinModel_RNAmods(PANORAMA = PANORAMA,
                              RNAmod = RNAmod_i,
                              n_thresholds = n_thresholds,
                              min_AUC = min_AUC,
                              min_MCC = min_MCC,
                              blacklist = blacklist,
                              rmGeneNoRNAmod = rmGeneNoRNAmod,
                              rmRNAmod_thrSel = rmRNAmod_thrSel,
                              removeOtherRNAmods = removeOtherRNAmods,
                              rmUnmodTRNA = rmUnmodTRNA,
                              rmOtherRNAmodKeepConf = rmOtherRNAmodKeepConf,
                              nucBias = nucBias,
                              max_NArate = max_NArate,
                              geneInfo = geneInfo)
    }) %>% magrittr::set_names(RNAmods)
}

#' Training linear model
#'
#' @param PANORAMA
#' @param RNAmod
#' @param n_thresholds
#' @param min_AUC
#' @param min_MCC
#' @param blacklist
#' @param rmGeneNoRNAmod
#' @param rmRNAmod_thrSel logical. Removes annotated modified nucleotides for
#' threshold selection
#' @param removeOtherRNAmods logical. Removes all other RNAmods
#' @param rmUnmodTRNA
#' @param rmOtherRNAmodKeepConf logical. Removes all other RNAmods but keeps confounding RNAmods
#' @param nucBias logical. Removes all positions for which the refSeq is not the same as the reference nucleotide
#' @param max_NArate logical. Maximum allowed NA rate for metric in modified sites
#' @param geneInfo
#'
#' @return list
#' @export
trainLinModel_RNAmods <- function(PANORAMA,
                                  RNAmod,
                                  n_thresholds = 100,
                                  min_AUC = 0.9,
                                  min_MCC = 0.5,
                                  blacklist = NULL,
                                  rmGeneNoRNAmod = FALSE,
                                  rmRNAmod_thrSel = TRUE,
                                  removeOtherRNAmods = FALSE,
                                  rmUnmodTRNA = TRUE,
                                  rmOtherRNAmodKeepConf = TRUE,
                                  nucBias = TRUE,
                                  max_NArate = 0.4,
                                  geneInfo = NULL){
    emptyModel <- list(list(model = NULL, model_type = NULL, thr = NULL, highThr = NULL,
                            lmScore_thrs = NULL, vars = NULL, MCCs = NULL,
                            AUCs = NULL, MCC = NULL, onNucs = NULL, plot_score = NULL,
                            plot_pairs = NULL))
    isModFun <- RNAModFunList[[RNAmod]]
    unModNuc <- RNAMod_nucRef[[RNAmod]]
    confRNAmods <- confRNAmods_list[[RNAmod]]
    baseNuc <- RNAMod_baseNuc[[RNAmod]]
    CALLS <- lapply(PANORAMA$CALLS, function(x) data.frame(x[[RNAmod]])) %>% do.call(what = "rbind")
    # Sites filtering
    if(rmUnmodTRNA){
        if(is.null(geneInfo)){stop("If unmodified tRNA positions are to be removed",
                                   ", geneInfo argument must have a table with gene's 'name' and 'biotype' columns")}
        CALLS <- add_biotype(CALLS, geneInfo)
        CALLS <- dplyr::filter(CALLS, !(nuc == unModNuc & gene_biotype == "tRNA"))
    }
    if(sum(isModFun(CALLS$nuc)) == 0){
        warning("Model for RNAmod: ", RNAmod, " was set to NULL as no known sites are present in dataset.")
        return(emptyModel)
    }
    if(!is.null(blacklist)){
        blacklist$pos <- paste(blacklist$gene, blacklist$txcoor, sep = ":")
        CALLS <- CALLS[!CALLS$pos %in% blacklist$pos,]
    }
    if(removeOtherRNAmods){ # Remove other RNA modifications
        CALLS$nuc <- as.character(CALLS$nuc)
        CALLS$nuc[CALLS$nuc == "Um"] <- "Um"
        CALLS <- CALLS[CALLS$nuc %in% c(unModNuc, RNAmod),]
        CALLS$nuc[CALLS$nuc == "Um"] <- "Um"
        CALLS$nuc <- as.factor(CALLS$nuc)
    }else if(rmOtherRNAmodKeepConf){
        CALLS <- dplyr::filter(CALLS, nuc %in% c(unModNuc, RNAmod, confRNAmods))
    }
    if(nucBias){
        CALLS <- dplyr::filter(CALLS, refSeq == baseNuc)
    }
    if(rmGeneNoRNAmod){
        if(rmOtherRNAmodKeepConf){
            keepGenes <- unique(as.character(CALLS$gene)[isModFun(CALLS$nuc) | CALLS$nuc %in% confRNAmods])
        }else{
            keepGenes <- unique(as.character(CALLS$gene)[isModFun(CALLS$nuc)])
        }
        CALLS <- CALLS[CALLS$gene %in% keepGenes, ]
    }
    FEATS <- CALLS[, c("chr", "gencoor", "strand", "gene", "txcoor", "refSeq", "nuc")]
    varNames <- intersect(names(CALLS), allMetrics)
    FEATS$isMod <- isModFun(FEATS$nuc)
    FEATS$isConfRNAmod <- FEATS$nuc %in% confRNAmods
    FEATS$isBaseNuc <- FEATS$nuc %in% unModNuc
    isMod <- isModFun(FEATS$nuc)
    if(sum(isMod) == 0){
        warning("Model for RNAmod: ", RNAmod, " was set to NULL as no known sites are present in data after filtering.")
        return(emptyModel)
    }
    # Variables filtering
    # Informative variables
    res_auc <- purrr::map_dbl(varNames, function(var_i){
        values <- CALLS[[var_i]]
        isModV <- FEATS$isMod %>% as.numeric() %>% factor()
        sel_v <- !is.na(values)
        values <- values[sel_v]
        isModV <- isModV[sel_v]
        AUC::auc(AUC::roc(values, isModV))
    }) %>% magrittr::set_names(varNames)
    maxMCCs <- purrr::map_dbl(CALLS[, varNames], function(x){
        maxMCC(x, isMod)}) %>% magrittr::set_names(varNames)
    infoVars <- varNames[res_auc >= min_AUC | res_auc <= abs(min_AUC - 1)]
    # Variables with NA rate less than max_NArate in modified sites
    modSites_NArate <- purrr::map_dbl(varNames, function(var_i){
        mean(is.na(dplyr::filter(CALLS, nuc == RNAmod)[[var_i]]))
    }) %>% set_names(varNames)
    passNArate_vars <- names(modSites_NArate)[modSites_NArate <= max_NArate]
    passVars <- intersect(passNArate_vars, infoVars)
    if(length(passVars) == 0){
        warning("Model for: ", RNAmod, " was not generated due to no metric passsing informative and maxNArate filters")
        return(emptyModel)
    }
    CALLS <- dplyr::select(CALLS, all_of(passVars))

    # Keep rows with no NA values
    selR <- rowSums(is.na(CALLS)) == 0
    CALLS <- CALLS[selR,]
    FEATS <- FEATS[selR,]

    # Training model
    trainData <- data.table::data.table(CALLS, isMod = as.numeric(FEATS$isMod)) %>%
        magrittr::set_names(c(passVars, "isMod"))
    model_lin <- stats::lm(isMod ~ ., data = trainData)
    model_lin$posData <- FEATS
    trainData$linear_Score <- stats::predict(model_lin, newdata = trainData)
    thresholds <- seq(from = stats::median(trainData$linear_Score[FEATS$isBaseNuc], na.rm = TRUE),
                      to = stats::median(trainData$linear_Score[FEATS$isMod], na.rm = TRUE),
                      length.out = n_thresholds) %>% round(digits = 4)

    MCCs <- purrr::map_dbl(thresholds, function(thr_i){
        if(rmRNAmod_thrSel){
            MCC(trainData$linear_Score[FEATS$isBaseNuc | FEATS$isMod], FEATS$isMod[FEATS$isBaseNuc | FEATS$isMod], thr_i)
        }else{
            MCC(trainData$linear_Score, FEATS$isMod, thr_i)
        }
    })
    # Selecting threshold
    selThr <- thresholds[which.max(MCCs)] # low-confidence threshold
    if(rmRNAmod_thrSel){
        MCC_lm <- MCC(trainData$linear_Score[FEATS$isBaseNuc | FEATS$isMod], FEATS$isMod[FEATS$isBaseNuc | FEATS$isMod], selThr) %>%
            magrittr::set_names("linear_Score")
        AUC_lm <- myROC(trainData$linear_Score[FEATS$isBaseNuc | FEATS$isMod], FEATS$isMod[FEATS$isBaseNuc | FEATS$isMod]) %>%
            AUC::auc() %>% magrittr::set_names("linear_Score")
    }else{
        MCC_lm <- MCC(trainData$linear_Score, FEATS$isMod, selThr) %>% magrittr::set_names("linear_Score")
        AUC_lm <- myROC(trainData$linear_Score, FEATS$isMod) %>% AUC::auc() %>% magrittr::set_names("linear_Score")
    }
    medMod <- stats::median(trainData[trainData$isMod == 1, ]$linear_Score, na.rm = TRUE)
    medUnMod <- stats::median(trainData[FEATS$nuc == unModNuc, ]$linear_Score, na.rm = TRUE)
    highMod <- stats::quantile(trainData[trainData$isMod == 1, ]$linear_Score, 0.95, na.rm = TRUE)
    higThr <- hlp_calcMidpoint(selThr, medMod) # High-confidence threshold
    lmScore_thrs <- c(medUnMod, selThr, higThr, medMod, highMod) %>%
        magrittr::set_names(c("median_baseNuc", "lowConfidence_thr", "highconfidence_thr", "median_modified", "modified_95pct"))
    if(MCC_lm < min_MCC){
        warning("Model for RNAmod: ", RNAmod, " was set to NULL as obtained MCC: ",
                round(MCC_lm, 4), " is less than minimum MCC: ",
                min_MCC)
        return(emptyModel)
    }
    maxMCCs <- c(maxMCCs, MCC_lm); res_auc <- c(res_auc, AUC_lm)
    tmpDT <- data.table::data.table(linear_Score = trainData$linear_Score,
                                    nuc = factor(FEATS$nuc),
                                    gene = FEATS$gene) %>% stats::na.omit()
    labs1 <- c(RNAmod, unModNuc, confRNAmods,  "Low-confidence", "High-confidence")
    altColors <- c("purple", "limegreen", "gold1", "deeppink")
    cols <- c("dodgerblue", "gray", altColors[seq(confRNAmods)],  "orange", "red") %>% set_names(labs1)
    tmpGG <- ggplot2::ggplot(tmpDT, ggplot2::aes(y = .data$linear_Score, x = .data$gene, color = .data$nuc)) +
        ggplot2::geom_jitter(alpha = 0.5) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = selThr, color = "Low-confidence")) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = higThr, color = "High-confidence")) +
        ggplot2::scale_colour_manual(values = cols, aesthetics = c("color", "fill"), breaks = labs1) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") +
        ggplot2::ggtitle(RNAmod)

    tmpGGpairs <- GGally::ggpairs(data = trainData, mapping = ggplot2::aes(color = factor(ifelse(isMod, "modified", "background"))),
                                  columns = which(names(trainData) != "isMod")) +
        ggplot2::ggtitle(paste("Pairs plot -", RNAmod)) + ggplot2::theme_bw()
    # OUT Model
    return(list(model = model_lin, model_type = "lm", thr = selThr, highThr = higThr,
                lmScore_thrs = lmScore_thrs, vars = passVars, MCCs = maxMCCs,
                AUCs = res_auc, MCC = MCC_lm, onNucs = baseNuc, plot_score = tmpGG,
                plot_pairs = tmpGGpairs))
}

#' Models summary table
#'
#' Generates a table summarizing the models variables, coefficients and
#' correlation to the fitted results used in training.
#'
#' @param models
#'
#' @return data.table
#' @export
modelSummary <- function(models){
    modelVars <- lapply(names(models), function(RNAmod_i){
        tmpV <- models[[RNAmod_i]]$vars
        if(length(tmpV) == 0){return(NULL)}
        tmpC <- purrr::map_dbl(tmpV, function(var_i){
            stats::cor(models[[RNAmod_i]]$model$model[,var_i], models[[RNAmod_i]]$model$fitted.values, use = "p")
        })
        data.table::data.table(RNAmod = RNAmod_i,
                               metric = tmpV,
                               coeff = signif(models[[RNAmod_i]]$model$coefficients[tmpV], 4),
                               pearson_R = round(tmpC, 4))
    }) %>% do.call(what = rbind)
    modelVars$RNAmod <- factor(modelVars$RNAmod, levels = PANORAMA::RNAmods_vec)
    modelVars
}

#' Call RNA modifications using linear regression models
#'
#' @param PANORAMA
#' @param models
#' @param nucBias
#'
#' @return PANORAMA object
#' @export
call_RNAmods_linRegMods <- function(PANORAMA, models, nucBias = TRUE){
    PANORAMA$CALLS <- lapply(names(PANORAMA$CALLS), function(set){
        RNAmods <- intersect(names(PANORAMA$CALLS[[set]]), names(models))
        tmp <- lapply(RNAmods, function(RNAmod){
            if(is.null(models[[RNAmod]]$model)){
                warning("There is no model for RNAmod:", RNAmod,
                        " predictions and linear",
                        " score were not calculated.")
                return(PANORAMA$CALLS[[set]][[RNAmod]])
            }
            if(all(models[[RNAmod]]$vars %in% colnames(PANORAMA$CALLS[[set]][[RNAmod]]))){
                nData <- data.frame(PANORAMA$CALLS[[set]][[RNAmod]])
                PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score <-
                    stats::predict(models[[RNAmod]]$model, newdata = nData)
                PANORAMA$CALLS[[set]][[RNAmod]]$pred <-
                    PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score >= models[[RNAmod]]$thr
                # High/low confidence
                PANORAMA$CALLS[[set]][[RNAmod]]$confidence <- factor(NA, levels = c("low-confidence", "high-confidence"))
                PANORAMA$CALLS[[set]][[RNAmod]]$confidence[PANORAMA$CALLS[[set]][[RNAmod]]$pred] <- "low-confidence"
                PANORAMA$CALLS[[set]][[RNAmod]]$confidence[PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score >= models[[RNAmod]]$highThr] <- "high-confidence"
                if(nucBias){
                    notBaseNuc <- !(PANORAMA$CALLS[[set]][[RNAmod]]$refSeq %in% models[[RNAmod]]$onNucs)
                    PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score[notBaseNuc] <- NA
                    PANORAMA$CALLS[[set]][[RNAmod]]$pred[notBaseNuc] <- NA
                    PANORAMA$CALLS[[set]][[RNAmod]]$confidence[notBaseNuc] <- NA
                }
                PANORAMA$CALLS[[set]][[RNAmod]]$pred[is.na(PANORAMA$CALLS[[set]][[RNAmod]]$linear_Score)] <- NA
                return(PANORAMA$CALLS[[set]][[RNAmod]])
            }else{
                warning("Some model metrics are not found in PANORAMA$CALLS ",
                        "for RNAmod:", RNAmod, " predictions and linear",
                        " score were not calculated.")
                return(PANORAMA$CALLS[[set]][[RNAmod]])
            }
            return(tmp)
        }) %>% magrittr::set_names(RNAmods)
    }) %>% magrittr::set_names(names(PANORAMA$CALLS))
    PANORAMA
}


#' Subset PANORAMA object by genes
#'
#' @param PANORAMA list
#' @param genes character
#'
#' @return list. PANORAMA object
#' @export
subsetPANORAMA_byGenes <- function(PANORAMA, genes){
    if(!is.null(PANORAMA$DATA)){
        SAMPLES <- names(PANORAMA$DATA)
        PANORAMA$DATA <- lapply(SAMPLES, function(samp){
            base::subset(PANORAMA$DATA[[samp]], gene %in% genes)
        }) %>% magrittr::set_names(SAMPLES)
    }
    if(!is.null(PANORAMA$RES)){
        PANORAMA$RES <- base::subset(PANORAMA$RES, gene %in% genes)
    }
    if(!is.null(PANORAMA$CALLS)){
        SETS <- names(PANORAMA$CALLS)
        PANORAMA$CALLS <- lapply(SETS, function(set_i){
            RNA_MODS <- names(PANORAMA$CALLS[[set_i]])
            lapply(RNA_MODS, function(RNAMOD_i){
                base::subset(PANORAMA$CALLS[[set_i]][[RNAMOD_i]], gene %in% genes)
            }) %>% magrittr::set_names(RNA_MODS)
        }) %>% magrittr::set_names(SETS)
    }
    PANORAMA
}

#' Confusion matrix
#'
#' Confusion matrix generation for known sites in PANORAMA object. Known sites are
#' added with
#'
#' @param PANORAMA
#' @param nucBias
#' @param removeOtherRNAmods
#' @param rmNotModGenes
#' @param blacklist
#'
#' @return data.table
#' @export
panorama_confMatrix <- function(PANORAMA, nucBias = TRUE, removeOtherRNAmods = TRUE,
                             rmNotModGenes = TRUE, blacklist = NULL){
    if(!"CALLS" %in% names(PANORAMA)){
        stop("CALLS is not an element of PANORAMA, use panorama_assignScores() to ",
             "assign metrics to the")
    }
    lapply(names(PANORAMA$CALLS), function(set_i){
        lapply(names(PANORAMA$CALLS[[set_i]]), function(RNAmod_i){
            RNAModFun <- RNAModFunList[[RNAmod_i]]
            selCALLS <- PANORAMA$CALLS[[set_i]][[RNAmod_i]]
            if(!is.null(blacklist)){
                blacklist$pos <- paste(blacklist$gene, blacklist$txcoor, sep = ":")
                selCALLS <- selCALLS[!selCALLS$pos %in% blacklist$pos,]
            }
            if(nucBias){
                selCALLS <- selCALLS[selCALLS$refSeq %in% RNAMod_baseNuc[[RNAmod_i]],]
            }
            if(rmNotModGenes){
                modGenes <- unique(selCALLS[RNAModFun(selCALLS$nuc),]$gene)
                selCALLS <- dplyr::filter(selCALLS, .data$gene %in% modGenes)
            }
            if(removeOtherRNAmods){
                selCALLS$nuc <- as.character(selCALLS$nuc)
                selCALLS$nuc[selCALLS$nuc == "Um"] <- "Um"
                selCALLS <- selCALLS[selCALLS$nuc %in% c(RNAMod_nucRef[[RNAmod_i]], RNAmod_i),]
                selCALLS$nuc[selCALLS$nuc == "Um"] <- "Um"
                selCALLS$nuc <- as.factor(selCALLS$nuc)
            }
            # Empty table in case there are no predictions or no modified sites
            if(sum(RNAModFun(selCALLS$nuc), na.rm = TRUE) < 1 | sum(selCALLS$pred, na.rm = TRUE) < 1){
                return(data.table::data.table(outcome = NA, refSeq = NA,
                                              RNAmod = RNAmod_i, set = set_i, freq = NA))
            }
            tmpT <- table(truth = factor(RNAModFun(selCALLS$nuc), levels = c("TRUE", "FALSE")),
                          pred = factor(selCALLS$pred, levels = c("TRUE", "FALSE")),
                          refSeq = selCALLS$refSeq, useNA = "always") %>%
                data.frame() %>% data.table::data.table() %>% hlpr_confusMatOutcome()
            return(data.table::data.table(outcome = tmpT$outcome,
                                          refSeq = tmpT$refSeq,
                                          RNAmod = RNAmod_i,
                                          set = set_i,
                                          freq = tmpT$Freq))
        }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)
}

panorama_confMatrix4plot <- function(PANORAMA, removeOtherRNAmods = TRUE,
                                  rmNotModGenes = TRUE, blacklist = NULL){
    lapply(names(PANORAMA$CALLS), function(set_i){
        lapply(names(PANORAMA$CALLS[[set_i]]), function(RNAmod_i){
            RNAModFun <- RNAModFunList[[RNAmod_i]]
            selCALLS <- PANORAMA$CALLS[[set_i]][[RNAmod_i]]
            if(!is.null(blacklist)){
                blacklist$pos <- paste(blacklist$gene, blacklist$txcoor, sep = ":")
                selCALLS <- selCALLS[!selCALLS$pos %in% blacklist$pos,]
            }
            if(rmNotModGenes){
                modGenes <- unique(selCALLS[RNAModFun(selCALLS$nuc),]$gene)
                selCALLS <- dplyr::filter(selCALLS, gene %in% modGenes)
            }
            if(removeOtherRNAmods){
                selCALLS$nuc <- as.character(selCALLS$nuc)
                selCALLS$nuc[selCALLS$nuc == "Um"] <- "Um"
                selCALLS <- selCALLS[selCALLS$nuc %in% c(RNAMod_nucRef[[RNAmod_i]], RNAmod_i),]
                selCALLS$nuc[selCALLS$nuc == "Um"] <- "Um"
                selCALLS$nuc <- as.factor(selCALLS$nuc)
            }
            # Empty table in case there are no predictions or no modified sites
            if(sum(RNAModFun(selCALLS$nuc), na.rm = TRUE) < 1 | sum(selCALLS$pred, na.rm = TRUE) < 1){
                return(NULL)
                # tmpT <- table(truth = c(0, 0, 0, 0), pred = c(0, 0, 0, 0), refSeq = c("A", "T", "G", "C")) %>%
                #     data.frame() %>% data.table::data.table() %>%
                #     hlpr_confusMatOutcome()
            }else{
                tmpT <- table(truth = selCALLS$nuc %>% RNAModFun(),
                              pred = selCALLS$pred, refSeq = selCALLS$refSeq) %>%
                    data.frame() %>% data.table::data.table() %>%
                    hlpr_confusMatOutcome()
            }
            data.table::data.table(outcome = rep(tmpT$outcome, tmpT$Freq),
                                   refSeq = rep(tmpT$refSeq, tmpT$Freq),
                                   RNAmod = RNAmod_i,
                                   set = set_i)
        }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)
}

#' Evaluate predictions
#'
#' Evaluate predictions using known modified nucleotides identity
#'
#' @param PANORAMA list. PANORAMA object
#' @param nucBias logical. Nucleotide bias: only cases in the base nucleotide
#' of the respective RNAmod are considered.
#' @param roundDig numeric. Number of digits to round the resulting specificity,
#' sensitivity, FDR, and query rate values.
#' @param removeOtherRNAmods logical. Remove other RNAmods from training step
#'
#' @return
#' @export
panorama_evalPreds <- function(PANORAMA, nucBias = TRUE, roundDig = 4,
                            removeOtherRNAmods = TRUE, rmNotModGenes = TRUE,
                            blacklist = NULL, summarizeSETS = TRUE){
    SETS <- unique(PANORAMA$META$set)
    confMat <- panorama_confMatrix(PANORAMA, nucBias = nucBias,
                                removeOtherRNAmods = removeOtherRNAmods,
                                rmNotModGenes = rmNotModGenes, blacklist = blacklist)
    if(summarizeSETS & (length(SETS) > 1)){
        sumName <- paste0("summary_", length(SETS), "_smps")
        confMat <- confMat %>% dplyr::group_by(outcome, refSeq, RNAmod) %>%
            dplyr::summarize(freq = sum(freq), .groups = "keep") %>%
            tibble::add_column(set = sumName, .after = "RNAmod") %>%
            data.table::data.table()
    }
    lapply(unique(confMat$set), function(set_i){
        RNAMODS <- intersect(RNAmods_vec, unique(confMat$RNAmod))
        lapply(RNAMODS, function(RNAmod_i){
            RNAModFun <- RNAModFunList[[RNAmod_i]]
            tmpC <- confMat[confMat$RNAmod == RNAmod_i & confMat$set == set_i,]
            tmpC <- tapply(tmpC$freq, tmpC$outcome, "sum")
            sens <- tmpC["TP"] / (tmpC["TP"] + tmpC["FN"])
            spec <- tmpC["TN"] / (tmpC["FP"] + tmpC["TN"])
            FDR <- tmpC["FP"] / (tmpC["FP"] + tmpC["TP"])
            nMod <- sum(RNAModFun(PANORAMA$CALLS[[set_i]][[RNAmod_i]]$nuc))
            nPred <- sum(tmpC["TP"], tmpC["TN"], tmpC["FP"], tmpC["FN"], na.rm =  TRUE)
            nTrial <- sum(tmpC["TP"], tmpC["TN"], tmpC["FP"], tmpC["FN"], tmpC["NA"], na.rm =  TRUE)
            data.frame(set = set_i,
                       RNAmod = RNAmod_i,
                       TP = tmpC["TP"],
                       FP = tmpC["FP"],
                       TN = tmpC["TN"],
                       FN = tmpC["FN"],
                       sens = round(sens, roundDig),
                       spec = round(spec, roundDig),
                       FDR = round(FDR, roundDig),
                       MCC = mcc(tp = tmpC["TP"], tn = tmpC["TN"], fp = tmpC["FP"], fn = tmpC["FN"]),
                       n_knownRNAmod = nMod,
                       queryRate = round(nPred / nTrial, roundDig)
            )
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind") %>% magrittr::set_rownames(NULL) %>% # TODO: Make output to be a PANORAMA object
        data.table::data.table()
}

#' Generate inspection table from PANORAMA object
#'
#' @param PANORAMA list. PANORAMA object
#' @param nCores numeric. Number of cores used for multi-thread processes
#'
#' @return list containing inspection table and sets vector
#' @export
panorama_insTable <- function(PANORAMA, nCores){
    PANORAMA <- rm_knownRNAmods(PANORAMA)
    minMETA <- miniMETA(PANORAMA)
    minMETA$id <- paste(minMETA$bioTreat, minMETA$replicate, sep = "_")
    # Summarize PANORAMA run results
    SETS <- minMETA$set
    RNAMODS <- intersect(RNAmods_vec, names(PANORAMA$CALLS[[1]]))
    METRICS <- unique(PANORAMA$RES$metric)

    # Putatively modified, any sample.
    tmp_txDT <- dplyr::select(PANORAMA$CALLS[[1]][[1]], all_of(panorama_baseCoorCols)) # initialize txDT
    bio_fct <- minMETA$bioTreat[match(SETS, minMETA$set)]

    RES_meaDet <- parallel::mclapply(mc.cores = nCores, RNAMODS, function(RNAmod_i){
        # lm scores
        if(is.null(PANORAMA$CALLS[[1]][[RNAmod_i]]$linear_Score)){return(NULL)}
        tmpData <- lapply(SETS, function(set_i){
            PANORAMA$CALLS[[set_i]][[RNAmod_i]]$linear_Score
        }) %>% do.call(what = cbind) %>% set_colnames(SETS)
        # summarize lm scores
        tmpV <- tmpData %>% round(4)
        tmpV[is.na(tmpV)] <- ""
        summVal <- apply(tmpV, MARGIN = 1, function(x){paste(x, collapse = ",")})
        maxVal <- apply(tmpData, MARGIN = 1, function(x) max(x, na.rm = TRUE)) %>%
            suppressWarnings()
        maxVal[is.infinite(maxVal)] <- NA
        # Putative sites = If detected in any sample
        tmpR <- lapply(SETS, function(set_i){
            PANORAMA$CALLS[[set_i]][[RNAmod_i]]$pred
        }) %>% do.call(what = cbind)
        det_times <- rowSums(tmpR, na.rm = TRUE)
        detected <- det_times > 0
        # Detected in WHICH bio-treatments
        tmpS <- lapply(seq(nrow(tmpR)), function(i){
            tapply(X = tmpR[i,], INDEX = bio_fct, FUN = any)
        }) %>% do.call(what = rbind) %>% data.table::data.table()
        det_bioT <- lapply(seq(ncol(tmpS)), function(j){
            tmp <- character(length(tmpS[[j]]))
            tmp[which(tmpS[[j]])] <- names(tmpS)[j]
            tmp
        }) %>% do.call(what = cbind) %>%
            apply(MARGIN = 1, function(x){paste(x, collapse = ",")})
        # Detected confidence
        tmpR <- lapply(SETS, function(set_i){
            as.character(PANORAMA$CALLS[[set_i]][[RNAmod_i]]$confidence)
        }) %>% do.call(what = cbind)
        tmpR[is.na(tmpR)] <- ""
        tmpR[tmpR == "high-confidence"] <- "h"
        tmpR[tmpR == "low-confidence"] <- "l"
        set_confidence <- apply(tmpR, MARGIN = 1, function(x){paste(x, collapse = ",")})
        max_confidence <- lapply(seq(nrow(tmpR)), function(i){
            if(any(tmpR[i,] == "h")){
                return("high-confidence")
            }else if(any(tmpR[i,] == "l")){
                return("low-confidence")
            }else {return("")}
        }) %>% unlist()
        # Dynamic sites: AnalysisOfVariance test
        put_index <- which(detected)
        res_aov <- as.numeric(rep(NA, nrow(tmpData)))
        whichDyn <- as.logical(rep(NA, nrow(tmpData)))
        det_direc <- character(length = nrow(tmpData))
        coherent_dir <- character(length = nrow(tmpData))
        if(length(put_index) > 0){
            res_aov[put_index] <- lapply(put_index, function(i){
                tmpDF <- data.frame(lmScore = tmpData[i, ], bioT = bio_fct)
                tmpAOV <- try(stats::aov(lmScore ~ bioT, data = tmpDF), silent = TRUE)
                if(class(tmpAOV)[1] == "try-error"){
                    return(NA)
                }else{
                    dummy <- summary(tmpAOV)
                    tmpOut <- dummy[[1]][["Pr(>F)"]][1]
                    if(is.null(tmpOut)){
                        return(NA)
                    }else{
                        return(tmpOut)
                    }
                }
            }) %>% do.call(what = c)
            whichDyn <- res_aov < 0.05 # TRUE if p-Val < 0.05
            bioTreat_mean <- lapply(put_index, function(i){
                tapply(X = tmpData[i,], INDEX = bio_fct,
                       FUN = function(x) mean(x, na.rm = TRUE)) %>%
                    diff() %>% is_greater_than(0) %>% ifelse("+", "-")
            }) %>% do.call(what = rbind) %>% data.frame()
            det_direc[put_index] <- apply(bioTreat_mean, MARGIN = 1,
                                          function(x){paste(x, collapse = ",")})
            selR_pos <- apply(bioTreat_mean, MARGIN = 1, FUN = function(x) all(x == "+"))
            selR_neg <- apply(bioTreat_mean, MARGIN = 1, FUN = function(x) all(x == "-"))
            coherent_dir[put_index][selR_pos] <- "positive"
            coherent_dir[put_index][selR_neg] <- "negative"
        }
        IT_RNAmodColNames <- c(paste("lmScore", RNAmod_i, sep = "_"),
                               paste("lmScore_max", RNAmod_i, sep = "_"),
                               paste("detected", RNAmod_i, sep = "_"),
                               paste("detected_times", RNAmod_i, sep = "_"),
                               paste("confidence", RNAmod_i, sep = "_"),
                               paste("confidence_max", RNAmod_i, sep = "_"),
                               paste("detected_biotreat", RNAmod_i, sep = "_"),
                               paste("change", RNAmod_i, sep = "_"),
                               paste("anova_pVal", RNAmod_i, sep = "_"),
                               paste("change_directions", RNAmod_i, sep = "_"),
                               paste("change_coherence", RNAmod_i, sep = "_"))
        data.table::data.table(summVal, maxVal, detected, det_times, set_confidence,
                               max_confidence, det_bioT, whichDyn, res_aov, det_direc,
                               coherent_dir) %>%
            set_colnames(IT_RNAmodColNames)
    }) %>% do.call(what = cbind)

    RES_meaDet <- cbind(tmp_txDT, RES_meaDet)

    # Add all metrics all samples (,,,,)
    RES_metrics <- parallel::mclapply(mc.cores = nCores, METRICS, function(metric_i){
        lapply(SETS, function(set_i){
            tmp <- subset(PANORAMA$RES, metric == metric_i & set == set_i)
            tmp <- tmp[match(tmp$pos, tmp_txDT$pos), ]
            round(tmp$score, 4)
        }) %>%
            do.call(what = cbind) %>%
            apply(MARGIN = 1, function(x) paste(x, collapse = ",")) %>%
            stringr::str_replace_all("NA", "")
    }) %>% do.call(what = cbind) %>% set_colnames(METRICS)

    RES_meaDet <- cbind(RES_meaDet, RES_metrics)

    # Detected RNAmods
    tmpD <- lapply(RNAMODS, function(RNAmod_i){
        ifelse(RES_meaDet[[paste("detected", RNAmod_i, sep = "_")]], RNAmod_i, "")
    }) %>%
        do.call(what = cbind) %>%
        apply(MARGIN = 1, function(x) paste(rm_empty(x), collapse = ","))
    RES_meaDet <- tibble::add_column(RES_meaDet, detected_RNAmods = tmpD, .after = "refSeq") %>%
        data.table::data.table()

    # 21nt seq
    RES_meaDet$seq10 <- ""
    tmp <- RES_meaDet[, 1:8]
    # tmp$detected <- tmp$detected_RNAmods != ""
    tmp$detected <- TRUE
    RES_meaDet$seq10[tmp$detected] <- txtools::tx_get_flankSequence(tmp, "detected", upFlank = 10, doFlank = 10) %>%
        suppressWarnings()
    rm(tmp)
    # biotype
    rRNA_biotype <- c("_16S", "_23S", "_5S", "_18s", "_28s", "_5.8s", "_5s") %>%
        paste0("$")
    RES_meaDet$gene_biotype <- ""
    tmpI <- lapply(rRNA_biotype, function(patt_i){
        which(stringr::str_detect(RES_meaDet$gene, patt_i))
    }) %>% unlist() %>% unique()
    RES_meaDet$gene_biotype[tmpI] <- "rRNA"
    RES_meaDet$gene_biotype[stringr::str_detect(RES_meaDet$gene, "tRNA") %>% which()] <- "tRNA"
    list(inspectionTable = RES_meaDet, sets = SETS)
}


#' Resolve RNAmod multicalling
#'
#' Using a percentiles version of the linear scores of the conflicting RNAmods
#' the RNAmod with the highest score according to their percentile will be chosen
#' and written in the column 'detected_RNAmod_resolve'.
#'
#' If still more than one RNAmod that reaches the maximum score then those that
#' reach the maximum score will be output.
#'
#' @param insTable data.table
#' @param models list
#'
#' @return data.table
#' @export
resolve_detectedRNAmods <- function(insTable, models){
    # Detect those positions with multi-calling that require resolution
    detected_ind <- which(insTable$detected_RNAmods != "")
    multiCall_ind <- which(insTable$detected_RNAmods %>% remove_Nm() %>% remove_trailCommas() %>% stringr::str_detect(","))
    multiCall_ind_ok <- setdiff(detected_ind, multiCall_ind)
    multiCall_ind_Nm <- multiCall_ind[lapply(c("Am", "Cm", "Gm", "Um"), function(Nm_i){
        str_detect(insTable$detected_RNAmods[multiCall_ind], Nm_i)
    }) %>% do.call(what = rbind) %>% apply(MARGIN = 2, any)]

    # Add pctScore - Percentile score, based on the median(UnmodScore) = 0, 95pct Modified = 100
    pctScores <- lapply(RNAmods_vec, function(RNAmod_i){
        THRS <- models[[RNAmod_i]]$lmScore_thrs
        selVar <- paste0("lmScore_max_", RNAmod_i)
        varName <- paste0("pctScore_", RNAmod_i)
        pctScore <- ((insTable[[selVar]] - THRS["median_baseNuc"])/(THRS["modified_95pct"] - THRS["median_baseNuc"])) %>%
            magrittr::multiply_by(100) %>% round(2) %>% limit_0_100()
        data.table::data.table(pctScore) %>% magrittr::set_names(varName)
    }) %>% do.call(what = cbind)
    insTable <- cbind(insTable, pctScores)

    # Final Selection Process
    insTable$detected_RNAmod_resolve <- ""
    insTable[multiCall_ind_ok,]$detected_RNAmod_resolve <- insTable[multiCall_ind_ok,]$detected_RNAmods
    if(length(multiCall_ind) > 0){
        for(pos_i in multiCall_ind){
            detRNAmods <- insTable[pos_i,]$detected_RNAmods %>% remove_Nm() %>% remove_trailCommas() %>%
                stringr::str_split(pattern = ",") %>% unlist()
            selVars <- paste0("pctScore_", detRNAmods)
            selScores <- insTable[pos_i, selVars, with = FALSE] %>% as.numeric()
            insTable[pos_i,]$detected_RNAmod_resolve <- paste(detRNAmods[which.max_all(selScores)], collapse = ",")
            if(pos_i %in% multiCall_ind_Nm){
                addNm <- stringr::str_replace_all(insTable$refSeq[pos_i], c(A = "Am", C = "Cm", G = "Gm", 'T' = "Um"))
                insTable[pos_i,]$detected_RNAmod_resolve <- paste(c(insTable[pos_i,]$detected_RNAmod_resolve, addNm), collapse = ",")
            }
        }
    }
    insTable
}

#' Resolve RNAmod multicalling by distance
#'
#' Resolves conflicting RNAmod multicalling by calculating the geometric median
#' of the training sites and the distance of the conflicting sites to them per
#' SET (sample).
#'
#' The space where the geometric median is calculated is in one with all the metrics
#' involved in estimating the conflicting RNAmods.
#'
#' Additionally it generates a consensus of the calledRNAmods per position for
#' all the samples. And a consensus of 2-O-methylations detection confidence.
#'
#' @param insTable
#' @param PANORAMA
#' @param PANORAMA_training
#' @param models
#' @param blacklist_train
#' @param nCores
#'
#' @return
#' @export
#'
resolve_RNAmodsByDist <- function(insTable, PANORAMA, PANORAMA_training, models,
                                  blacklist_train = NULL, nCores = 1){
    SETS <- names(PANORAMA$CALLS)
    # Resolve multicallings by distance to geometric mean (analog to centroid but instead of using the mean value use the median)
    allSetsResolveByDist <- parallel::mclapply(mc.cores = nCores, SETS, function(set_i){
        tmpO <- data.table::data.table(pos = PANORAMA$CALLS[[set_i]][[1]]$pos, resolveNucByDist = "")
        RNAMODS <- check_whichRNAmods_pred(PANORAMA, set_i) %>% base::setdiff(c("Am", "Cm", "Gm", "Um"))
        tmpR <- lapply(RNAMODS, function(RNAmod_i){
            PANORAMA$CALLS[[set_i]][[RNAmod_i]]$pred
        }) %>% do.call(what = cbind) %>% magrittr::set_colnames(RNAMODS) %>%
            tibble::as_tibble()
        tmpR[is.na(tmpR)] <- FALSE
        tmpRow <- which(rowSums(tmpR, na.rm = TRUE) == 1)
        tmpO$resolveNucByDist[tmpRow] <- apply(
            tmpR[tmpRow,],
            MARGIN = 1, function(x){RNAMODS[which(x)]}, simplify = FALSE) %>% unlist()
        tmpRow <- which(rowSums(tmpR, na.rm = TRUE) > 1)
        tmpconfRNAmods <- apply(tmpR[tmpRow,], MARGIN = 1, function(x){
            paste(as.character(na.omit(RNAMODS[x])), collapse = ",")}, simplify = FALSE) %>% unlist()
        confRows <- data.table::data.table(row = tmpRow, set = set_i,
                                           pos = PANORAMA$CALLS[[set_i]][[1]][tmpRow, ]$pos,
                                           refSeq = insTable[tmpRow, ]$refSeq,
                                           confRNAmods = tmpconfRNAmods)
        confGROUPS <- unique(dplyr::select(confRows, c("refSeq", "confRNAmods")))

        # Calculate closest centroid to sites with multiple calling
        centroidRES <- lapply(seq(nrow(confGROUPS)), function(row_i){
            RNAmodgroup_i <- confGROUPS$refSeq[row_i]
            confNucs <- confGROUPS$confRNAmods[row_i]
            confRNAmods <- confNucs %>% stringr::str_split(",") %>% unlist()
            selRNAnucs <- confRNAmods
            selVars <- lapply(selRNAnucs, function(RNAmod_i){
                models[[RNAmod_i]]$vars
            }) %>% unlist() %>% unique()
            NArateDT <- NArate_metricNuc(PANORAMA_training, metrics = selVars, nucs = selRNAnucs) %>%
                dplyr::group_by(metric) %>% dplyr::summarize(maxNArate = max(NArate), minMeasured = min(measuredSites))
            selVars <- dplyr::filter(NArateDT, minMeasured > 1)$metric # Select metrics with more than one measured site per RNAmod
            test_Pos <- dplyr::filter(confRows, confRNAmods == confNucs)$pos
            testDATA <- dplyr::filter(PANORAMA$RES, set == set_i & pos %in% test_Pos & metric %in% selVars) %>%
                tibble::add_column(trainGroup = "test")
            if(!"nuc" %in% names(testDATA)){
                testDATA <- tibble::add_column(testDATA, nuc = "", .after = "refSeq")
            }
            trainDATA <- dplyr::filter(PANORAMA_training$RES, metric %in% selVars & nuc %in% selRNAnucs) %>%
                tibble::add_column(trainGroup = "train")
            if(!is.null(blacklist_train)){
                trainDATA <- dplyr::filter(trainDATA, !pos %in% blacklist_train$pos)
            }
            bothDATA <- rbind(trainDATA, testDATA)
            transDATA <- bothDATA %>% dplyr::group_by(metric) %>% dplyr::mutate(Zscore = as.numeric(scale(score))) %>%
                data.table::data.table() %>% dplyr::select(!"score") %>% tidyr::pivot_wider(names_from = metric, values_from = Zscore) %>%
                data.table::data.table()
            # Calculate centroids of selected nucleotides in training data
            centroids <- lapply(selRNAnucs, function(nuc_i){
                tmpDT <- dplyr::filter(transDATA, nuc %in% nuc_i, trainGroup == "train") %>%
                    dplyr::select(all_of(selVars)) %>% stats::na.omit()
                RES_centroids <- suppressWarnings(pracma::geo_median(as.matrix(tmpDT)))[["p"]] %>%
                    data.frame(row.names = selVars) %>% t() %>% data.table::data.table()
                data.table::data.table(nuc = nuc_i, RES_centroids)
            }) %>% do.call(what = rbind)
            centroids <- dplyr::filter(centroids, nuc %in% confRNAmods)

            # Calculate distance of testdata to centroids
            testPos <- dplyr::filter(transDATA, trainGroup == "test") %>% dplyr::select(tidyselect::all_of("pos"))
            dist_mat <- stats::dist(rbind(dplyr::filter(transDATA, trainGroup == "test") %>% dplyr::select(tidyselect::all_of(selVars)),
                                          centroids[,-1])) %>% as.matrix() %>% data.table::data.table() %>% tibble::tibble()
            dist_mat_nuc <- dist_mat[1:(nrow(dist_mat) - nrow(centroids)),
                                     (ncol(dist_mat) - nrow(centroids)+1):ncol(dist_mat)] %>%
                data.table::data.table() %>% magrittr::set_colnames(paste0("distScore_", confRNAmods))
            # plot(dist_mat_nuc$dist_A, dist_mat_nuc$dist_m66A, col = fct_drop(selRES$nuc), pch = 16)
            closestNuc <- apply(dist_mat_nuc, MARGIN = 1, function(x){
                selRNAnucs[which.min(x)]
            })
            data.table::data.table(testPos, resolveNucByDist = closestNuc)
        }) %>% do.call(what = rbind)
        tmpO$resolveNucByDist[match(centroidRES$pos, tmpO$pos)] <- centroidRES$resolveNucByDist
        tmpO
    })

    RNAMODS <- check_whichRNAmods(PANORAMA) %>% setdiff(c("Am", "Cm", "Gm", "Um"))
    tmpSRES <- lapply(allSetsResolveByDist, function(x) x$resolveNucByDist) %>%
        do.call(what = cbind)
    tmpRES_OUT1 <- parallel::mclapply(seq(nrow(tmpSRES)), mc.cores = nCores, function(i){paste(tmpSRES[i,], collapse = ",")}) %>% unlist()

    timesRNAmodDet <- lapply(RNAMODS, function(RNAmod_i){
        rowSums(tmpSRES == RNAmod_i)
    }) %>% do.call(what = cbind) %>% magrittr::set_colnames(RNAMODS)

    tmpRES_OUT2 <- rep("", nrow(tmpSRES))
    selRows <- which(rowSums(timesRNAmodDet) > 0)
    tmpRES_OUT2[selRows] <- lapply(selRows, function(i){
        RNAMODS[which.max_all(timesRNAmodDet[i,])]
    }) %>% lapply(function(x) paste(x, collapse = ",")) %>% unlist()

    noDetsVec <- paste(rep(",", length(SETS) - 1), collapse = "")
    detectedNmBySets <- rep(noDetsVec, nrow(insTable))
    tmpCols <- paste("confidence", c("Am", "Cm", "Gm", "Um"), sep = "_")
    for(col_i in tmpCols){
        selRows <- which(insTable[[col_i]] != noDetsVec & !is.na(insTable[[col_i]]))
        detectedNmBySets[selRows] <- insTable[[col_i]][selRows]
    }

    tmpNmConf <- stringr::str_split(detectedNmBySets, ",", simplify = TRUE)

    timesNmConf <- lapply(c("l", "h"), function(conf_i){
        rowSums(tmpNmConf == conf_i)
    }) %>% do.call(what = cbind) %>% magrittr::set_colnames(c("l", "h"))

    # Nm confidence consensus
    tmpRES_OUT3 <- rep("", nrow(insTable))
    selRows <- which(rowSums(timesNmConf) > 0)
    tmpRES_OUT3[selRows] <- lapply(selRows, function(i){
        c("l", "h")[which.max_all(timesNmConf[i,])]
    }) %>% lapply(function(x) paste(x, collapse = ",")) %>% unlist()
    tmpRES_OUT3 <- stringr::str_replace_all(tmpRES_OUT3, c('l,h' = "tie", l = "low", h = "high"))

    NmReplace <- c(A = "Am", C = "Cm", G = "Gm", 'T' = "Um")
    selRowsX <- which(tmpRES_OUT2 != "")
    selRowsNm <- which(tmpRES_OUT3 != "")
    selRowsXm <- intersect(selRowsNm, selRowsX)
    tmpCombDet <- data.table::data.table(tmpRES_OUT2[selRowsXm],
                                         stringr::str_replace_all(insTable$refSeq[selRowsXm], NmReplace))
    tmpCombDetOut <- lapply(seq(nrow(tmpCombDet)), function(row_i){
        paste(tmpCombDet[row_i,], collapse = ",")
    }) %>% unlist()

    combined_detectedRNAmod <- rep("", nrow(insTable))
    combined_detectedRNAmod[selRowsX] <- tmpRES_OUT2[selRowsX]
    combined_detectedRNAmod[selRowsNm] <- stringr::str_replace_all(insTable$refSeq[selRowsNm], NmReplace)
    combined_detectedRNAmod[selRowsXm] <- tmpCombDetOut

    OUT_consens <- cbind(select(PANORAMA$CALLS[[1]][[1]], "pos"),
                         detectedRNAmodBySets = tmpRES_OUT1,
                         detectedRNAmodConsensus = tmpRES_OUT2,
                         detectedNmConfidConsens = tmpRES_OUT3,
                         combined_detected_RNAmod = combined_detectedRNAmod)

    dplyr::left_join(insTable, OUT_consens, by = "pos")
}
