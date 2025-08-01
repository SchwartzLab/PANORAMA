#' Create PANORAMA META - Experiment design matrix
#'
#' PANORAMA META constructor from RDS file names
#'
#' @param fileNames character. Full path to read1-FASTQ files, file names must
#' contain the "_R1" suffix.
#' @param varsList list. A list of char vectors that define a sample and their
#' belonging to a group, and set of samples. Each element of the list must have
#' a name that specifies it.
#' @param groupVars character. The names of variables used to define a group
#' of samples. By default the library treatment ("libTreat") and reverse
#' transcriptase ("RTase") in that order. Changing this only allows for flexibility
#' on the definition of comparable samples. **Should not be changed unless having**
#' **a deep understanding of the PANORAMA package and require to develop new**
#' ** functionality**
#' @param setVars character. Vector of variables from varsList, which combination
#' define the belonging of a group of FASTQ read files to the same biological
#' sample
#' @param idVars character. Vector of variable names from varsList, which
#' combination define the identity of each FASTQ read file. None combination
#' of variables should result in a duplicated id.
#' @param outDir character. Path to directory where all processed files
#' (BAM and RDS) will be written to.
#' @param bisTreats character. Vector with names of library treatments that
#' require to be aligned with a bisulfite version
#'
#' @return
#' @export
#'
panorama_META <- function(fileNames, varsList = PANORAMA::vList,
                       groupVars = c("libTreat", "RTase"),
                       setVars = c("libNum", "organism", "bioTreat", "replicate"),
                       idVars = c("libNum", "organism", "RTase", "libTreat", "bioTreat", "replicate"),
                       outDir,
                       bisTreats = c("m5C", "RBSseqHeatMg")){
    if(magrittr::not("libTreat" %in% names(varsList))){
        stop("libTreat must be one of the element names in varsList argument")
    }
    DT <- data.table::data.table(FASTQ = fileNames)
    procFileNames <- DT$FASTQ %>% stringr::str_split(pattern = "/") %>%
        lapply(function(x) utils::tail(x, 1)) %>% unlist
    DT <- lapply(names(varsList), function(i){
        patt <- paste0("(", paste(varsList[[i]], collapse = "|"), ")")
        tmp <- stringr::str_extract_all(string = procFileNames, pattern = patt)
        lapply(tmp, function(x) paste(x, collapse = ".")) %>% unlist %>%
            data.table::data.table() %>% magrittr::set_names(i)
    }) %>% do.call(what = cbind) %>% cbind(DT, .)
    DT <- lapply(DT, function(x){
        if(!all(x == "")){return(x)}
    }) %>% do.call(what = "cbind") %>% data.table::data.table()
    setIn <- setVars %in% colnames(DT)
    idIn <- idVars %in% colnames(DT)
    groupIn <- groupVars %in% colnames(DT)
    setVars <- setVars[setIn]
    idVars <- idVars[idIn]
    groupVars <- groupVars[groupIn]
    DT$group <-  data.frame(DT[, groupVars, with = FALSE]) %>%
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% factor()
    DT$set <- data.frame(DT[, setVars, with = FALSE]) %>%
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% factor()
    DT$id <- data.frame(DT[, idVars, with = FALSE]) %>%
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% as.character()
    if(sum(duplicated(DT$id)) != 0){
        warning("Generated ids are not unique per sample")
    }
    if(sum(duplicated(paste(DT$group, DT$set)))){
        warning("Combinations of group and set should be unique between samples")
    }
    if(sum(DT$libTreat == "") != 0){
        warning("Some samples lack libTreat")
    }
    rootNames <- DT$FASTQ %>% lapply(function(x){strsplit(x, split = "/") %>%
            unlist %>% utils::tail(1)}) %>% unlist %>% gsub(pattern = "_R1(.)+", replacement = "")
    BAM <- file.path(outDir, paste0(rootNames, "_Aligned.out.sorted.bam"))
    RDS <- file.path(outDir, paste0(rootNames, "_Aligned.out.sorted.txDT.rds"))
    bis <- ifelse(DT$libTreat %in% bisTreats, yes = TRUE, no = FALSE)
    DT <- tibble::add_column(.data = DT, .after = "FASTQ", BAM = BAM) %>%
        tibble::add_column(.after = "BAM", RDS = RDS) %>%
        tibble::add_column(bisAlign = bis)
    return(DT)
}

#' BAM to txDataTable
#'
#' From BAM file to count data in a table
#'
#' @param BAMfile
#' @param geneAnnot
#' @param genome
#' @param dtType
#' @param paired
#' @param outDir
#' @param remL
#' @param minR
#' @param nCores
#' @param ySize
#' @param verb
#'
#' @return
#' @export
#'
#' @examples
bam2TxDT <- function(BAMfile, geneAnnot, genome, dtType, paired = TRUE,
                     outDir, remL = 1000, minR = 50, nCores = 2,
                     ySize = 100000, verb = TRUE){
    t0 <- Sys.time() # Start time
    # Process arguments
    # Load sequence or not
    if(dtType == "cov"){
        lSeq <- FALSE
    }else if(dtType == "covNuc"){
        lSeq <- TRUE
    }else{
        stop("dtType argument must be one of the options: 'cov' or 'covNuc'")
    }
    # Check yieldSize to be integer
    txtools:::check_integer_arg(ySize, "ySize")
    # Set OUTPUT filename
    outName <- strsplit(BAMfile, split = "/") %>% unlist %>% utils::tail(1) %>%
        gsub(x = ., pattern = ".bam$", replacement = ".txDT.rds", perl = T)
    # Output all genes even with no reads overlapping
    if(minR == 0){makeFULL <- TRUE}else{makeFULL <- FALSE}
    if(minR == 0){minR <- 1}
    # MAIN program
    # Load gene annotation
    gA <- geneAnnot
    if(verb){
        cat("Gene annotation loaded with", length(gA), "gene models.\n")
    }
    GENOME <- genome
    # Load BAM file
    bam <- txtools::tx_load_bam(file = BAMfile,
                                yieldSize = ySize,
                                scanFlag = "default",
                                loadSeq = lSeq,
                                verbose = verb,
                                pairedEnd = paired)
    t1 <- Sys.time() # Loading BAM file time
    # Convert to transcriptomic
    txReads <- txtools::tx_reads(reads = bam,
                                 geneAnnot = gA,
                                 nCores = nCores,
                                 minReads = minR,
                                 withSeq = lSeq,
                                 verbose = verb)
    # Filter by length
    if(!is.na(remL)){
        txReads <- txtools::tx_filter_maxWidth(x = txReads, thr = remL, nCores = nCores)
    }

    # Data.table generation
    if(dtType == "cov"){
        if(verb){
            cat("Generating coverage data.table")
        }
        OUT <- txtools::tx_makeDT_coverage(x = txReads,
                                           geneAnnot = gA,
                                           nCores = nCores,
                                           fullDT = makeFULL,
                                           genome = GENOME)
    }else if(dtType == "covNuc"){
        if(verb){
            cat("Generating coverage and nucleotide frequency data.table. \n")
        }
        OUT <- txtools::tx_makeDT_covNucFreq(x = txReads,
                                             geneAnnot = gA,
                                             nCores = nCores,
                                             fullDT = makeFULL,
                                             genome = GENOME)
    }else{print("This message should not be printed, let the maintainer know.")}
    t2 <- Sys.time() # Creating data.table time
    # NOTE: In practice using too many cores made for longer processing times
    # newNCores <- min(nCores, 10)

    # Saving file as .rds
    saveRDS(object = OUT, file = file.path(outDir, outName))
    # Report
    timeBam <- t1 - t0 # Total time taken
    timePrc <- t2 - t1 # Total time taken
    timeTot <- t2 - t0 # Total time taken
    reportName <- strsplit(BAMfile, split = "/") %>% unlist %>% utils::tail(1) %>%
        gsub(x = ., pattern = ".bam$", replacement = ".txDT.log", perl = T)
    readsInOut <- parallel::mclapply(mc.cores = nCores, txReads, names) %>% unlist
    uniqReadsInOut <- unique(readsInOut)
    report <- c("BAM file name:", BAMfile,
                "Paired-end reads in BAM file:", length(bam),
                "Output contains:", " ",
                "Number of genes:", length(unique(OUT$gene)),
                "Number of reads in output:", length(readsInOut),
                "Number of unique reads in output:", length(uniqReadsInOut),
                "Fraction of total reads in output:", round(length(uniqReadsInOut)/length(bam), 4),
                "Loading BAM time:", paste(round(timeBam, 2), units(timeBam), sep = " "),
                "Processing time:", paste(round(timePrc, 2), units(timePrc), sep = " "),
                "Total time taken:", paste(round(timeTot, 2), units(timeTot), sep = " ")) %>%
        matrix(ncol = 2, byrow =T)
    utils::write.table(x = report,
                       file = file.path(outDir, reportName),
                       sep = "\t",
                       quote = FALSE,
                       row.names = FALSE,
                       col.names = FALSE)
    file.path(outDir, outName)
}

#' Add known RNA modifications nucleotide identity
#'
#' Add the nucleotide identity based on a table with columns "pos" and "nuc
#' establishing the position and nucleotide identity
#' Useful to annotate already known RNA modifications. The function will add
#' a column called "nuc" in the "RES" and "CALLS" tables.
#'
#' @param PANORAMA PANORAMA object
#' @param RNAmods data.frame with two columns "pos" for the positions in the
#' same format as the tables inside the PANORAMA object as PANORAMA$RES, a combination
#' in the form of "gene:txcoor"; And "nuc" for the identity of the nucleotide
#' on said position.
#'
#' @return
#' @export
add_knownRNAmods <- function(PANORAMA, RNAmods){
    if("RES" %in% names(PANORAMA)){
        if("refSeq" %in% names(PANORAMA$RES)){
            tmp <- RNAmods$nuc[match(PANORAMA$RES$pos, RNAmods$pos)]
            tmp[which(is.na(tmp))] <- PANORAMA$RES$refSeq[which(is.na(tmp))]
            tmp[tmp == "T"] <- "U"
            PANORAMA$RES <- hlp_removeColumnIfPresent(PANORAMA$RES, "nuc")
            PANORAMA$RES <- tibble::add_column(nuc = factor(tmp),
                .after = "refSeq",
                .data = PANORAMA$RES)
        }else{
            stop("PANORAMA$RES, requires to have the refSeq column.")
        }
    }else{
        stop("RES is not an element of the PANORAMA object. Metrics have not been",
             " calculated yet.")
    }
    if("CALLS" %in% names(PANORAMA)){
        for(set in unique(PANORAMA$META$set)){
            for(i in seq_along(PANORAMA$CALLS[[set]])){
                if("refSeq" %in% names(PANORAMA$CALLS[[set]][[i]])){
                    tmp <- RNAmods$nuc[match(PANORAMA$CALLS[[set]][[i]]$pos, RNAmods$pos)]
                    tmp[which(is.na(tmp))] <- PANORAMA$CALLS[[set]][[i]]$refSeq[which(is.na(tmp))]
                    tmp[tmp == "T"] <- "U"
                    PANORAMA$CALLS[[set]][[i]] <- hlp_removeColumnIfPresent(PANORAMA$CALLS[[set]][[i]], "nuc")
                    PANORAMA$CALLS[[set]][[i]] <- tibble::add_column(
                        nuc = tmp,
                        .after = "refSeq",
                        .data = PANORAMA$CALLS[[set]][[i]])
                }else{
                    stop("PANORAMA$CALLS[[set]][[RNAmod]], requires to have the refSeq column.")
                }
            }
        }
    }
    PANORAMA
}

#' Remove known RNAmods ('nuc' column)
#'
#' @param PANORAMA PANORAMA object
#'
#' @return PANORAMA object
#' @export
#'
rm_knownRNAmods <- function(PANORAMA){
    if("RES" %in% names(PANORAMA)){
        if("nuc" %in% names(PANORAMA$RES)){
            PANORAMA$RES$nuc <- NULL
        }
    }else{
        stop("PANORAMA object should at least contain the 'RES' element so that removing known RNAmods makes sense")
    }
    if("CALLS" %in% names(PANORAMA)){
        SETS <- names(PANORAMA$CALLS)
        for(set_i in SETS){
            for(RNAmod_i in names(PANORAMA$CALLS[[set_i]])){
                if("nuc" %in% names(PANORAMA$CALLS[[set_i]][[RNAmod_i]])){
                    PANORAMA$CALLS[[set_i]][[RNAmod_i]]$nuc <- NULL
                }
            }
        }
    }
    if("EVAL" %in% names(PANORAMA)){
        PANORAMA$EVAL <- NULL
    }
    PANORAMA
}

#' BAM 2 txDT from META
#'
#' @param META
#' @param nCores
#' @param geneAnnot
#' @param genome
#' @param outDir
#' @param force
#' @param verb
#' @param remL
#' @param minR
#'
#' @return
#' @export
#'
#' @examples
bam2txDT_META <- function(META, nCores, geneAnnot, genome, outDir, remL = 1000,
                          minR = 50, force = FALSE, verb = TRUE){
    if(!all(file.exists(META$BAM))){
        stop("Not all BAM files in META exist, make sure alignment ",
             "has been done.")
    }else if(!all(file.exists(META$RDS)) | force){
        misBAM <- META$BAM[which(!file.exists(META$RDS))]
        if(force){
            misBAM <- META$BAM
        }
        lapply(seq_along(misBAM), function(i){
            bam2TxDT(BAMfile = misBAM[i],
                     geneAnnot = geneAnnot,
                     genome = genome,
                     dtType = "covNuc",
                     outDir = outDir,
                     nCores = nCores,
                     remL = remL,
                     minR = minR,
                     verb = verb)
            # Fix genom COors
            if(verb){cat(i, "out of", length(misBAM),
                         "BAM files processed...\n\n")}
        }) %>% invisible()

    }else{
        cat("All RDS files already exist\n")
    }
    if(verb){cat("All bam files procesed!\n")}
}
