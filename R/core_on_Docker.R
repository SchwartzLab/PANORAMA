#' Make STAR genome
#'
#' Wrapper for STAR genome generation
#'
#' @param fastaGenome
#' @param bedAnnotation
#' @param outDir
#' @param nCores
#' @param maxReadLength
#'
#' @return
#' @export
#'
#' @examples
mkSTARgenome_docker <- function(fastaGenome, bedAnnotation = NULL, outDir = NULL,
                         nCores = 2, maxReadLength = 101){
    outName <- paste0(strsplit(fastaGenome, split = "/") %>% unlist %>% utils::tail(1),
                      ".STAR")
    if(is.null(outDir)){
        mkTmpDir()
        outDir <- file.path(getwd(), "PANORAMA_tmpDir", outName)
    }else if(!is.null(outDir)){
        outDir <- outDir
    }
    genome <- txtools::tx_load_genome(fastaGenome)
    lGen <- sum(Biostrings::width(genome))
    nRef <- length(Biostrings::width(genome))
    if(nRef > 1000){
        genomChrBinNbits <- floor(min(18,log2(max(lGen/nRef, maxReadLength))))
    }else{
        genomChrBinNbits <- 18
    }
    genomeindexNb <- floor(min(14, log2(lGen)/2 - 1))
    if(is.null(bedAnnotation)){
        com <- paste("/STAR-2.7.10b/source/STAR",
                     "--runMode genomeGenerate",
                     "--runThreadN", nCores,
                     "--genomeDir", outDir,
                     "--genomeFastaFiles", fastaGenome,
                     "--genomeChrBinNbits", genomChrBinNbits,
                     "--genomeSAindexNbases", genomeindexNb)
        system(com)
    }else{
        tmpGTF <- mkGTF(bedAnnotation = bedAnnotation, outName = tempfile())
        com <- paste("/STAR-2.7.10b/source/STAR",
                     "--runMode genomeGenerate",
                     "--runThreadN", nCores,
                     "--genomeDir", outDir,
                     "--genomeFastaFiles", fastaGenome,
                     "--sjdbOverhang", maxReadLength - 1,
                     "--sjdbGTFfile", tmpGTF,
                     "--genomeChrBinNbits", genomChrBinNbits,
                     "--genomeSAindexNbases", genomeindexNb)
        system(com)
        invisible(file.remove(tmpGTF))
    }
    outDir
}

#' STAR alignment
#'
#' Wrapper to perform alignment of sequencing reads to a reference genome
#' using STAR (Dobin) and sorting and indexing using Samtools.
#'
#' @param read1Files character. Path to R1 FASTQ files
#' @param STARgenomeDir character. Path to STAR genome to map to
#' @param pairedEnd logical. Indicating if to expect a Read2 FASTQ file
#' File will be automatically looked for but both need to include "_R1", and "_R2"
#' in their respective file names.
#' @param zipped logical. TRUE as default for zipped FASTQ files, will be read with zcat.
#' If set to FALSE FASTQ files will be read with cat, as not zipped.
#' @param nCores numeric. Number of cores used for alignment and sorting processes.
#' @param outFilterMultimapNmax numeric. Threshold for which a read will have multiple
#' mappings. Default is 10. For unique alignment change the value to 1.
#' @param outDir character. Path to output directory
#' @param alignIntronMax numeric. maximum intron size, if 0, max intron size will be determined
#' by default as (2ˆwinBinNbits)*winAnchorDistNbins (See STAR manual for more info)
#' @param alignEndsType character. type of read ends alignment (See STAR manual for more info)
#' @param otherSTARparams character. Additional parameters not covered by the arguments
#' of this function can be added here, separated by spaces.
#' @param dry logical. If set to TRUE the alignment is not performed, only output
#' are the paths to expected output files.
#' @param tmpDir character. Path to directory to be used to generate
#' intermediate files.
#' @param logSumm logical. Set to FALSE so final logs are not appended to
#' *mappingSummary.txt*.
#'
#' @return
#' @export
alignSTAR_docker <- function(read1Files, STARgenomeDir, pairedEnd = TRUE, zipped = TRUE,
                      nCores = 4, alignEndsType = "Local", alignIntronMax = 0,
                      outFilterMultimapNmax = 10, outDir, otherSTARparams = "",
                      dry = FALSE, tmpDir = NULL, logSumm = TRUE){
    if(!all(grepl(pattern = "_R1", read1Files))){stop("All read1Files must contain the string '_R1'")}
    rootNames <- lapply(read1Files, function(x){strsplit(x, split = "/") %>% unlist %>%
            utils::tail(1)}) %>% unlist %>% gsub(pattern = "_R1(.)+", replacement = "")
    if(!dry & is.null(tmpDir)){mkTmpDir()}
    if(is.null(tmpDir)){
        tmpDir <- "PANORAMA_tmpDir"
    }
    bamFiles <- file.path(tmpDir, paste0(rootNames, "_Aligned.out.bam"))
    BAM <- c(gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam"),
             gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam.bai"))
    outBAM <- lapply(gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam"),
                     function(x){strsplit(x, split = "/") %>% unlist %>%
                             utils::tail(1)}) %>% unlist %>% file.path(outDir, .)
    if(dry){return(outBAM)}
    if(!dir.exists(outDir)){dir.create(outDir)}
    if(!dir.exists(tmpDir)){dir.create(tmpDir)}
    if(zipped){rFCom <- "zcat"}else if(!zipped){rFCom <- "cat"}
    for(read1F in read1Files){
        if(pairedEnd){
            read2F <- gsub(read1F, pattern = "_R1", replacement = "_R2")
            if(!file.exists(read2F)){stop(read2F, " does not exist.")}
        }else if(!pairedEnd){
            read2F <- ""
        }else{stop("pairedEnd must be logical either TRUE or FALSE")}
        outFPrefix <- strsplit(read1F, split = "/") %>% unlist %>% utils::tail(1) %>%
            gsub(pattern = "R1.fastq.gz", replacement = "") %>%
            gsub(pattern = "R1.fastq", replacement = "")
        if(is.null(tmpDir)){
            outFPrefix <- file.path("PANORAMA_tmpDir", outFPrefix)
        }else{
            outFPrefix <- file.path(tmpDir, outFPrefix)
        }
        com <- paste0("/STAR-2.7.10b/source/STAR",
                      " --runMode alignReads",
                      " --runThreadN ", nCores,
                      " --genomeDir ", STARgenomeDir,
                      " --readFilesCommand ", rFCom,
                      " --readFilesIn ", read1F, " ", read2F,
                      " --outFileNamePrefix ", outFPrefix,
                      " --outSAMtype ", "BAM Unsorted",
                      " --outFilterMultimapNmax ", outFilterMultimapNmax,
                      " --alignEndsType ", alignEndsType,
                      " --alignIntronMax ", alignIntronMax,
                      " ", otherSTARparams)
        system(com)
    }
    # Alignment Summary Report
    if(logSumm){
        logFiles <- file.path(tmpDir, paste0(rootNames, "_Log.final.out"))
        RES <- lapply(logFiles, function(x){
            utils::read.delim(file = x, header = FALSE, stringsAsFactors = FALSE)
        })
        # Merge in one table
        summary <- lapply(seq_along(RES), function(x){
            RES[[x]][,2]
        }) %>% do.call(what = cbind)
        rownames(summary) <- RES[[1]][,1]
        colnames(summary) <- rootNames
        outReport <- file.path(outDir, "mappingSummary.txt")
        if(file.exists(outReport)){ # Add columns to existing summary report
            tmp <- data.table::fread(outReport, header = TRUE) %>%
                tibble::column_to_rownames("V1")
            utils::write.table(x = cbind(tmp, summary), file = outReport,
                               sep = "\t", quote = F, col.names = NA)
        }else{
            utils::write.table(x = summary, file = outReport, sep = "\t", quote = F,
                               col.names = NA)
        }
    }
    # Sort and index with samtools
    for(file in bamFiles){
        system(paste0("/usr/bin/samtools sort -o ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam"),
                      " ", file, " -@", nCores))
        system(paste0("/usr/bin/samtools index ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam")))
        system(paste0("rm ", file))
    }
    # Remove all garbage
    if(logSumm){
        garbageSuffix <- c("_SJ.out.tab", "_Log.progress.out", "_Log.out",
                           "_Log.final.out")
    }else{
        garbageSuffix <- c("_SJ.out.tab", "_Log.progress.out", "_Log.out")
    }
    garbage <- file.path(tmpDir, lapply(rootNames, function(x){
        paste0(x, garbageSuffix)
    }) %>% unlist)
    invisible(file.remove(garbage))
    # Move files to output dir
    for(file in BAM){
        system(paste("mv", file, outDir))
    }
    outBAM
}



#' Library complexity report
#'
#' @param META
#' @param maxExtrapolation
#' @param steps
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
libComplexReport_docker <- function(META, maxExtrapolation = 2.01e6, steps = 1e5, verbose = FALSE){
    if(all(file.exists(META$BAM))){
        for(file in META$BAM){
            com <- paste0("/root/miniconda3/envs/preseq/bin/preseq lc_extrap -P -B ",
                          "-e ", maxExtrapolation, " -s ", steps, " ", file, " -o ",
                          gsub(pattern = ".bam$", replacement = ".lce.txt", x = file),
                          " &")
            system(com)
        }
        Sys.sleep(time = 20)
        lastBam <- gsub(pattern = ".bam$", replacement = ".lce.txt", x = META$BAM)
        while(min(difftime(Sys.time(), file.info(lastBam)$mtime, units = "secs"), na.rm = TRUE) < 60){
            Sys.sleep(time = 2)
        }
    }else{
        stop("Files ", paste(META$BAM[!file.exists(META$BAM)], collapse = " "),
             " do not exist")
    }
    if(verbose){cat("DONE: Library complexity reports.")}
}

# Helpers ######
# Make GTF from BED
mkGTF_docker <- function(bedAnnotation, outName = NULL, source = "user"){
    if(is.null(outName)){
        mkTmpDir()
        tmpDir <- file.path(getwd(), "PANORAMA_tmpDir")
        outName <- file.path(tmpDir, "tmp_geneAnnot.gtf")
    }
    com <- paste0("/root/miniconda3/envs/ucsctools/bin/bedToGenePred ",
                  bedAnnotation, " /dev/stdout | /root/miniconda3/envs/ucsctools/bin/genePredToGtf file /dev/stdin ", outName)
    system(com)
    tmp <- utils::read.delim(outName, header = FALSE)
    tmp[,2] <- source
    utils::write.table(tmp, file = outName, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    outName
}
