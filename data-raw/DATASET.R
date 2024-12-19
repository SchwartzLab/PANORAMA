## code to create PANORAMA's in-package objects

# Other objects ################################################################
#' @importFrom rlang .data

RNAMod_baseNuc <- list(Y = "T", Am = "A", Cm = "C", Gm = "G", Um = "T", m5C = "C",
                       ac4C = "C", m1A = "A", m7G = "G", m3U = "T", m66A = "A",
                       m1acp3Y = "T", D = "T", ho5C = "C", m1G = "G", m2A = "A",
                       m2G = "G", m3C = "C", m3Y = "T", m4Cm = "C", m22G = "G",
                       Nm = c("A", "C", "G", "T"), m5U = "T", m6A = "A", m2A = "A")

RNAMod_nucRef <- list(Y = "U", Am = "A", Cm = "C", Gm = "G", Um = "U", m5C = "C",
                      ac4C = "C", m1A = "A", m7G = "G", m3U = "U", m66A = "A",
                      m1acp3Y = "U", D = "U", ho5C = "C", m1G = "G", m2A = "A",
                      m2G = "G", m3C = "C", m3Y = "U", m4Cm = "C", m22G = "G",
                      Nm = c("A", "C", "G", "U"), m5U = "U", m6A = "A", m2A = "A")

panorama_baseCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                    "refSeq", "set", "nuc")

panorama_baseCoorCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                        "refSeq")

RNAmods_vec <- c("m1A", "m66A", "m3C", "m5C", "ac4C", "m1G", "m22G", "m7G", "D",
                 "Y", "m3U", "m1acp3Y", "Am", "Cm", "Gm", "Um")

RNAmods_vec_ext <- c("m1A", "m66A", "m3C", "m5C", "ac4C", "m1G", "m22G", "m7G", "D",
                     "Y", "m3U", "m1acp3Y", "m3Y",
                     paste0(c("m1A", "m66A", "m3C", "m5C", "ac4C", "m1G", "m22G",
                              "m7G", "D", "Y", "m3U", "m1acp3Y", "m3Y"), "m"),
                     "Am", "Cm", "Gm", "Um")

RNAmods_vec_plus <- c(RNAmods_vec_ext, "ho5C", "m2G", "m4Cm", "m6A", "m2A", "m5U")

# META default variables list
vList <- list(libNum = c("499", "524", "553", "697"),
              organism = sort(c("Human_Yeast", "PyroAbyss", "TherAcid", "Yeast",
                                "Human", "Hs_Sc", "ThermoKoda", "ClosTherm", "Mouse",
                                "PlanoHalo", "HsSc", "BSubtilis", "PyroFurio", "EColi",
                                "Tacid", "Hvolc", "PAO", "SL1344", "Sacid"), decreasing = TRUE),
              RTase = c("SSIII", "SSIV", "TGIRT", "RTHIV"),
              libTreat = c("Ac4C", "CMC", "DeacetylatedAc4C", "MocklowdNTPs",
                           "Mock", "Dimroth", "m5C", "AlkBmix", "RBSseqHeatMg",
                           "NaBH4HydBiotin", "BIDseq"),
              bioTreat = sort(c("10C", "15C", "25C", "30C", "37C", "42C", "50C",
                                "60C", "65deg", "75deg",
                                "78deg", "80deg", "85deg", "95deg", "100deg",
                                "102deg", "AcidpH1", "AcidpH2", "AcidpH3",
                                "pH1", "pH2", "pH3", "pH1p4_42C", "pH2p4_30C",
                                "pH2p4_42C", "pH2p4_50C", "pH4_42C", "OD1_25C",
                                "OD1_30C", "OD1_37C", "OD2_25C", "OD2_30C",
                                "OD2_37C", "pH2p5_75C", "pH3_55C", "pH3p5_75C",
                                "pH5_75C", "pHunk_85C", "WT_10mth", "WT_30mth",
                                "CR_30mth", "RM_30mth", "3T3"), decreasing = TRUE),
              replicate = c(paste0("rep", 1:6), "S5", "S6"))



usethis::use_data(vList, overwrite = TRUE)
usethis::use_data(RNAMod_baseNuc, overwrite = TRUE)
usethis::use_data(RNAMod_nucRef, overwrite = TRUE)
usethis::use_data(panorama_baseCols, overwrite = TRUE)
usethis::use_data(panorama_baseCoorCols, overwrite = TRUE)
usethis::use_data(RNAmods_vec, overwrite = TRUE)
usethis::use_data(RNAmods_vec_ext, overwrite = TRUE)
usethis::use_data(RNAmods_vec_plus, overwrite = TRUE)

# RNAmod METRICS ###############################################################
# Metrics to be assigned to specific RNAmods
Y_metrics <- c("SRD1bpDS_CMC.TGIRT_Mock.TGIRT",
               "SRD1bpDS_CMC.SSIII_Mock.SSIII",
               "SRlog2FCh1bpDS_CMC.TGIRT_Mock.TGIRT",
               "SRlog2FCh1bpDS_CMC.SSIII_Mock.SSIII")
# "DRD_BIDseq.SSIII_Mock.SSIII" - not used

Nm_metrics <- c("NmStopScore_MocklowdNTPs.SSIII_Mock.SSIII",
                "ScoreA3p_Mock.TGIRT",
                "ScoreA3p_Mock.SSIII",
                "ScoreA3p_Mock.RTHIV",
                "ScoreC3p_Mock.RTHIV",
                "ScoreC3p_Mock.SSIII",
                "ScoreC3p_Mock.TGIRT")

ac4C_metrics <- c("CtoT.MRD_Ac4C.TGIRT_Mock.TGIRT",
                  "CtoT.MRD_Ac4C.TGIRT_DeacetylatedAc4C.TGIRT",
                  "CtoT.MRD_Ac4C.SSIII_Mock.SSIII",
                  "CtoT.MRD_Ac4C.SSIII_DeacetylatedAc4C.SSIII")

m1A_metrics <- c("SRD1bpDS_Mock.SSIII_Dimroth.SSIII",
                 "SRlog2FCh1bpDS_Mock.SSIII_Dimroth.SSIII",
                 "MRD_Mock.TGIRT_Dimroth.TGIRT",
                 "MRD_Mock.SSIII_RBSseqHeatMg.SSIII",
                 "MRD_Mock.TGIRT_RBSseqHeatMg.TGIRT",
                 "MRD_Mock.TGIRT_AlkBmix.TGIRT")

m7G_metrics <- c("MRD_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                 "SRD1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                 "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                 "SRlog2FCh1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                 "SRlog2FCh1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                 "MRD_NaBH4HydBiotin.RTHIV_Mock.RTHIV",
                 "MRD_NaBH4HydBiotin.TGIRT_Mock.TGIRT")

m5C_metrics <- c("CytPer_m5C.TGIRT_Mock.TGIRT",
                 "CytPer_RBSseqHeatMg.SSIII_Mock.SSIII",
                 "CytPer_RBSseqHeatMg.TGIRT_Mock.TGIRT")

m3U_metrics <- c("MRD_Mock.TGIRT_AlkBmix.TGIRT",
                 "MRD_Mock.TGIRT_Mock.SSIII",
                 "MRD_Mock.RTHIV_Mock.TGIRT")

m1acp3Y_metrics <- c("SRD1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT", # Metric suggested by analysis
                     "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII", # Metric suggested by analysis
                     "SRD1bpDS_Mock.SSIII_Dimroth.SSIII", # Metric suggested by analysis
                     "MRD_Mock.SSIII_RBSseqHeatMg.SSIII", # Metric suggested by analysis
                     "MRD_Mock.TGIRT_RBSseqHeatMg.TGIRT",
                     "SR1bpDS_Mock.SSIII",
                     "MR_Mock.RTHIV")

m66A_metrics <- c("MRD_Mock.SSIII_RBSseqHeatMg.SSIII", # Metric suggested by analysis
                  "MRD_Mock.TGIRT_RBSseqHeatMg.TGIRT", # Metric suggested by analysis
                  "DRD_RBSseqHeatMg.TGIRT_Mock.TGIRT", # Metric suggested by analysis
                  "MRDAtoG_Mock.TGIRT_m5C.TGIRT", # Metric suggested by analysis
                  "MRDAtoG_Mock.SSIII_m5C.SSIII") # Metric suggested by analysis

# New RNAmods from tRNAs
D_metrics <- c("MRD_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
               "MRD_NaBH4HydBiotin.RTHIV_Mock.RTHIV",
               "MRD_Mock.TGIRT_Dimroth.TGIRT")

m1G_metrics <- c("MRD_Mock.TGIRT_AlkBmix.TGIRT")

m3C_metrics <- c("MRD_Mock.TGIRT_AlkBmix.TGIRT",
                 "MR_Mock.TGIRT")

m22G_metrics <- c("MR_Mock.TGIRT",
                  # "MR_Mock.RTHIV", # Removed due to low query rate in PANORAMA-3 (only one site measured) and collinearity with the other MR
                  "MR_Mock.SSIII")

# All metrics list
metricsList <- list(m1A = m1A_metrics,
                    m66A = m66A_metrics,
                    m3C = m3C_metrics,
                    m5C = m5C_metrics,
                    ac4C = ac4C_metrics,
                    m1G = m1G_metrics,
                    m22G = m22G_metrics,
                    m7G = m7G_metrics,
                    D = D_metrics,
                    Y = Y_metrics,
                    m3U = m3U_metrics,
                    m1acp3Y = m1acp3Y_metrics,
                    Am = Nm_metrics,
                    Cm = Nm_metrics,
                    Gm = Nm_metrics,
                    Um = Nm_metrics)

# All RNAmod is.X functions
RNAModFunList <- list(Y = is.pseudoU, Am = is.Am, Cm = is.Cm, Gm = is.Gm,
                      Um = is.Um, m5C = is.m5C, ac4C = is.ac4C, m1A = is.m1A,
                      m7G = is.m7G, m3U = is.m3U, m66A = is.m66A,
                      m1acp3Y = is.m1acp3Y, D = is.D, m1G = is.m1G, m3C = is.m3C,
                      m22G = is.m22G, m2G = is.m2G, m3Y = is.m3Y, m4Cm = is.m4Cm,
                      Nm = is.Nm, m6A = is.m6A, m5U = is.m5U, m2A = is.m2A)

# Metrics table
RNAmod_metrics <- lapply(names(metricsList), function(RNAmod){
    cbind(RNAmod, metric = metricsList[[RNAmod]])
}) %>% do.call(what = "rbind") %>% data.table::data.table()
RNAmod_metrics$baseNuc <- unlist(RNAMod_baseNuc[RNAmod_metrics$RNAmod])
RNAmod_metrics$metricFunction <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[1]))


RNAmod_metrics$RTase_A <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[2])) %>%
    stringr::str_extract(pattern = "(TGIRT|SSIII|RTHIV)")
RNAmod_metrics$RTase_B <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[3])) %>%
    stringr::str_extract(pattern = "(TGIRT|SSIII|RTHIV)")
RNAmod_metrics$libTreat_A <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[2])) %>%
    stringr::str_extract(pattern = paste0("(", paste(vList$libTreat, collapse = "|"), ")"))
RNAmod_metrics$libTreat_B <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[3])) %>%
    stringr::str_extract(pattern = paste0("(", paste(vList$libTreat, collapse = "|"), ")"))
RNAmod_metrics$group_A <- paste(RNAmod_metrics$libTreat_A, RNAmod_metrics$RTase_A, sep = ".")
RNAmod_metrics$group_B <- paste(RNAmod_metrics$libTreat_B, RNAmod_metrics$RTase_B, sep = ".")
RNAmod_metrics$group_B[is.na(RNAmod_metrics$libTreat_B)] <- NA
RNAmod_metrics <- RNAmod_metrics[,c(1, 3, 2, 5, 6, 7, 8, 9, 10, 4)]
RNAmod_metrics$RNAmod <- factor(RNAmod_metrics$RNAmod, levels = RNAmods_vec)
RNAmod_metrics <- RNAmod_metrics[order(RNAmod),]

# All metrics vector
allMetrics <- unique(unlist(metricsList))

usethis::use_data(Y_metrics, overwrite = TRUE)
usethis::use_data(Nm_metrics, overwrite = TRUE)
usethis::use_data(ac4C_metrics, overwrite = TRUE)
usethis::use_data(m1A_metrics, overwrite = TRUE)
usethis::use_data(m7G_metrics, overwrite = TRUE)
usethis::use_data(m5C_metrics, overwrite = TRUE)
usethis::use_data(m3U_metrics, overwrite = TRUE)
usethis::use_data(m1acp3Y_metrics, overwrite = TRUE)
usethis::use_data(m66A_metrics, overwrite = TRUE)
usethis::use_data(D_metrics, overwrite = TRUE)
usethis::use_data(m1G_metrics, overwrite = TRUE)
usethis::use_data(m3C_metrics, overwrite = TRUE)
usethis::use_data(m22G_metrics, overwrite = TRUE)
usethis::use_data(metricsList, overwrite = TRUE)
usethis::use_data(RNAModFunList, overwrite = TRUE)
usethis::use_data(RNAmod_metrics, overwrite = TRUE)
usethis::use_data(allMetrics, overwrite = TRUE)

# Taoka Sc sites ###############################################################
rRNAmods_Sc_Taoka <- readRDS("/home/labs/schwartzlab/miguelg/BIGDATA/RNAmod_Annot/Taoka/rib_mods_Sc.rds")
usethis::use_data(rRNAmods_Sc_Taoka, overwrite = TRUE)

# RNAMod functions list
confRNAmods_list <- list(m1A = c("m66A", "Am"),
                         m66A = c("m1A", "Am"),
                         m3C = c("m5C", "ac4C", "Cm"),
                         m5C = c("m3C", "ac4C", "Cm"),
                         ac4C = c("m3C", "m5C", "Cm"),
                         m1G = c("m22G", "m7G", "Gm"),
                         m22G = c("m1G", "m7G", "Gm"),
                         m7G = c("m22G", "m1G", "Gm"),
                         D = c("Y", "m3U", "m1acp3Y", "Um"),
                         Y = c("D", "m1acp3Y", "m3U", "Um"),
                         m3U = c("D", "m1acp3Y", "Y", "Um"),
                         m1acp3Y = c("D", "Y", "m3U", "Um"),
                         Am = c("m66A", "m1A"),
                         Cm = c("m3C", "m5C", "ac4C"),
                         Gm = c("m1G", "m22G", "m7G"),
                         Um = c("D", "Y", "m3U", "m1acp3Y"))

usethis::use_data(confRNAmods_list, overwrite = TRUE)

################################################################################
################################################################################
################################################################################
# Make sample PANORAMA object. ####################################################

# Genome and Gene Annotation
fastaGenome <- "/home/labs/schwartzlab/miguelg/BIGDATA/genome_references/ribosomal/Sc_ribosomal_seqs.fa"
bedAnnotation <- "/home/labs/schwartzlab/miguelg/BIGDATA/gene_annotations/ribosomal/Sc_ribosomal_seqs.bed"
r1Files <- grep(list.files("/home/labs/schwartzlab/miguelg/WORKSPACE_wexac/panorama_seq/data/testDATA",
                           full.names = TRUE), pattern = "R1", value = TRUE)
OUTDIR <- "/home/labs/schwartzlab/miguelg/WORKSPACE_wexac/panorama_seq/data/test_OUT"
EXP_NAME <- "S.cerevisiae_SK1-lib524"
NCORES <- 10

# Experimental design table
META <- panorama_META(fileNames = r1Files,
                   varsList = vList,
                   outDir = OUTDIR)

# Make transcriptome, make new gene annotation for transcriptome
fastaTxOme <- mkTranscriptome(fastaGenome, bedAnnotation, nCores = NCORES, outFile = "data/fastaTxOme")
bedTxOme <- mkBedFromFastaTxOme(fastaTxOme, outFile = "data/bedTxOme")

# Creating bisulphite transcriptome
bisTxPath <- bisGenome(fastaTxOme, outFile = "data/bisTxPath")

# Create STAR genomes
STARGenome <- mkSTARgenome(fastaTxOme, bedTxOme, outDir = "data/STARGenome")
STARGenome_bis <- mkSTARgenome(bisTxPath, bedTxOme, outDir = "data/STARGenome_bis")

#Loading transcriptome and annotation
GENOME <- txtools::tx_load_genome(fastaTxOme)
TXOME <- txtools::tx_load_bed(bedTxOme)

if(!all(file.exists(META$BAM))){
    # STAR Alignment
    alignSTAR(read1Files = META[libTreat != "m5C" & libTreat != "RBSseqHeatMg", FASTQ],
              nCores = NCORES,
              zipped = TRUE,
              STARgenomeDir = STARGenome,
              alignEndsType = "Local",
              outDir = OUTDIR)
    alignSTAR(read1Files = META[libTreat == "m5C" | libTreat == "RBSseqHeatMg", FASTQ],
              nCores = NCORES,
              zipped = TRUE,
              STARgenomeDir = STARGenome_bis,
              alignEndsType = "Local",
              outDir = OUTDIR)
}

if(!all(file.exists(META$RDS))){
    rdsFiles <- parallel::mclapply(mc.cores = NCORES, seq_along(META$FASTQ), function(i){
        bam2TxDT(BAMfile = META$BAM[i],
                 geneAnnot = TXOME,
                 genome = GENOME,
                 dtType = "covNuc",
                 outDir = OUTDIR,
                 nCores = 1,
                 remL = 1000,
                 minR = 0)
    })
}
yeast_PANORAMA <- panorama_PANORAMA(META, GENOME, TXOME, nCores = 1)
# Save object to .rda file
usethis::use_data(yeast_PANORAMA, overwrite = TRUE)
# Remove tmp files and dirs
tmpDirs <- c(STARGenome, STARGenome_bis)
tmpFiles <- c(fastaTxOme, bisTxPath, bedTxOme)
rmdFiles <- file.remove(tmpFiles)
unlink(c(STARGenome, STARGenome_bis), recursive = TRUE)
