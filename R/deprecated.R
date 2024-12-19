
#' Add scoreA 5prime
#'
#' @param PANORAMA
#' @param group_A
#' @param newColName
#' @param onNucs
#' @param minMedCov
#'
#' @export
panorama_add_scoreA5p <- function(PANORAMA, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    lifecycle::deprecate_warn("0.0.6", "panorama_add_scoreA5p()", "add_scoreA5p()")
    add_scoreA5p(PANORAMA, group_A, newColName, onNucs, minMedCov)
}

#' Add scoreA 3prime
#'
#' @param PANORAMA
#' @param group_A
#' @param newColName
#' @param onNucs
#' @param minMedCov
#'
#' @export
panorama_add_scoreA3p <- function(PANORAMA, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    lifecycle::deprecate_warn("0.0.6", "panorama_add_scoreA3p()", "add_scoreA3p()")
    add_scoreA3p(PANORAMA, group_A, newColName, onNucs, minMedCov)
}

#' Make calls (deprecated)
#'
#' @param PANORAMA list. PANORAMA object
#' @export
panorama_makeCalls <- function(PANORAMA){
    lifecycle::deprecate_warn("0.0.6", "panorama_makeCalls()", "panorama_assignScores()")
    panorama_assignScores(PANORAMA)
}


#' Add all default metrics for Pan-Mod-seq v1
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated as the Pan-Mod-seq library treatments were expanded.
#' This functions applied only for the library treatments of lib 499.
#'
#' @param PANORAMA list. PANORAMA object as output by \code{\link{panorama_PANORAMA}}
#' @export
add_default_metrics_v1 <- function(PANORAMA){
    lifecycle::deprecate_warn("0.0.6", "add_default_metrics_v1()", "add_default_metrics()")
        PANORAMA %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_scoreA3p("Mock.TGIRT") %>%
        add_scoreA3p("Mock.SSIII") %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = TRUE) %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII") %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRD1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT")
}

#' is. 2-O-methylation
#'
#' @param nucs
#'
#' @export
is.2Ome <- function(nucs){
    lifecycle::deprecate_warn("0.1.5", "is.2Ome()", "is.Nm()")
    nucs %in% c("Am", "Gm", "Um", "Cm", "Ym")
} # Any 2-O-methylation


#' Add known RNAmods to PANORAMA object
#'
#' Adds a column 'nuc' to RES and CALLS, using the column 'pos' as reference
#' Reference RNAmod annotation should have at least two columns:
#'  - 'pos': Position given by 'gene:txcoor'
#'  - 'nuc': RNA modification's short name.
#'
#' @param PANORAMA
#' @param RNAmods
#'
#' @export
addKnownRNAmods <- function(PANORAMA, RNAmods){
    lifecycle::deprecate_warn("0.1.5", "addKnownRNAmods()", "add_knownRNAmods()")
    add_knownRNAmods(PANORAMA = PANORAMA, RNAmods = RNAmods)
}
