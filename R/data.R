#' PANORAMA base columns
#'
#' Used internally by PANORAMA but may be useful
#'
#' @format character
#'
"panorama_baseCols"

#' PANORAMA base coordinate columns
#'
#' Used internally by PANORAMA but may be useful
#'
#' @format character
#'
"panorama_baseCoorCols"

#' RNAmods vector
#'
#' RNAmods that Pan-Mod-seq will attempt to detect
#'
#' @format character
#'
"RNAmods_vec"

#' RNAmods vector (extended)
#'
#' As RNAmods_vec, but includes RNAmods that are combinations of RNAmods detected
#' by PANORAMA
#'
#' @format character
#'
"RNAmods_vec_ext"

#' RNAmods vector (Plus)
#'
#' As RNAmods_vec_ext, but includes RNAmods that we have encountered along PANORAMA
#'
#' @format character
#'
"RNAmods_vec_plus"

#' Pseudouridine scores
#'
#' @format character
#' @aliases Y_scores
#'
"Y_metrics"

#' 2-O-methylation scores
#'
#' @format character
#' @aliases Nm_scores
"Nm_metrics"

#' ac4C scores
#'
#' @format character
#' @aliases ac4C_scores
#'
"ac4C_metrics"

#' m1A scores
#'
#' @format character
#' @aliases m1A_scores
#'
"m1A_metrics"

#' m7G scores
#'
#' @format character
#' @aliases m7G_scores
#'
"m7G_metrics"

#' m5C scores
#'
#' @format character
#' @aliases m5C_scores
#'
"m5C_metrics"

#' m3U scores
#'
#' @format character
#' @aliases m3U_scores
#'
"m3U_metrics"

#' m1acp3Y metrics
#'
#' @format character
#'
"m1acp3Y_metrics"

#' m66A metrics
#'
#' @format character
#'
"m66A_metrics"

#' D_metrics
#'
#' @format character
#'
"D_metrics"

#' m1G_metrics
#'
#' @format character
#'
"m1G_metrics"

#' m3C_metrics
#'
#' @format character
#'
"m3C_metrics"

#' m22G_metrics
#'
#' @format character
#'
"m22G_metrics"

#' All default RNAmods' metrics
#'
#' Mostly used for internal PANORAMA working
#'
#' @format list
#'
"metricsList"

#' All is.RNAmod functions
#'
#' Mostly used for internal PANORAMA working
#'
#' @format list
#'
"RNAModFunList"


#' Yeast PANORAMA object
#'
#' Example PANORAMA object derived from Pan-Mod-seq data of yeast.
#'
#' @format list
#'
"yeast_PANORAMA"

#' Yeast rRNA modifications reference
#'
#' rRNA modifications in yeast according to Taoka et al.,
#' 2016 (Nucleic acids research)
#'
#' @format data.frame
#'
"rRNAmods_Sc_Taoka"

#' Variables list for experimental design table creation
#'
#' @format list
#'
"vList"

#' RNAmod metrics table
#'
#' Table with all RNAmod metrics used in PANORAMA with their respective RNAmod,
#' groups, RTase, and function used to generate it.
#'
#' @format data.frame
#'
"RNAmod_metrics"

#' All metrics vector
#'
#' A vector containing all PANORAMA metrics
#'
#' @format character
#'
"allMetrics"

#' RNAmods' base nucleotides
#'
#' A list containing the base nucleotides for different RNAmods used in PANORAMA
#'
#' @format list
#'
"RNAMod_baseNuc"

#' RNAmods' base nucleobases
#'
#' A list containing the base nucleobases for different RNAmods used in PANORAMA
#'
#' @format list
#'
"RNAMod_nucRef"

#' List of confounding RNAmods
#'
#' A list containing the confounding RNAmods per RNAmod identified by PANORAMA
#'
#' @format list
#'
"confRNAmods_list"
