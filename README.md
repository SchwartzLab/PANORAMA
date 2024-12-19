
# PANORAMA <img src='man/figures/logo.png' align="right" height="150" /></a>

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-1.2.5-blue.svg)](https://github.com/SchwartzLab/PANORAMA)
<!-- badges: end -->

## Description

**PANORAMA** is a package that processes Pan-Mod-seq-based data, making
use of the [**txtools**](https://github.com/AngelCampos/txtools)
package.

## Installation

You can install the development version from
[GitHub](https://github.com/SchwartzLab/PANORAMA) typing in the
following commands in the R console:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("SchwartzLab/PANORAMA", build_vignettes = FALSE)
```

If you want to build the vignette, please install with the following
command (only works while connected to WEXAC).

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("SchwartzLab/PANORAMA", build_vignettes = TRUE)
```

## Further documentation

To open the “Using PANORAMA at WEXAC” vignette run the following
command.

``` r
vignette("Using_PANORAMA_at_WEXAC", package = "PANORAMA")
```

## Current limitations

## Licence

PANORAMA is developed under the [Artistic Licence
2.0](https://opensource.org/licenses/Artistic-2.0).
