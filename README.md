
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *aurora* <a href="https://dalimilbujdos.github.io/aurora/"><img src="man/figures/logo.png" align="right" height="138" alt = ""/></a>

<!-- badges: start -->
<!-- badges: end -->

Description: A primary goal of microbial Genome Wide Association Studies
is identifying genomic variants associated with habitats. Existing tools
fail to identify known causal variants if the analyzed trait (like
habitat adaptation) shaped the phylogeny. Due to the presence of
allochthonous strains or metadata errors the stated sources of microbial
strains in public databases are often incorrect. We developed a new tool
aurora (<https://github.com/DalimilBujdos/aurora>) that identifies
autochthonous strains and the genes associated with habitats while
acknowledging the potential role of the trait in shaping the phylogeny.
We successfully validated aurora by applying it to simulated datasets,
and multiple species with various traits.

## Installation

You can install the development version of aurora from
[GitHub](https://github.com/DalimilBujdos/aurora/) with:

``` r
# install.packages("devtools")
devtools::install_github("DalimilBujdos/aurora")
```

For instructions on how to use the package, please refer to
[vignettes](https://dalimilbujdos.github.io/aurora/articles/panGWAS.html)
