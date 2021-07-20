# CNVRG

A package to provide a user-friendly wrapper for 'RStan' to implement Dirichlet-multinomial modeling of count data.

This package is managed by Joshua Harrison.

A github repo for this package exists at: https://github.com/JHarrisonEcoEvo/CNVRG
Please check that repo out for new versions of this software that are under development. Also, if you discover issues or have questions please post them there so that others may benefit from answers.

To install from source. Download repo and unzip it then use the following commands from within the 'R' interpreter. First, you will need to change the path to the 'CNVRG' directory as appropriate. 

install.packages("./CNVRG/", repos = NULL, type = "source")

Alternatively, one can download and install from GitHub using devtools like so: 
devtools::install_github("https://github.com/JHarrisonEcoEvo/CNVRG")

To install from CRAN, open R and type:
install.packages("CNVRG")

You may see a prompt asking if you want to compile the package, type "Yes" if you see this. 

VIGNETTE

The vignette can be [HERE](https://rpubs.com/harrisonjg/792276) and the Rmd can be found at the GitHub repo linked above. Additionally, a video providing the intuition behind the CNVRG model can be found here: https://use.vg/OSVhFJ

For those interested in the details of the model see our paper in Molecular Ecology Resources called, ["Dirichlet‐multinomial modelling outperforms alternatives for analysis of microbiome and other ecological count data"](https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13128)


<!-- badges: start -->
[![Travis build status](https://travis-ci.com/JHarrisonEcoEvo/CNVRG.svg?branch=master)](https://travis-ci.com/JHarrisonEcoEvo/CNVRG)
<!-- badges: end -->
