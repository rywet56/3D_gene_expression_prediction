#! /usr/lib/R/bin/Rscript

# create directory to install R libraries into
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE) 
.libPaths(Sys.getenv("R_LIBS_USER"))

# enable the creation of a R jupyter lab notebook
install.packages("IRkernel")
IRkernel::installspec()

# install R packages
# install.packages("Seurat")
