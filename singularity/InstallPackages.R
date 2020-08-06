pkgs <- read.table(file="/ListOfPackages.txt",stringsAsFactors=FALSE)
pkgs <- pkgs$V1
.libPaths(c("/usr/local/lib/R/library/",.libPaths()))
# For some reason XML doesn't install, trying to add it manual with a specific repo
install.packages("XML", repos = "http://www.omegahat.net/R")
install.packages(c("BiocManager","stringi"),repos="https://cloud.r-project.org/")
BiocManager::install(update=TRUE, ask=FALSE, version="3.10")
BiocManager::install(pkgs)
