pkgs <- read.table(file="/ListOfPackages.txt",stringsAsFactors=FALSE)
pkgs <- pkgs$V1
.libPaths(c("/usr/local/lib/R/library/",.libPaths()))
install.packages(c("BiocManager","stringi"),repos="https://cloud.r-project.org/")
BiocManager::install(update=TRUE, ask=FALSE)
BiocManager::install(pkgs)
