# Singularity Container

This folder contains the recipe and helper files to build the singularity image.

## Requirements

- Singularity v3.50

## Software 

The image contains the following software
- R 3.6.3 with Bioconductor Version 3.10 and all packages specified in ListOfPackages.txt
- CellRanger 3.1.0

## Building the image

1. Download cellranger-3.1.0.tar.gz from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and place in this folder
2. Redefine the paths of the files in l8-10 of CovidPBMC.def (absolute paths required)
3. Run 'sudo singularity build CovidPMBC.sif CovidPBMC.def'


## Using the image

Once the image has been built, it can be used by calling:

```
singularity exec CovidPBMC.sif Rscript MyScript.R
```

### Interaction with the host environment

Singularity mounts the host $HOME directory which means that some .files are read, potentially leading to unexpected results.
In particular for this image that means a .Renviron file in the user's $HOME directory can overwrite R environment variables, for instance R_LIBS, which will lead to the singularity instance trying to load packages from the host.
Invoking singularity with '-c' prevents all interactions with the host environment, which requires specifying a $HOME directory with '-H'. Note that all files that singularity needs to access need to be contained within subdirectories of the path specified in -H. 

```
singularity exec -c -H ~/Cambridge/CovidPBMC/ CovidPBMC.sif Rscript Myscript.R
```

### EBI Cluster

As this image is build using singularity v3.50, the respective module will have to be loaded on the EBI cluster.

```
module load singularity/3.5.0
```
