Running Seurat on Stanford Sherlock

1. LOGIN TO STANFORD SHERLOCK
ssh <username>@login.sherlock.stanford.edu

2. RETRIEVE DATA
ml python/3.9.0
pip install gdown # for downloading big Google drive files
gdown --id 1sx-kWG-cBc45mItsFBFj5xQwDjwCiaaC -O TotalTissue.combined.harmony.RData # the id was taken from the shareable Google drive link https://drive.google.com/file/d/<id>/view

3. SET UP FOR R SESSION
sh_dev -c 4 # initiate session with 4 CPUs
ml glpk # required to properly install igraph
ml R/4.3.2 # load R
ml harfbuzz freetype fribidi # needed for installing some later packages

4. INSTALL REQUIRED PACKAGES
R # start R session, choose CRAN mirror 69 (USA:OR)
install.packages("igraph", Ncpus=4)
install.packages("leiden", Ncpus=4)
install.packages("Seurat", Ncpus=4) # igraph and leiden are required before installing Seurat
packages <- c("dplyr", "magrittr", "data.table", "Matrix", "devtools", "RcppArmadillo", "Rcpp", "scales", "pheatmap", "gplots", "ggplot2", "cowplot", "tibble", "data.table", "remotes") # install the rest
install.packages(packages, Ncpus=4)
remotes::install_github("immunogenomics/harmony") # need batch correction package Harmony
remotes::install_github("immunogenomics/presto") # need presto to speed up marker finding in Seurat

5. RUN R SCRIPT
# Create R script in sbatch file and upload to desired directory on Sherlock
cd /scratch/users/<username>/<foldername> # change with desired directory
ls # check that sbatch file is there
sbatch R_main.sbatch
squeue -u <username> # check status of run
htop # see compute usage
# Look at slurm.out file for report on run