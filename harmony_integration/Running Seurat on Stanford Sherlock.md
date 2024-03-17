Running Seurat on Stanford Sherlock <br>

1. LOGIN TO STANFORD SHERLOCK <br>
ssh <username>@login.sherlock.stanford.edu <br>

2. RETRIEVE DATA <br>
ml python/3.9.0 <br>
pip install gdown # for downloading big Google drive files <br>
gdown --id 1sx-kWG-cBc45mItsFBFj5xQwDjwCiaaC -O TotalTissue.combined.harmony.RData # the id was taken from the shareable Google drive link https://drive.google.com/file/d/<id>/view <br>

3. SET UP FOR R SESSION <br>
sh_dev -c 4 # initiate session with 4 CPUs <br>
ml glpk # required to properly install igraph <br>
ml R/4.3.2 # load R <br>
ml harfbuzz freetype fribidi # needed for installing some later packages <br>

4. INSTALL REQUIRED PACKAGES FOR SEURAT ENVIRONMENT <br>
R # start R session, choose CRAN mirror 69 (USA:OR) <br>
install.packages("igraph", Ncpus=4) <br>
install.packages("leiden", Ncpus=4) <br>
install.packages("Seurat", Ncpus=4) # igraph and leiden are required before installing Seurat <br>
packages <- c("dplyr", "magrittr", "data.table", "Matrix", "devtools", "RcppArmadillo", "Rcpp", "scales", "pheatmap", "gplots", "ggplot2", "cowplot", "tibble", "data.table", "remotes") # install the rest <br>
install.packages(packages, Ncpus=4) <br>
remotes::install_github("immunogenomics/harmony") # need batch correction package Harmony <br>
remotes::install_github("immunogenomics/presto") # need presto to speed up marker finding in Seurat <br>

5. RUN R SCRIPT <br>
#Create R script in sbatch file and upload to desired directory on Sherlock <br>
cd /scratch/users/<username>/<foldername> # change with desired directory <br>
ls # check that sbatch file is there <br>
sbatch R_main.sbatch <br>
squeue -u <username> # check status of run <br>
htop # see compute usage <br>
#Look at slurm.out file for report on run
