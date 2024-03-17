This folder contains the scripts necessary to load the single-cell PDAC data and perform batch integration using Harmony. <br>

### Contents

1. [`R_pdac.R`]
- R script from Steele et al. (2020) for loading and combining the filtered feature barcode matrices of the 17 PDAC scRNA-seq sample.
- The final output of the script is a Seurat object "TotalTissue.combined.harmony.RData" which is loaded into "R_harmony.sbatch.
- The next steps are very compute heavy and need to be run using the Stanford Sherlock computing cluster.
2. [`Running Seurat on Stanford Sherlock.md`]
- Instructions on how to set up an R environment for Seurat and submit a Slurm job using Sherlock.
3. [`R_harmony.sbatch`]
- Bash script for Slurm job submission in Sherlock.
- The final output is a Harmony-integrated Seurat object "TotalTissue.combined.harmony_scaled.RData" which is used for all downstream analyses.

### References

Steele, N. G., Carpenter, E. S., Kemp, S. B., Sirihorachai, V. R., The, S., Delrosario, L., ... & Pasca di Magliano, M. (2020). Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer. Nature Cancer, 1(11), 1097-1112.
