This folder contains the scripts necessary to run the main analyses in this project. <br>

### Contents

1. [`somatic_bacterial_analysis.Rmd`]
- The main R markdown script for performing the analyses in our project.
- Requires the input data "TotalTissue.combined.harmony_scaled.RData" that was generated from the Harmony integration on Stanford Sherlock.
- Requires the input data "steele_microbiome.RDS" which is provided on the SAHMI Github repository (Ghaddar et al., 2022.
- The inferCNV analysis requires the files in the `cnv_analysis` directory.
- Manual cell type annotation information was acquired in the output file "final_cluster_labeling.csv".
- Manual reactome pathway annotation information was acquired in the output file "reactome_corrected.csv".

### References

Ghaddar, B., Biswas, A., Harris, C., Omary, M. B., Carpizo, D. R., Blaser, M. J., & De, S. (2022). Tumor microbiome links cellular programs and immunity in pancreatic cancer. Cancer Cell, 40(10), 1240-1253.
