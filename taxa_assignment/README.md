This folder contains the scripts necessary to assign taxa labels to cell barcodes in the PDAC data. <br>

### Contents

1. [`barcode_overlap_matrix.ipynb`]
- Jupyter notebook for calculating the number of overlapping barcodes between each anonymized PDAC sample and the processed SAHMI samples.
- Output data in "cell_barcodes.xlsx".
2. [`celltype_x_taxa_counts_matrix.ipynb`]
- Jupyter notebook for assigning taxa counts to cell types based on their matching cell barcodes, and for calculating Shannon Diversity Index.
- Output data in "taxa_x_celltypes.xlsx"
