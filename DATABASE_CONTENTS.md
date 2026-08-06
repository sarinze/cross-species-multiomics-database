# Database Contents

The SQLite database `pseudomonas_new.db` contains the integrated multi-omics data generated during this MSc project.

The database brings together genomic, transcriptomic, ribosome-profiling, transposon-sequencing and ortholog-mapping information from *Pseudomonas aeruginosa* and *Klebsiella pneumoniae*.

## Main Data Categories

The database includes tables related to:

- Genome annotation
- RNA-seq datasets
- Ribo-seq datasets
- Tn-seq datasets
- Ortholog and locus-tag mapping
- Cross-species linker tables
- Integrated multi-omics matrices
- Quality-control outputs
- Normalised analytical data
- Hierarchical clustering outputs
- PCA outputs
- Conserved and divergent response analyses

## Main Integrated Outputs

Key integrated datasets include:

- Master ortholog backbone
- Master transcriptomics matrix
- Master functional genomics matrix
- Master multi-omics matrix
- Normalised multi-omics matrix
- Final clustering input
- Final cluster assignments
- Biological cluster characterisation
- Final gene catalogue
- Conserved and divergent response outputs

## Final Analytical Dataset

The completed integrated database contained 14,497 gene records.

Following quality-control and filtering, 6,970 genes with sufficient multi-omics coverage were retained for downstream analysis.

These genes were analysed using hierarchical clustering, principal component analysis and cross-species comparative response analysis.

## Accessing the Database

The database can be opened using:

- DB Browser for SQLite
- SQLite command line tools
- Python using the `sqlite3` package
- SQLite-compatible extensions in Visual Studio Code

The database is provided as:

`pseudomonas_new.db`