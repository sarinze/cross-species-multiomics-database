# Cross-Species Multi-Omics Database

This repository contains the integrated database and Python workflows developed as part of my MSc research project in Biomedical and Molecular Sciences at the University of Dundee. UK.

The project focused on integrating and comparing publicly available multi-omics datasets from *Pseudomonas aeruginosa* and *Klebsiella pneumoniae* to investigate bacterial responses to antibiotic stress.

## Project Overview

The project brought together multiple biological data types that had originally been generated using different experimental platforms. These included:

- genome annotation data
- RNA sequencing
- ribosome profiling
- transposon sequencing
- ortholog mapping across strains and species

The datasets were harmonised into a single SQLite database and analysed using Python-based workflows.

## Main Database

`pseudomonas_new.db`

This SQLite database contains the integrated multi-omics resource developed during the project.

The database includes information from several *Pseudomonas aeruginosa* and *Klebsiella pneumoniae* strains and was used to support cross-species comparison, quality control, clustering, principal component analysis and conserved/divergent response analysis.

## Analysis Workflow

The numbered Python scripts reflect the main stages of the analytical workflow, including:

1. construction of the ortholog backbone
2. integration of transcriptomic and functional genomic datasets
3. addition of *Klebsiella pneumoniae* datasets
4. quality control and data filtering
5. data normalisation
6. hierarchical clustering
7. cluster evaluation and characterisation
8. heatmap generation
9. principal component analysis
10. conserved and divergent response analysis
11. import of final analytical outputs into the database

Additional scripts were used for data inspection, identifier mapping, validation and troubleshooting during database development.

## Final Analytical Dataset

The original integrated database contained 14,497 gene records.

Following quality-control and filtering steps, the final analytical dataset contained 6,970 genes with sufficient multi-omics coverage for downstream analysis.

Hierarchical clustering was performed using Euclidean distance and Ward linkage, and a 20-cluster solution was retained for biological interpretation.

Principal component analysis was also used to examine major patterns of variation across the integrated dataset.

## Cross-Species Analysis

Comparative analysis was used to investigate conserved and divergent antibiotic-response patterns between *Pseudomonas aeruginosa* and *Klebsiella pneumoniae*.

The project included both:

- a strict comparison based on directly comparable colistin-response measurements
- a broader comparison incorporating additional available measurements

## Tools and Technologies

The project was developed using:

- Python
- pandas
- NumPy
- SciPy
- scikit-learn
- Matplotlib
- SQLite
- NCBI BLAST+
- Visual Studio Code

## Repository Contents

- `pseudomonas_new.db` — integrated SQLite multi-omics database
- numbered `.py` files — main data integration and analysis workflow
- `import_*.py` — dataset import scripts
- `create_*.py` — ortholog and identifier-linking workflows
- `inspect_*.py` and `check_*.py` — data inspection and validation scripts

## Data Sources

This project was based on publicly available genomic, transcriptomic, ribosome-profiling and transposon-sequencing datasets.

Original datasets remain subject to the terms, licences and citation requirements of their respective publications and repositories.

## Purpose of This Repository

This repository is intended to demonstrate the computational workflow developed during my MSc research, including multi-omics data integration, database construction, comparative genomics, quality control, statistical analysis and biological interpretation.

It also serves as a reproducible record of the code used to construct and analyse the integrated database.

## Author

Chiamaka Sandra Arinze  
MSc Biomedical and Molecular Sciences  
University of Dundee