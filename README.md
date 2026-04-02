TREM1 enables fatty acid uptake to drive lipid-associated macrophage differentiation in colorectal cancer

This repository contains the complete analysis pipeline for the multi-omics study investigating how TREM1 enables fatty acid uptake to drive lipid-associated macrophage (LAM) differentiation in colorectal cancer (CRC).

📋 Overview
This project integrates three major omics approaches:
Single-cell RNA sequencing (scRNA-seq): Characterization of tumor microenvironment and LAM identification
ATAC-seq: Chromatin accessibility profiling of macrophage subpopulations
Spatial transcriptomics: Spatial mapping of TREM1+ macrophages and lipid metabolism
🗂️ Repository Structure
plain
复制
├── scRNA-seq/
│   ├── 01_preprocessing.R              # Cell Ranger output processing & QC
│   ├── 02_clustering_annotation.R      # Clustering, cell type annotation
│   ├── 03_lam_identification.R         # Lipid-associated macrophage identification
│   ├── 04_trem1_analysis.R             # TREM1+ macrophage characterization
│   ├── 05_trajectory_analysis.R        # Pseudotime trajectory (Monocle3/Slingshot)
│   ├── 06_cellchat_analysis.R          # Cell-cell communication analysis
│   └── utils/
│       └── scRNA_functions.R           # Custom helper functions
│
├── ATAC-seq/
│   ├── 01_preprocessing.sh             # FastQC, trimming, alignment (Bowtie2)
│   ├── 02_peak_calling.sh              # MACS2 peak calling
│   ├── 03_differential_accessibility.R # DiffBind/edgeR differential analysis
│   ├── 04_motif_analysis.R             # HOMER/ChIPseeker motif enrichment
│   ├── 05_integration_with_scRNA.R     # scRNA-ATAC integration (Signac/Seurat)
│   └── utils/
│       └── atac_config.yaml            # Configuration parameters
│
├── spatial_transcriptomics/
│   ├── 01_spaceranger_processing.sh    # 10x Visium data processing
│   ├── 02_qc_normalization.R           # Spot QC, SCTransform normalization
│   ├── 03_spatial_clustering.R         # BayesSpace/Seurat spatial clustering
│   ├── 04_deconvolution.R              # RCTD/SPOTlight cell type deconvolution
│   ├── 05_trem1_spatial_mapping.R      # TREM1+ LAM spatial distribution
│   ├── 06_neighborhood_analysis.R      # Spatial neighborhood analysis (Giotto)
│   └── utils/
│       └── spatial_viz_functions.R
│
├── data/
│   ├── metadata/
│   │   ├── sample_info.xlsx            # Sample metadata and clinical info
│   │   └── cell_type_markers.xlsx      # Marker gene lists for annotation
│   └── reference/
│       └── refdata-gex-GRCh38-2020-A/  # Cell Ranger reference (user-provided)
│
├── results/
│   ├── figures/                        # Publication-ready figures
│   └── tables/                         # Supplementary tables
│
├── docker/
│   └── Dockerfile                      # Containerized environment
│
└── renv/                               # R environment lock files
🚀 Quick Start
Prerequisites
Option 1: Using Docker (Recommended for Reproducibility)
bash
复制
# Build the container
docker build -t trem1-crc-analysis .

# Run interactive session
docker run -it -v $(pwd):/workspace trem1-crc-analysis
Option 2: Local Installation
System Requirements:
R >= 4.2.0
Python >= 3.9
64GB RAM recommended (scRNA-seq), 32GB (ATAC-seq), 16GB (Spatial)
R Dependencies:
r
复制
# Install core packages
install.packages(c("Seurat", "Signac", "monocle3", "CellChat", "DESeq2", "edgeR"))

# Install from GitHub/Bioconductor
remotes::install_github("satijalab/seurat")
remotes::install_github("cole-trapnell-lab/monocle3")
BiocManager::install(c("GenomicRanges", "rtracklayer", "ChIPseeker"))
Python Dependencies:
bash
复制
pip install scanpy squidpy cell2location
Data Download
Due to data sharing restrictions, raw sequencing data are available through:
GEO Accession: GSEXXXXXX (scRNA-seq & ATAC-seq)
ArrayExpress: E-MTAB-XXXXX (Spatial transcriptomics)
Zenodo DOI: 10.5281/zenodo.XXXXXXX (Processed Seurat objects)
Note: Users must download data independently and place in data/raw/ following the directory structure in docs/DATA_STRUCTURE.md.
🔬 Analysis Workflows
1. Single-Cell RNA-seq Analysis
bash
复制
# Navigate to scRNA-seq directory
cd scRNA-seq/

# Run complete pipeline
Rscript 01_preprocessing.R --sample_list ../data/metadata/samples.csv
Rscript 02_clustering_annotation.R --input ../data/processed/seurat_raw.rds
Rscript 03_lam_identification.R --markers ../data/metadata/lam_markers.csv
Rscript 04_trem1_analysis.R
Rscript 05_trajectory_analysis.R
Rscript 06_cellchat_analysis.R
Key Analyses:
QC & Filtering: nFeature > 200, nFeature < 7500, percent.mt < 15%
Doublet Removal: Scrublet (Python) or DoubletFinder (R)
Integration: Harmony for batch correction across samples
LAM Identification: Based on lipid metabolism gene signatures (LPL, FABP5, CD36, TREM1)
Trajectory: Monocle3 pseudotime from monocytes → LAMs
CellChat: TREM1-TYROBP signaling pathway analysis
2. ATAC-seq Analysis
bash
复制
cd ATAC-seq/

# Step 1: Preprocessing (bash)
bash 01_preprocessing.sh --fastq_dir ../data/raw/atac/

# Step 2-5: Downstream analysis (R)
Rscript 03_differential_accessibility.R \
    --peaks ../data/processed/consensus_peaks.bed \
    --metadata ../data/metadata/atac_samples.csv
Key Analyses:
Peak Calling: MACS2 with --nomodel --shift -100 --extsize 200
Differential Accessibility: DiffBind (DESeq2) for TREM1+ vs TREM1- LAMs
Motif Analysis: HOMER for TF motifs in TREM1-associated peaks
Footprinting: HINT-ATAC for TF binding dynamics
Multi-omics Integration: Link peaks to genes using Signac
3. Spatial Transcriptomics Analysis
bash
复制
cd spatial_transcriptomics/

# Process all samples
bash 01_spaceranger_processing.sh --samples CRC1,CRC2,CRC3,CRC4

# R analysis pipeline
Rscript 02_qc_normalization.R
Rscript 03_spatial_clustering.R --method BayesSpace
Rscript 04_deconvolution.R --sc_ref ../scRNA-seq/results/seurat_final.rds
Rscript 05_trem1_spatial_mapping.R
Rscript 06_neighborhood_analysis.R
Key Analyses:
Platform: 10x Genomics Visium
Deconvolution: RCTD using scRNA-seq as reference
Spatial Domains: BayesSpace clustering (k=10)
Niche Analysis: Cell neighborhood analysis with Giotto
Lipid Metabolism: Spatial mapping of fatty acid uptake genes
📊 Key Findings & Figures
表格
Figure	Description	Script
Fig 1	scRNA-seq landscape of CRC TME	scRNA-seq/02_clustering_annotation.R
Fig 2	LAM identification and TREM1 characterization	scRNA-seq/03_lam_identification.R
Fig 3	TREM1 drives lipid uptake via CD36	scRNA-seq/04_trem1_analysis.R
Fig 4	ATAC-seq reveals TREM1 regulatory network	ATAC-seq/03_differential_accessibility.R
Fig 5	Spatial distribution of TREM1+ LAMs	spatial_transcriptomics/05_trem1_spatial_mapping.R
Fig 6	TREM1-TYROBP signaling in tumor progression	scRNA-seq/06_cellchat_analysis.R
🧬 Key Gene Signatures
Lipid-Associated Macrophage (LAM) Signature:
r
复制
lam_markers <- c("TREM1", "TREM2", "CD9", "FABP5", "LPL", "CD36", 
                 "SPP1", "APOE", "C1QA", "C1QB", "GPNMB", "LGALS3")
Fatty Acid Uptake Pathway:
r
复制
fatty_acid_genes <- c("CD36", "FABP1", "FABP4", "FABP5", "SLC27A1", 
                      "ACSL1", "ACSL4", "FASN", "SCD1")
🛠️ Troubleshooting
Common Issues
Issue 1: Memory errors in Seurat integration
r
复制
# Solution: Use sketch-based integration or increase memory
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB
Issue 2: ATAC-seq fragment file reading
bash
复制
# Ensure tabix index exists
tabix -p bed fragments.tsv.gz
Issue 3: Spatial deconvolution fails
r
复制
# RCTD requires matching gene sets
# Run: CheckGeneSymbols() before deconvolution
📚 Citation
If you use this code or data, please cite:
bibtex
复制
@article{trem1_crc_2024,
  title={TREM1 enables fatty acid uptake to drive lipid-associated macrophage differentiation in colorectal cancer},
  author={[Author List]},
  journal={[Journal Name]},
  year={2024},
  doi={[DOI]}
}
📧 Contact
For technical questions about the analysis pipeline, please open an Issue.
For questions regarding the biological findings, contact: [corresponding.author@institution.edu]
📄 License
This project is licensed under the MIT License - see LICENSE file.
🙏 Acknowledgments
10x Genomics for spatial transcriptomics tools
Satija Lab for Seurat and Signac
Trapnell Lab for Monocle3
Greenleaf Lab for chromVAR
Last Updated: April 2026
