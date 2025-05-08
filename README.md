Pan-cancer mechanobiology analysis (30 cell lines spanning 8 tissue types on 7 different substrates):

**Paper reference to be added**

- data: contains the physical feature data used in the study
- Figures:
  - main and supplementary figures in the paper
  - cell trajectory plots, for which the cell tracking data is publically available at 10.6084/m9.figshare.28916828. Note that data/motility.tsv is also extracted from this cell tracking data.
- substrate_sensitivity: analysis for how physical properties of cell lines change in response to mechanosensitive and chemosensitive changes in ECM conditions (Results Section I of the paper)
- compare_normal_cancer_cell_lines: within tissue type comaprisons between normal and cancer cell lines on each of the substrates (Results Section II of the paper)
- consensus_clustering: unsupervised machine learning identify phenotypic classes for each physical behavior (Results Section III of the paper)
- migratory_persistence: using cell tracking data (publically available at 10.6084/m9.figshare.28916828) to calculate directional autocorrelation and then estimate decorrelation time as the measure of migratory persistence for cell lines on different substrates (Supp. Fig. 15-21).
- Supplementary_Information: LaTex files for generating the supplementary file associated with the paper

### Version of R packages used
tidyverse: 2.0.0\
foreach: 1.5.2\
doParallel: 1.0.17\
RColorBrewer: 1.1.3\
ComplexHeatmap: 2.20.0\
circlize: 0.4.16\
ggsignif: 0.6.4\
maotai: 0.2.5\
geostats: 1.6\
ConsensusClusterPlus: 1.68.0\
fpc: 2.2.12\
pracma: 2.4.4\
cowplot: 1.1.3\
corrplot: 0.92\
Hmisc: 5.1.2\
kableExtra: 1.4.0\
patchwork: 1.2.0\
plotrix: 3.8.4\
celltrackR: 1.2.1
