# devnetwork

Prior to running scripts dealing with figures, multiple separate analyses need to be done**:
  1) differential expression (src/analysis_R/DEanalysis.R)
  2) phylostratigraphy (phylostratigraphy/src/snakefile_blastPS)
  3) Aligning drosophila data (src/Drosophila_data/snakefile_calcExpression) and called sex-bias (src/Drosophila_data/dmel_sex.R)
  4) Aligning Jasper et al honey bee tissue data (src/jasper_data/snakefile_calcExpression) and calculate tau (src/jasper_data/jasper_tissue.R)
  5) Plaid clustering (src/analysis_R/plaid.R)
  
**This list is detailed and implemented in the file complete.snake

For the moment, data used for phylostratigraphy and for aligning reads (i.e. genomes and fastq) are not included.
