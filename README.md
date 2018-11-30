# devnetwork

Prior to running scripts dealing with figures (src/analysis_R/submission_figures.R), multiple separate analyses need to be done**:
  1) differential expression (src/analysis_R/DEanalysis.R)
  2) phylostratigraphy (src/evolution/snakefile_blastPS)
  3) Aligning drosophila data (src/Drosophila_data/snakefile_calcExpression) and calculate sex-bias (src/Drosophila_data/dmel_sex.R)
  4) Aligning Jasper et al honey bee tissue data (src/jasper_data/snakefile_calcExpression) and calculate tau (src/jasper_data/jasper_tissue.R)
  5) Plaid clustering (src/analysis_R/plaid.R)
  6) evolutionary rate analyses (src/evolution)
  
**This list is detailed and implemented in the file complete.snake

Evolutionary data can be found here: https://www.dropbox.com/s/ug1jsyv5zm5xy2g/evolution_data.tar.gz?dl=0. 
Download, unpack, and put the evolution/ directory in data/. 

Software required:
