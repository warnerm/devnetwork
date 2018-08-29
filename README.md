# devnetwork

Prior to running scripts dealing with figures, multiple separate analyses need to be done**:
  1) differential expression (src/analysis_R/DEanalysis.R)
  2) phylostratigraphy (phylostratigraphy/src/snakefile_blastPS)
  3) Aligning drosophila data (src/Drosophila_data/snakefile_calcExpression) and called sex-bias (src/Drosophila_data/dmel_sex.R)
  4) Aligning Jasper et al honey bee tissue data (src/jasper_data/snakefile_calcExpression) and calculate tau (src/jasper_data/jasper_tissue.R)
  5) Plaid clustering (src/analysis_R/plaid.R)
  
**This list is detailed and implemented in the file complete.snake

For the moment, data used for phylostratigraphy and for aligning reads (i.e. genomes and fastq) are not included.

_Further descriptions_

The file "results/DEtests.RData" contains a number of lists of differential expression results as well as summarized data

Lists of DE results:
  antTests, beeTests, antTests_oneLarv, beeTests_oneLarv:
    These four lists are lists of differential expression results (table containing logFC, FDR, etc for all genes) at each         stage/tissue separately. The lists are named by tissue/stage. the "oneLarv" lists contain results for genes differentially     expressed between queen and worker-destined larvae overall (i.e. ~caste + stage + colony for all larval samples)
    Negative logFC indicates queen-biased
    
  antSocial,beeSocial:
    These are lists of differential expression results for nurses vs foragers at each tissue. 
    Negative logFC indicates nurse-biased

  antDevel2,beeDevel2:
    These are results for differential expression across larval stages, using all larvae with the model ~stage+colony+caste

Summarized Data:
  ant_sexDE,bee_sexDE,ant_VM,bee_VM:
    These are lists of genes differentially expressed between queens and males (sexDE), and queens and gynes (VM). 
    Negative logFC indicates queen-biased
    
  antRes,antRes_allstage,antSocRes,beeRes,beeRes_allstage,beeSocRes:
    These are some lists of dataframes I put together synthesizing DE results, where each list corresponds to the DE result       list (i.e. antRes is antTests_oneLarv, antRes_allstage is antTests, and antSocRes is antSocial)
    The first dataframe in each list is logFC, second is differential expression calls (given FDR < 0.1), and the the third is     the FDR for each comparison
    
  
