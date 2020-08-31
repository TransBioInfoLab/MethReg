
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MethReg

<!-- badges: start -->

[![codecov](https://codecov.io/gl/tiagochst/methtf/branch/%5Cx6d6173746572/graph/badge.svg?token=NESBYPVF64)](https://codecov.io/gl/tiagochst/methtf)
[![license](https://img.shields.io/badge/license-GPL%20\(%3E%3D%202\)-blue)]()
<!-- badges: end -->

`MethReg` can be used to generate testable hypothesis on the synergistic
interaction of DMRs and TFs in gene regulation. `MethReg` can be used
either to evaluate regulatory potentials of candidate regions or to
search for methylation coupled TF regulatory processes in the entire
genome.

## Installation

You can install the released version of MethReg from GitHub with:

``` r
devtools::install_github("TransBioInfoLab/MethReg")
```

## Example

This is a basic example which shows you how to use the package:

``` r
library(MethReg)
#---------------------------------------
# Data input
#---------------------------------------
# 1) Gene expression matrix
# 2) DNA methylation
# With same column names
data("dna.met.chr21")
data("gene.exp.chr21.log2")
all(colnames(dna.met.chr21) == colnames(gene.exp.chr21.log2))
#> [1] TRUE

# Since we are working with regions we need to map our 450k array to regions
dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
#> o Creating a SummarizedExperiment from DNA methylation input
#> oo Fetching probes metadata
#> oo Removing masked probes
#> oo Preparing SummarizedExperiment object
```

``` r
#---------------------------------------
# Mapping regions
#---------------------------------------
# For each region get target gene and predicted TF biding to the regions
# get_triplet incorporates two other functions:
# 1) get_region_target_gene
# 2) get_tf_in_region
triplet <- create_triplet_distance_based(
    region = rownames(dna.met.chr21),
    motif.search.window.size = 50,
    motif.search.p.cutoff = 10^-3,
    target.method = "genes.promoter.overlap",
    genome = "hg19",
    cores = 1
)
#> Finding target genes
#> Mapping regions to the closest gene
#> Looking for TFBS
#> 
#> Evaluating 554 JASPAR Human TF motifs
#> This may take a while...
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> Preparing output
#> Joining, by = "regionID"
```

``` r
#---------------------------------------
# Evaluate two models: 
#---------------------------------------
# 1) target gene ~ TF + DNAm + TF * DNAm
# 2) target gene ~ TF + DNAm_group + TF * DNAm_group 
# where DNAm_group is a binary indicator if the sample belongs to: Q4 or Q1
results <- interaction_model(
    triplet = triplet, 
    dnam = dna.met.chr21,
    exp = gene.exp.chr21.log2
)
#> Removing genes with RNA expression equal to 0 for all samples from triplets
#> Removing triplet with no DNA methylation information for more than 25% of the samples
#> Filtering results to have interaction, TF or DNAm significant
#> Filtering results to wilcoxon test TF Q1 vs Q4 not significant
```

``` r
head(results)
#>                  regionID target_gene_name          target
#> 1 chr21:30372219-30372220          RPL23P2 ENSG00000176054
#> 2 chr21:30430511-30430512       AF129075.5 ENSG00000231125
#> 3 chr21:33109780-33109781       AP000255.6 ENSG00000273091
#> 4 chr21:40692859-40692860            BRWD1 ENSG00000185658
#> 5 chr21:43983587-43983588       AP001625.6 ENSG00000235772
#> 6 chr21:45468403-45468404          H2AFZP1 ENSG00000213440
#>   TF_external_gene_name              TF TF_symbol target_symbol   pval_met
#> 1                  ETS2 ENSG00000157557      ETS2       RPL23P2 0.01373419
#> 2                 BACH1 ENSG00000156273     BACH1    AF129075.2 0.61134838
#> 3                 GABPA ENSG00000154727     GABPA    AP000255.1 0.00000000
#> 4                 GABPA ENSG00000154727     GABPA         BRWD1 0.67481637
#> 5                 GABPA ENSG00000154727     GABPA    AP001625.2 0.61368294
#> 6                PKNOX1 ENSG00000160199    PKNOX1        H2AZP1 0.27934091
#>    pval_rna.tf pval_met:rna.tf estimate_met estimate_rna.tf estimate_met:rna.tf
#> 1 8.927789e-02      0.01263826   -94.880189      -1.8366801           4.5153185
#> 2 9.635728e-03      0.67327726    -9.775541       1.5154069           0.4727370
#> 3 2.341924e-08      0.00000000   258.708769       1.6970768         -15.4781016
#> 4 3.783765e-02      0.62747811     6.141972       1.6780249          -0.4139619
#> 5 6.170294e-01      0.57380835  -167.791692      -9.3709854          10.9969876
#> 6 8.209816e-01      0.28503982   -91.586347      -0.9169839           5.6172434
#>     Model.interaction met.q4_minus_q1 quant_pval_metGrp quant_pval_rna.tf
#> 1 Robust Linear Model      0.18256876      3.953828e-05      5.958024e-03
#> 2 Robust Linear Model      0.31037921      7.437666e-01      6.570509e-04
#> 3 Robust Linear Model      0.08016092      2.208774e-03      5.803942e-01
#> 4 Robust Linear Model      0.04063833      3.823862e-01      7.142112e-08
#> 5 Robust Linear Model      0.01417579      5.619040e-01      1.305771e-03
#> 6 Robust Linear Model      0.14299567      1.524877e-03      3.270442e-01
#>   quant_pval_metGrp:rna.tf quant_estimate_metGrp quant_estimate_rna.tf
#> 1             3.768097e-05            -83.041219            -2.0395999
#> 2             7.945988e-01             -3.533136             1.3437155
#> 3             2.305241e-03            -30.858474             0.2627746
#> 4             4.999473e-01             -3.386762             1.3876278
#> 5             5.234533e-01             -5.525037             1.2125286
#> 6             1.727783e-03           -182.019620             1.8203862
#>   quant_estimate_metGrp:rna.tf      Model.quantile Wilcoxon_pval_tf_q4_vs_q1
#> 1                    3.8764266 Robust Linear Model                 0.5787417
#> 2                    0.1642009 Robust Linear Model                 0.4812509
#> 3                    1.8528803 Robust Linear Model                 0.9117972
#> 4                    0.1522536 Robust Linear Model                 0.7959363
#> 5                    0.3594288 Robust Linear Model                 0.2798610
#> 6                   11.0029788 Robust Linear Model                 0.5288489
#>   % 0 target genes (All samples) % of 0 target genes (Q1 and Q4)
#> 1                         2.63 %                             0 %
#> 2                         5.26 %                             5 %
#> 3                        21.05 %                            20 %
#> 4                            0 %                             0 %
#> 5                         7.89 %                             5 %
#> 6                        10.53 %                            10 %
#>   Max_interaction_pval  fdr_met   fdr_rna.tf fdr_met:rna.tf quant_fdr_metGrp
#> 1          0.012638262 0.576836 7.055858e-01      0.4967681      0.002807218
#> 2          0.794598832 1.000000 2.665885e-01      0.9993441      0.999999995
#> 3          0.002305241 0.000000 1.943797e-06      0.0000000      0.039205740
#> 4          0.627478107 1.000000 6.281050e-01      0.9993441      0.999999995
#> 5          0.573808352 1.000000 1.000000e+00      0.9993441      0.999999995
#> 6          0.285039821 1.000000 1.000000e+00      0.9993441      0.039205740
#>   quant_fdr_rna.tf quant_fdr_metGrp:rna.tf Wilcoxon_fdr_tf_q4_vs_q1
#> 1     8.016700e-02             0.002675349                0.7885356
#> 2     1.598824e-02             0.999999995                0.7388219
#> 3     9.999978e-01             0.040918021                0.9376028
#> 4     5.213742e-06             0.999999995                0.9119969
#> 5     2.383031e-02             0.999999995                0.5546336
#> 6     9.947595e-01             0.040918021                0.7685937
```

# Session information

``` r
sessionInfo()
#> R version 4.0.2 (2020-06-22)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Catalina 10.15.6
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] BSgenome.Hsapiens.UCSC.hg19_1.4.3 BSgenome_1.56.0                  
#>  [3] rtracklayer_1.48.0                Biostrings_2.56.0                
#>  [5] XVector_0.28.0                    GenomicRanges_1.40.0             
#>  [7] GenomeInfoDb_1.24.2               IRanges_2.22.2                   
#>  [9] S4Vectors_0.26.1                  sesameData_1.6.0                 
#> [11] ExperimentHub_1.14.2              AnnotationHub_2.20.2             
#> [13] BiocFileCache_1.12.1              dbplyr_1.4.4                     
#> [15] BiocGenerics_0.34.0               MethReg_0.1.0                    
#> 
#> loaded via a namespace (and not attached):
#>   [1] colorspace_1.4-1              ggsignif_0.6.0               
#>   [3] ellipsis_0.3.1                rio_0.5.16                   
#>   [5] ggpubr_0.4.0                  bit64_4.0.2                  
#>   [7] interactiveDisplayBase_1.26.3 AnnotationDbi_1.50.3         
#>   [9] codetools_0.2-16              R.methodsS3_1.8.0            
#>  [11] doParallel_1.0.15             pscl_1.5.5                   
#>  [13] knitr_1.29                    Rsamtools_2.4.0              
#>  [15] seqLogo_1.54.3                annotate_1.66.0              
#>  [17] broom_0.7.0                   GO.db_3.11.4                 
#>  [19] png_0.1-7                     R.oo_1.23.0                  
#>  [21] sfsmisc_1.1-7                 shiny_1.5.0                  
#>  [23] BiocManager_1.30.10           readr_1.3.1                  
#>  [25] compiler_4.0.2                httr_1.4.2                   
#>  [27] backports_1.1.9               assertthat_0.2.1             
#>  [29] Matrix_1.2-18                 fastmap_1.0.1                
#>  [31] later_1.1.0.1                 htmltools_0.5.0              
#>  [33] prettyunits_1.1.1             tools_4.0.2                  
#>  [35] gtable_0.3.0                  glue_1.4.1                   
#>  [37] TFMPvalue_0.0.8               GenomeInfoDbData_1.2.3       
#>  [39] reshape2_1.4.4                dplyr_1.0.2                  
#>  [41] rappdirs_0.3.1                Rcpp_1.0.5                   
#>  [43] carData_3.0-4                 Biobase_2.48.0               
#>  [45] cellranger_1.1.0              vctrs_0.3.2                  
#>  [47] iterators_1.0.12              xfun_0.16                    
#>  [49] CNEr_1.24.0                   stringr_1.4.0                
#>  [51] openxlsx_4.1.5                mime_0.9                     
#>  [53] lifecycle_0.2.0               poweRlaw_0.70.6              
#>  [55] gtools_3.8.2                  rstatix_0.6.0                
#>  [57] XML_3.99-0.5                  zlibbioc_1.34.0              
#>  [59] MASS_7.3-52                   scales_1.1.1                 
#>  [61] hms_0.5.3                     promises_1.1.1               
#>  [63] SummarizedExperiment_1.18.2   JASPAR2020_0.99.10           
#>  [65] yaml_2.2.1                    curl_4.3                     
#>  [67] memoise_1.1.0                 ggplot2_3.3.2                
#>  [69] stringi_1.4.6                 RSQLite_2.2.0                
#>  [71] BiocVersion_3.11.1            foreach_1.5.0                
#>  [73] caTools_1.18.0                zip_2.1.0                    
#>  [75] BiocParallel_1.22.0           rlang_0.4.7                  
#>  [77] pkgconfig_2.0.3               matrixStats_0.56.0           
#>  [79] bitops_1.0-6                  pracma_2.2.9                 
#>  [81] evaluate_0.14                 lattice_0.20-41              
#>  [83] purrr_0.3.4                   GenomicAlignments_1.24.0     
#>  [85] bit_4.0.4                     tidyselect_1.1.0             
#>  [87] plyr_1.8.6                    magrittr_1.5                 
#>  [89] R6_2.4.1                      generics_0.0.2               
#>  [91] DelayedArray_0.14.1           DBI_1.1.0                    
#>  [93] pillar_1.4.6                  haven_2.3.1                  
#>  [95] foreign_0.8-80                KEGGREST_1.28.0              
#>  [97] abind_1.4-5                   RCurl_1.98-1.2               
#>  [99] tibble_3.0.3                  crayon_1.3.4                 
#> [101] car_3.0-9                     rmarkdown_2.3                
#> [103] progress_1.2.2                TFBSTools_1.26.0             
#> [105] grid_4.0.2                    readxl_1.3.1                 
#> [107] data.table_1.13.0             blob_1.2.1                   
#> [109] forcats_0.5.0                 digest_0.6.25                
#> [111] xtable_1.8-4                  tidyr_1.1.1                  
#> [113] httpuv_1.5.4                  R.utils_2.9.2                
#> [115] munsell_0.5.0                 DirichletMultinomial_1.30.0  
#> [117] motifmatchr_1.10.0
```
