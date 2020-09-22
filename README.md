
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MethReg

<!-- badges: start -->

[![codecov](https://codecov.io/gl/tiagochst/methtf/branch/%5Cx6d6173746572/graph/badge.svg?token=NESBYPVF64)](https://codecov.io/gl/tiagochst/methtf)
[![license](https://img.shields.io/badge/license-GPL%20\(%3E%3D%202\)-blue)]()
<!-- badges: end -->

`MethReg` can be used to generate testable hypothesis on the synergistic
interaction of DMRs and TFs in gene regulation.

`MethReg` can be used either to evaluate regulatory potentials of
candidate regions or to search for methylation coupled TF regulatory
processes in the entire genome.

## Installation

You can install the MethReg from Bioconductor with:

``` r
BiocManager::install("MethReg")
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
dna.met.chr21 <- make_dnam_se(dna.met.chr21)
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
```

``` r
head(results)
#>                  regionID target_gene_name          target
#> 1 chr21:30372219-30372220          RPL23P2 ENSG00000176054
#> 2 chr21:30430511-30430512       AF129075.5 ENSG00000231125
#> 3 chr21:33109780-33109781       AP000255.6 ENSG00000273091
#> 4 chr21:40692859-40692860            BRWD1 ENSG00000185658
#> 5 chr21:43982646-43982647       AP001625.6 ENSG00000235772
#> 6 chr21:43983587-43983588       AP001625.6 ENSG00000235772
#>   TF_external_gene_name              TF TF_symbol target_symbol     met.IQR
#> 1                  ETS2 ENSG00000157557      ETS2       RPL23P2 0.182568764
#> 2                 BACH1 ENSG00000156273     BACH1    AF129075.2 0.310379208
#> 3                 GABPA ENSG00000154727     GABPA    AP000255.1 0.080160919
#> 4                 GABPA ENSG00000154727     GABPA         BRWD1 0.040638333
#> 5                 GABPA ENSG00000154727     GABPA    AP001625.2 0.008670064
#> 6                 GABPA ENSG00000154727     GABPA    AP001625.2 0.014175786
#>   quant_pval_metGrp quant_fdr_metGrp quant_pval_rna.tf quant_fdr_rna.tf
#> 1      3.953828e-05     0.0001186148      5.958024e-03     1.787407e-02
#> 2      7.437666e-01     0.7437665914      6.570509e-04     6.570509e-04
#> 3      2.208774e-03     0.0022087741      5.803942e-01     5.803942e-01
#> 4      3.823862e-01     0.3823861582      7.142112e-08     7.142112e-08
#> 5      4.434057e-01     0.4434057288      4.977765e-02     4.977765e-02
#> 6      5.619040e-01     0.9778889092      1.305771e-03     2.611541e-03
#>   quant_pval_metGrp:rna.tf quant_fdr_metGrp:rna.tf quant_estimate_metGrp
#> 1             3.768097e-05            0.0001130429            -83.041219
#> 2             7.945988e-01            0.7945988318             -3.533136
#> 3             2.305241e-03            0.0023052406            -30.858474
#> 4             4.999473e-01            0.4999473203             -3.386762
#> 5             4.663539e-01            0.4663538741            -10.643704
#> 6             5.234533e-01            0.9782246238             -5.525037
#>   quant_estimate_rna.tf quant_estimate_metGrp:rna.tf      Model.quantile
#> 1            -2.0395999                    3.8764266 Robust Linear Model
#> 2             1.3437155                    0.1642009 Robust Linear Model
#> 3             0.2627746                    1.8528803 Robust Linear Model
#> 4             1.3876278                    0.1522536 Robust Linear Model
#> 5             1.1156152                    0.5884195 Robust Linear Model
#> 6             1.2125286                    0.3594288 Robust Linear Model
#>   Wilcoxon_pval_target_q4met_vs_q1met Wilcoxon_pval_tf_q4met_vs_q1met
#> 1                          0.42735531                       0.5707504
#> 2                          0.47267559                       0.4726756
#> 3                          0.49466484                       0.9097219
#> 4                          0.03763531                       0.7913368
#> 5                          0.44952133                       0.2730363
#> 6                          1.00000000                       0.2730363
#>   % of 0 target genes (Q1 and Q4)
#> 1                             0 %
#> 2                             5 %
#> 3                            20 %
#> 4                             0 %
#> 5                            10 %
#> 6                             5 %
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
#> [11] ExperimentHub_1.15.3              AnnotationHub_2.20.2             
#> [13] BiocFileCache_1.12.1              dbplyr_1.4.4                     
#> [15] BiocGenerics_0.34.0               MethReg_0.99.11                  
#> 
#> loaded via a namespace (and not attached):
#>   [1] readxl_1.3.1                  backports_1.1.10             
#>   [3] plyr_1.8.6                    BiocParallel_1.22.0          
#>   [5] ggplot2_3.3.2                 TFBSTools_1.26.0             
#>   [7] digest_0.6.25                 foreach_1.5.0                
#>   [9] htmltools_0.5.0               GO.db_3.11.4                 
#>  [11] magrittr_1.5                  memoise_1.1.0                
#>  [13] doParallel_1.0.15             sfsmisc_1.1-7                
#>  [15] openxlsx_4.1.5                readr_1.3.1                  
#>  [17] annotate_1.66.0               matrixStats_0.56.0           
#>  [19] R.utils_2.10.1                JASPAR2020_0.99.10           
#>  [21] prettyunits_1.1.1             colorspace_1.4-1             
#>  [23] blob_1.2.1                    rappdirs_0.3.1               
#>  [25] haven_2.3.1                   xfun_0.17                    
#>  [27] dplyr_1.0.2                   crayon_1.3.4                 
#>  [29] RCurl_1.98-1.2                TFMPvalue_0.0.8              
#>  [31] iterators_1.0.12              glue_1.4.2                   
#>  [33] gtable_0.3.0                  sesame_1.6.0                 
#>  [35] zlibbioc_1.34.0               DelayedArray_0.14.1          
#>  [37] car_3.0-9                     wheatmap_0.1.0               
#>  [39] Rhdf5lib_1.10.1               HDF5Array_1.16.1             
#>  [41] abind_1.4-5                   scales_1.1.1                 
#>  [43] pscl_1.5.5                    DBI_1.1.0                    
#>  [45] rstatix_0.6.0                 Rcpp_1.0.5                   
#>  [47] xtable_1.8-4                  progress_1.2.2               
#>  [49] foreign_0.8-80                bit_4.0.4                    
#>  [51] preprocessCore_1.50.0         httr_1.4.2                   
#>  [53] RColorBrewer_1.1-2            ellipsis_0.3.1               
#>  [55] pkgconfig_2.0.3               XML_3.99-0.5                 
#>  [57] R.methodsS3_1.8.1             DNAcopy_1.62.0               
#>  [59] tidyselect_1.1.0              rlang_0.4.7                  
#>  [61] reshape2_1.4.4                later_1.1.0.1                
#>  [63] AnnotationDbi_1.50.3          munsell_0.5.0                
#>  [65] BiocVersion_3.11.1            cellranger_1.1.0             
#>  [67] tools_4.0.2                   DirichletMultinomial_1.30.0  
#>  [69] generics_0.0.2                RSQLite_2.2.0                
#>  [71] broom_0.7.0                   evaluate_0.14                
#>  [73] stringr_1.4.0                 fastmap_1.0.1                
#>  [75] yaml_2.2.1                    knitr_1.29                   
#>  [77] bit64_4.0.5                   zip_2.1.1                    
#>  [79] caTools_1.18.0                purrr_0.3.4                  
#>  [81] randomForest_4.6-14           KEGGREST_1.28.0              
#>  [83] mime_0.9                      R.oo_1.24.0                  
#>  [85] poweRlaw_0.70.6               pracma_2.2.9                 
#>  [87] compiler_4.0.2                curl_4.3                     
#>  [89] png_0.1-7                     interactiveDisplayBase_1.26.3
#>  [91] ggsignif_0.6.0                tibble_3.0.3                 
#>  [93] stringi_1.5.3                 forcats_0.5.0                
#>  [95] lattice_0.20-41               CNEr_1.24.0                  
#>  [97] Matrix_1.2-18                 vctrs_0.3.4                  
#>  [99] pillar_1.4.6                  lifecycle_0.2.0              
#> [101] BiocManager_1.30.10           data.table_1.13.0            
#> [103] bitops_1.0-6                  httpuv_1.5.4                 
#> [105] R6_2.4.1                      promises_1.1.1               
#> [107] rio_0.5.16                    codetools_0.2-16             
#> [109] MASS_7.3-53                   gtools_3.8.2                 
#> [111] assertthat_0.2.1              seqLogo_1.54.3               
#> [113] rhdf5_2.32.2                  SummarizedExperiment_1.18.2  
#> [115] GenomicAlignments_1.24.0      Rsamtools_2.4.0              
#> [117] GenomeInfoDbData_1.2.3        hms_0.5.3                    
#> [119] motifmatchr_1.10.0            grid_4.0.2                   
#> [121] tidyr_1.1.2                   rmarkdown_2.3                
#> [123] carData_3.0-4                 ggpubr_0.4.0                 
#> [125] Biobase_2.48.0                shiny_1.5.0
```
