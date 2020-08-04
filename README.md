
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coMethTF

<!-- badges: start -->

[![codecov](https://codecov.io/gl/tiagochst/methtf/branch/%5Cx6d6173746572/graph/badge.svg?token=NESBYPVF64)](https://codecov.io/gl/tiagochst/methtf)
[![license](https://img.shields.io/badge/license-GPL%20\(%3E%3D%202\)-blue)]()
<!-- badges: end -->

`coMethTF` can be used to generate testable hypothesis on the
synergistic interaction of DMRs and TFs in gene regulation. `coMethTF`
can be used either to evaluate regulatory potentials of candidate
regions or to search for methylation coupled TF regulatory processes in
the entire genome.

## Installation

You can install the released version of coMethTF from GitHub with:

``` r
devtools::install_github("TransBioInfoLab/coMethTF")
```

## Example

This is a basic example which shows you how to use the package:

``` r
library(coMethTF)
#---------------------------------------
# Data input
#---------------------------------------
# 1) Gene expression matrix
# 2) DNA methylation
# With same column names
data("dna.met.chr21")
data("gene.exp.chr21")
all(colnames(dna.met.chr21) == colnames(gene.exp.chr21))
#> [1] TRUE

# Since we are working with regions we need to map our 450k array to regions
dnam.regions <- map_probes_to_regions(dna.met.chr21)
```

``` r
#---------------------------------------
# Mapping regions
#---------------------------------------
# For each region get target gene and predicted TF biding to the regions
# get_triplet incorporates two other functions:
# 1) get_region_target_gene
# 2) get_tf_in_region
triplet <- get_triplet(
    region = rownames(dnam.regions),
    motif.search.window.size = 50,
    motif.search.p.cutoff = 10^-3,
    target.method = "closest.gene",
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
    dnam = dnam.regions,
    exp = gene.exp.chr21
)
#> Removing genes with RNA expression equal to 0 for all samples from triplets
#> Removing triplet with no DNA methylation information for more than 25% of the samples
```

``` r
head(results)
#>                  regionID target_gene_name          target
#> 1 chr21:37258041-37258042         PPP1R2P2 ENSG00000234008
#> 2 chr21:37258041-37258042         PPP1R2P2 ENSG00000234008
#> 3 chr21:45553263-45553264         C21orf33 ENSG00000160221
#> 4 chr21:45973002-45973003        KRTAP10-2 ENSG00000205445
#> 5 chr21:32412214-32412215        KRTAP19-8 ENSG00000206102
#> 6 chr21:31873903-31873904        KRTAP19-5 ENSG00000186977
#>   TF_external_gene_name              TF TF_symbol target_symbol   pval_met
#> 1                 OLIG2 ENSG00000205927     OLIG2      PPP1R2P2 0.41736883
#> 2                 OLIG1 ENSG00000184221     OLIG1      PPP1R2P2 0.79914752
#> 3                  ETS2 ENSG00000157557      ETS2        GATD3A 0.04510031
#> 4                 OLIG1 ENSG00000184221     OLIG1     KRTAP10-2 0.23056884
#> 5                 BACH1 ENSG00000156273     BACH1     KRTAP19-8        NaN
#> 6                 BACH1 ENSG00000156273     BACH1     KRTAP19-5        NaN
#>   pval_rna.tf pval_met:rna.tf  estimate_met estimate_rna.tf estimate_met:rna.tf
#> 1  0.41569327      0.35888262 -2.638193e+01   -2.267866e-01        3.005770e+00
#> 2  0.77910914      0.75178558  2.210383e+00   -2.583537e-02        4.028465e-01
#> 3  0.04657082      0.04320122 -3.945954e+02   -1.122433e+01        1.865710e+01
#> 4  0.64005823      0.45547429  5.915981e-01    1.689216e-02       -9.129480e-02
#> 5         NaN             NaN -1.554984e-10   -4.194612e-08       -2.813108e-09
#> 6  0.78394162             NaN -3.353047e+00    1.900536e+01       -3.207120e+00
#>                       Model.interaction met.q4_minus_q1 quant_pval_metGrp
#> 1 Zero-inflated Negative Binomial Model      0.04123867               NaN
#> 2 Zero-inflated Negative Binomial Model      0.04123867        0.95384060
#> 3                   Robust Linear Model      0.18105392        0.03621953
#> 4 Zero-inflated Negative Binomial Model      0.21609131        0.22705049
#> 5 Zero-inflated Negative Binomial Model      0.01910411               NaN
#> 6 Zero-inflated Negative Binomial Model      0.01505539                NA
#>   quant_pval_rna.tf quant_pval_metGrp:rna.tf quant_estimate_metGrp
#> 1               NaN                      NaN         -3.500503e-08
#> 2        0.99999776               0.99999931          1.453565e+01
#> 3        0.04858627               0.03489433         -1.444184e+02
#> 4        0.98050022               0.98097527          3.342399e-01
#> 5               NaN                      NaN         -5.982758e-10
#> 6                NA                       NA         -5.816523e+01
#>   quant_estimate_rna.tf quant_estimate_metGrp:rna.tf
#> 1         -3.317330e-07                -2.542329e-07
#> 2         -1.910617e+00                -5.854055e-01
#> 3         -5.433864e+00                 6.762088e+00
#> 4         -1.488147e+00                 1.451887e+00
#> 5         -1.158967e-08                -1.158964e-08
#> 6         -1.895244e+01                 2.911858e+01
#>                          Model.quantile % 0 target genes (All samples)
#> 1 Zero-inflated Negative Binomial Model                        84.21 %
#> 2 Zero-inflated Negative Binomial Model                        84.21 %
#> 3                   Robust Linear Model                         2.63 %
#> 4 Zero-inflated Negative Binomial Model                        65.79 %
#> 5 Zero-inflated Negative Binomial Model                        97.37 %
#> 6 Zero-inflated Negative Binomial Model                        94.74 %
#>   % of 0 target genes (Q1 and Q4) Max_interaction_pval
#> 1                            95 %           0.35888262
#> 2                            95 %           0.99999931
#> 3                             0 %           0.04320122
#> 4                            65 %           0.98097527
#> 5                            95 %                   NA
#> 6                            90 %                   NA
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
#> [11] ExperimentHub_1.14.0              AnnotationHub_2.20.0             
#> [13] BiocFileCache_1.12.0              dbplyr_1.4.4                     
#> [15] BiocGenerics_0.34.0               coMethTF_0.1.0                   
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
#>  [27] backports_1.1.8               assertthat_0.2.1             
#>  [29] Matrix_1.2-18                 fastmap_1.0.1                
#>  [31] later_1.1.0.1                 htmltools_0.5.0              
#>  [33] prettyunits_1.1.1             tools_4.0.2                  
#>  [35] gtable_0.3.0                  glue_1.4.1                   
#>  [37] TFMPvalue_0.0.8               GenomeInfoDbData_1.2.3       
#>  [39] reshape2_1.4.4                dplyr_1.0.1                  
#>  [41] rappdirs_0.3.1                Rcpp_1.0.5                   
#>  [43] carData_3.0-4                 Biobase_2.48.0               
#>  [45] cellranger_1.1.0              vctrs_0.3.2                  
#>  [47] iterators_1.0.12              CNEr_1.24.0                  
#>  [49] xfun_0.16                     stringr_1.4.0                
#>  [51] openxlsx_4.1.5                mime_0.9                     
#>  [53] lifecycle_0.2.0               poweRlaw_0.70.6              
#>  [55] gtools_3.8.2                  rstatix_0.6.0                
#>  [57] XML_3.99-0.5                  zlibbioc_1.34.0              
#>  [59] MASS_7.3-51.6                 scales_1.1.1                 
#>  [61] hms_0.5.3                     promises_1.1.1               
#>  [63] SummarizedExperiment_1.18.2   JASPAR2020_0.99.10           
#>  [65] yaml_2.2.1                    curl_4.3                     
#>  [67] memoise_1.1.0                 ggplot2_3.3.2                
#>  [69] biomaRt_2.44.1                stringi_1.4.6                
#>  [71] RSQLite_2.2.0                 BiocVersion_3.11.1           
#>  [73] foreach_1.5.0                 caTools_1.18.0               
#>  [75] zip_2.0.4                     BiocParallel_1.22.0          
#>  [77] rlang_0.4.7                   pkgconfig_2.0.3              
#>  [79] matrixStats_0.56.0            bitops_1.0-6                 
#>  [81] pracma_2.2.9                  evaluate_0.14                
#>  [83] lattice_0.20-41               purrr_0.3.4                  
#>  [85] GenomicAlignments_1.24.0      bit_4.0.3                    
#>  [87] tidyselect_1.1.0              plyr_1.8.6                   
#>  [89] magrittr_1.5                  R6_2.4.1                     
#>  [91] generics_0.0.2                DelayedArray_0.14.1          
#>  [93] DBI_1.1.0                     pillar_1.4.6                 
#>  [95] haven_2.3.1                   foreign_0.8-80               
#>  [97] KEGGREST_1.28.0               abind_1.4-5                  
#>  [99] RCurl_1.98-1.2                tibble_3.0.3                 
#> [101] crayon_1.3.4                  car_3.0-8                    
#> [103] rmarkdown_2.3                 progress_1.2.2               
#> [105] TFBSTools_1.26.0              grid_4.0.2                   
#> [107] readxl_1.3.1                  data.table_1.13.0            
#> [109] blob_1.2.1                    forcats_0.5.0                
#> [111] digest_0.6.25                 xtable_1.8-4                 
#> [113] tidyr_1.1.1                   httpuv_1.5.4                 
#> [115] R.utils_2.9.2                 openssl_1.4.2                
#> [117] munsell_0.5.0                 DirichletMultinomial_1.30.0  
#> [119] motifmatchr_1.10.0            askpass_1.1
```
