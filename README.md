
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
triplet <- get_triplet(
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
```

``` r
head(results)
#>                  regionID target_gene_name          target
#> 1 chr21:18245122-18245123       AF212831.2 ENSG00000232886
#> 2 chr21:18245122-18245123       AF212831.2 ENSG00000232886
#> 3 chr21:18245122-18245123       AF212831.2 ENSG00000232886
#> 4 chr21:25562143-25562144       AP000470.2 ENSG00000224018
#> 5 chr21:25562370-25562371       AP000470.2 ENSG00000224018
#> 6 chr21:25562573-25562574       AP000470.2 ENSG00000224018
#>   TF_external_gene_name              TF TF_symbol target_symbol pval_met
#> 1                 OLIG2 ENSG00000205927     OLIG2    AF212831.1      NaN
#> 2                 OLIG1 ENSG00000184221     OLIG1    AF212831.1      NaN
#> 3                 GABPA ENSG00000154727     GABPA    AF212831.1      NaN
#> 4                 GABPA ENSG00000154727     GABPA    AP000470.1      NaN
#> 5                 OLIG2 ENSG00000205927     OLIG2    AP000470.1      NaN
#> 6                 GABPA ENSG00000154727     GABPA    AP000470.1      NaN
#>   pval_rna.tf pval_met:rna.tf  estimate_met estimate_rna.tf estimate_met:rna.tf
#> 1         NaN             NaN  9.558688e-01   -6.802002e-02        1.844318e-04
#> 2         NaN             NaN  8.669732e-01   -1.321329e-02        6.725235e-04
#> 3         NaN             NaN  1.113898e+00   -1.714299e-01        6.999789e-06
#> 4         NaN             NaN -3.294481e-10   -2.628426e-07       -5.099006e-09
#> 5         NaN             NaN -2.077405e-09   -4.182182e-07       -1.863170e-08
#> 6         NaN             NaN -1.297804e-09   -2.552366e-07       -1.900931e-08
#>                       Model.interaction met.q4_minus_q1 quant_pval_metGrp
#> 1 Zero-inflated Negative Binomial Model     0.165474975               NaN
#> 2 Zero-inflated Negative Binomial Model     0.165474975               NaN
#> 3 Zero-inflated Negative Binomial Model     0.165474975               NaN
#> 4 Zero-inflated Negative Binomial Model     0.003748681               NaN
#> 5 Zero-inflated Negative Binomial Model     0.012488103               NaN
#> 6 Zero-inflated Negative Binomial Model     0.058289052               NaN
#>   quant_pval_rna.tf quant_pval_metGrp:rna.tf quant_estimate_metGrp
#> 1               NaN                      NaN         -1.075548e-08
#> 2               NaN                      NaN         -1.704966e-08
#> 3               NaN                      NaN         -6.584724e-09
#> 4               NaN                      NaN         -7.239205e-09
#> 5               NaN                      NaN         -1.075846e-08
#> 6               NaN                      NaN         -7.239199e-09
#>   quant_estimate_rna.tf quant_estimate_metGrp:rna.tf
#> 1         -1.321712e-07                -5.132464e-08
#> 2         -2.148548e-07                -1.083325e-07
#> 3         -1.633070e-07                -1.212368e-07
#> 4         -1.757191e-07                -1.344448e-07
#> 5         -1.140343e-07                -8.120024e-08
#> 6         -1.742115e-07                -1.333088e-07
#>                          Model.quantile % 0 target genes (All samples)
#> 1 Zero-inflated Negative Binomial Model                        92.11 %
#> 2 Zero-inflated Negative Binomial Model                        92.11 %
#> 3 Zero-inflated Negative Binomial Model                        92.11 %
#> 4 Zero-inflated Negative Binomial Model                        97.37 %
#> 5 Zero-inflated Negative Binomial Model                        97.37 %
#> 6 Zero-inflated Negative Binomial Model                        97.37 %
#>   % of 0 target genes (Q1 and Q4) Max_interaction_pval
#> 1                            95 %                   NA
#> 2                            95 %                   NA
#> 3                            95 %                   NA
#> 4                            95 %                   NA
#> 5                            95 %                   NA
#> 6                            95 %                   NA
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
#>  [47] iterators_1.0.12              xfun_0.16                    
#>  [49] CNEr_1.24.0                   stringr_1.4.0                
#>  [51] openxlsx_4.1.5                mime_0.9                     
#>  [53] lifecycle_0.2.0               poweRlaw_0.70.6              
#>  [55] gtools_3.8.2                  rstatix_0.6.0                
#>  [57] XML_3.99-0.5                  zlibbioc_1.34.0              
#>  [59] MASS_7.3-51.6                 scales_1.1.1                 
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
#>  [85] bit_4.0.3                     tidyselect_1.1.0             
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
