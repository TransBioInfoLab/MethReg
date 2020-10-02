test_that("get_region_target_gene works", {
   regions.gr <- data.frame(
      chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
      start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
      end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
      stringsAsFactors = FALSE)  %>%
      makeGRangesFromDataFrame

   # map to closest gene
   region.closest <- get_region_target_gene(
      regions.gr = regions.gr,
      genome = "hg19",
      method = "genes.promoter.overlap"
   )

   # map to all gene within region +- 250kbp
   region.window <- get_region_target_gene(
      regions.gr = regions.gr,
      genome = "hg19",
      method = "window"
   )

   # map to all gene within region +- 250kbp
   region.nearby.genes <- get_region_target_gene(
      regions.gr = regions.gr[2:3],
      genome = "hg38",
      method = "nearby.genes",
      num.flanking.genes = 5
   )

   expect_true("target" %in% colnames(region.window))
   expect_true("target" %in% colnames(region.closest))
   expect_true("target" %in% colnames(region.nearby.genes))
   expect_true("regionID" %in% colnames(region.window))
   expect_true("regionID" %in% colnames(region.closest))
   expect_true("regionID" %in% colnames(region.nearby.genes))

})

test_that("get_distance_region_target labels correctly for genes in + strand", {

   # PRDM16
   # ENSG00000142611
   # chr: 1
   # Start 3,069,168
   # End 3,438,621

   region.target.after <- data.frame(
      regionID = c( "chr1:3120352-3120353"),
      target = c( "ENSG00000142611"), # PRDM16
      stringsAsFactors = FALSE
   )

   ###################################################
   #                     -->
   #  region.before      |    region.after
   ###################################################

   result <- get_distance_region_target(
      region.target = region.target.after,
      genome = "hg38"
   )
   expect_true("distance_region_target_tss" %in% colnames(result))
   expect_true("target_tss_pos_in_relation_to_region" %in% colnames(result))
   expect_true("region_pos_in_relation_to_gene_tss" %in% colnames(result))
   expect_true(result$distance_region_target_tss > 0)
   expect_true(result$target_tss_pos_in_relation_to_region == "left")
   expect_true(result$region_pos_in_relation_to_gene_tss == "downstream")

   region.target.before <- data.frame(
      regionID = c( "chr1:3020352-3020353"),
      target = c( "ENSG00000142611"), # PRDM16
      stringsAsFactors = FALSE)


   result <- get_distance_region_target(
      region.target = region.target.before,
      genome = "hg38"
   )

   expect_true("distance_region_target_tss" %in% colnames(result))
   expect_true("target_tss_pos_in_relation_to_region" %in% colnames(result))
   expect_true("region_pos_in_relation_to_gene_tss" %in% colnames(result))
   expect_true(result$distance_region_target_tss < 0)
   expect_true(result$target_tss_pos_in_relation_to_region == "right")
   expect_true(result$region_pos_in_relation_to_gene_tss == "upstream")

})



test_that("get_distance_region_target labels correctly for genes in - strand", {
   # CASKP1
   # ENSG00000234059
   # chr: Y
   # Start 12,929,856 (hg38) 12930165
   # End 12,948,185 (hg38)

   ###################################################
   #                   <--
   #  region.before      |    region.after
   ###################################################
   region.target.after <- data.frame(
      regionID = c( "chrY:12830800-12830801"),
      target = c( "ENSG00000234059"), # PRDM16
      stringsAsFactors = FALSE
   )

   result <- get_distance_region_target(
      region.target = region.target.after,
      genome = "hg38"
   )
   expect_true("distance_region_target_tss" %in% colnames(result))
   expect_true("target_tss_pos_in_relation_to_region" %in% colnames(result))
   expect_true("region_pos_in_relation_to_gene_tss" %in% colnames(result))
   expect_true(result$distance_region_target_tss > 0)
   expect_true(result$target_tss_pos_in_relation_to_region == "right")
   expect_true(result$region_pos_in_relation_to_gene_tss == "downstream")

   region.target.before <- data.frame(
      regionID = c( "chrY:13930800-13930801"),
      target = c( "ENSG00000234059"), # PRDM16
      stringsAsFactors = FALSE
   )

   result <- get_distance_region_target(
      region.target = region.target.before,
      genome = "hg38"
   )

   expect_true("distance_region_target_tss" %in% colnames(result))
   expect_true("target_tss_pos_in_relation_to_region" %in% colnames(result))
   expect_true("region_pos_in_relation_to_gene_tss" %in% colnames(result))
   expect_true(result$distance_region_target_tss < 0)
   expect_true(result$target_tss_pos_in_relation_to_region == "left")
   expect_true(result$region_pos_in_relation_to_gene_tss == "upstream")

})

