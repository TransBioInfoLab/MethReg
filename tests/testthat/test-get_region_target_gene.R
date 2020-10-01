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
