test_that("get_promoter_avg throws error if not promoter regions/probes", {
    data("dna.met.chr21")
    skip_on_bioc()
    expect_error(
        promoter.avg <- get_promoter_avg(
            dnam = dna.met.chr21[c("cg02858594",
                                   "cg02858594"),],
            genome = "hg19",
            arrayType = "450k"
        ),
        regexp = "No overlap found between promoter regions and DNA methylation array found")
})

test_that("get_promoter_avg works correctly merges regions overlapping a promoter region", {
    skip_on_bioc()

    # promoter.regions <- get_promoter_regions("hg38")
    # chr1 203727581-203731580      * | ENSG00000221643            SNORA77
    dnam <- t(matrix(c(rep(0,10), rep(1,10)),ncol = 2))
    rownames(dnam) <- c("chr1:203727581-203728580","chr1:203722580-203731580")
    promoter.avg <- get_promoter_avg(
        dnam = dnam,
        genome = "hg38",
        arrayType = "450k"
    )
    expect_true(all(assay(promoter.avg) == 0.5))
})

test_that("get_promoter_avg works correctly merges regions overlapping a promoter region", {
    skip_on_bioc()

    # promoter.regions <- get_promoter_regions("hg38")
    # chr1 203727581-203731580      * | ENSG00000221643            SNORA77
    dnam <- t(matrix(c(rep(0,10), rep(1,10)),ncol = 2))
    rownames(dnam) <- c("chr1:203727581-203728580","chr2:203722580-203731580")
    colnames(dnam) <- paste0("Samples",1:10)
    promoter.avg <- get_promoter_avg(
        dnam = dnam,
        genome = "hg38",
        arrayType = "450k"
    )
    expect_true(all(assay(promoter.avg) == 0.0))
})


test_that("get_promoter_avg works correctly merges regions overlapping a promoter region", {
    skip_on_bioc()

    # promoter.regions <- get_promoter_regions("hg38")
    # chr1 203727581-203731580      * | ENSG00000221643            SNORA77
    dnam <- t(matrix(c(rep(0,10), rep(1,10)),ncol = 2))
    rownames(dnam) <- c("chr3:203727581-203728580","chr1:203722580-203731580")
    colnames(dnam) <- paste0("Samples",1:10)
    promoter.avg <- get_promoter_avg(
        dnam = dnam,
        genome = "hg38",
        arrayType = "450k"
    )
    expect_true(all(assay(promoter.avg) == 1.0))
})

