test_that("get_tf_in_region is able to scan", {

    skip_if_offline()
    regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
    region.tf <- get_tf_in_region(
        region = regions.names,
        genome = "hg38"
    )
    # human.tfs <- get_human_tfs()
    # expect_true(all(region.tf$TF %in% human.tfs$ensembl_gene_id))
    expect_type(region.tf$TF_external_gene_name,"character")
    expect_type(region.tf$regionID,"character")
    expect_type(region.tf$TF,"character")
})

test_that("get_tf_in_region accepts granges as input", {
    regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
    region.tf <- get_tf_in_region(
        region = regions.names %>% make_granges_from_names(),
        genome = "hg19"
    )
})


test_that("Wrong genome throws an error", {
    regions.names <- c("chr3:189631389-189632889","chr4:43162098-43163498")
    expect_error(
        get_tf_in_region(
            region = regions.names %>% make_granges_from_names(),
            genome = "xxx"
        )
    )
})


