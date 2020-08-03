test_that("make_granges_from_names returns a GRanges", {
    library(GenomicRanges)
    gr <- make_granges_from_names("chr22:11-20")
    expect_true(is(gr, "GRanges"))
    expect_equal(width(gr), 10)
})

test_that("make_names_from_granges returns a string", {
    library(GenomicRanges)
    regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
    regions.gr <- make_granges_from_names(regions.names)
    names <- make_names_from_granges(regions.gr)
    expect_equal(names, regions.names)
    expect_true(is(names,"character"))
})

test_that("map_probes_to_regions changes probes cpgs to regions", {
    data(dna.met.chr21)
    dna.met.chr21.with.region.name <- map_probes_to_regions(dna.met.chr21[1:2,])
    expect_true(all(grepl("chr21",rownames(dna.met.chr21.with.region.name))))
    expect_true(is(dna.met.chr21.with.region.name,"matrix"))
})

test_that("map_ensg_to_symbol returns the correct mapping", {
    data(dna.met.chr21)
    symbol <- map_ensg_to_symbol("ENSG00000073282")
    expect_equal("TP63",symbol)
    expect_true(is(symbol,"character"))
})

test_that("map_ensg_to_symbol returns the correct mapping", {
    data(dna.met.chr21)
    ensg <- map_symbol_to_ensg("TP63")
    expect_equal("ENSG00000073282",ensg)
    expect_true(is(ensg,"character"))
})


test_that("get_gene_information returns Granges or dataframes", {
    data(dna.met.chr21)
    gene.info <- get_gene_information(as.granges = TRUE)
    expect_true(is(gene.info,"GRanges"))
    gene.info <- get_gene_information(as.granges = FALSE)
    expect_true(is(gene.info,"data.frame"))
})


test_that("make_se_from_dnam_probes returns a SE with regions", {
    dna.met.chr21 <- get(data("dna.met.chr21"))
    se <- make_se_from_dnam_probes(dna.met.chr21)
    expect_true(is(se,"SummarizedExperiment"))
    expect_true(all(grepl("chr21",rownames(se))))
})

test_that("make_se_from_dnam_regions returns a SE with regions", {
    dna.met.chr21 <- get(data("dna.met.chr21"))
    dna.met.chr21.regions <- map_probes_to_regions(dna.met.chr21)
    se <- make_se_from_dnam_regions(dna.met.chr21.regions)
    expect_true(is(se,"SummarizedExperiment"))
    expect_true(all(grepl("chr21",rownames(se))))
})

test_that("make_se_from_dnam_regions returns a SE with ENSG as rownames", {
    gene.exp.chr21 <- get(data("gene.exp.chr21"))
    se <- make_se_from_gene_matrix(gene.exp.chr21)
    expect_true(is(se,"SummarizedExperiment"))
    expect_true(all(grepl("ENSG",rownames(se))))
    #expect_equal("ENSG00000142156",rowData(se)$ensembl_gene_id[rowData(se)$external_gene_name == "COL6A1"])
})


