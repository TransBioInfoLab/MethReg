context("Map motif probes to regions")
library(GenomicRanges)
genome <- "hg38"
arrayType <- "450k"
motifs.probes <- get(
    data(list = paste0("Probes.motif.",genome,".",toupper(arrayType)),
         package = "ELMER.data",
         envir = environment())
)

test_that("Mapping motifs only to regions overlapping single probes", {

    regions.gr <- get_met_probes_info(genome = genome,arrayType = arrayType)
    regions.gr <- regions.gr[!regions.gr$MASK_general]
    regions.gr <- regions.gr[1:3]

    motif.mapped <- map_motif_probes_to_regions(motifs.probes,
                                                genome,
                                                arrayType,
                                                regions.gr)
    expect_equal(nrow(motif.mapped),3)
    expect_true(all(make_names_from_granges(regions.gr) %in% rownames(motif.mapped)))
    expect_equal(motifs.probes[names(regions.gr)[1],colnames(motif.mapped)]  %>% as.logical(),
                 motif.mapped[make_names_from_granges(regions.gr[1]),]  %>% as.logical())
    expect_equal(motifs.probes[names(regions.gr)[2],colnames(motif.mapped)] %>% as.logical(),
                 motif.mapped[make_names_from_granges(regions.gr[2]),] %>% as.logical())
    expect_equal(motifs.probes[names(regions.gr)[3],colnames(motif.mapped)] %>% as.logical(),
                 motif.mapped[make_names_from_granges(regions.gr[3]),] %>% as.logical())

})

test_that("Mapping motifs only to regions overlapping multiple probes", {

    regions.gr <- get_met_probes_info(genome = genome,arrayType = arrayType)
    regions.gr <- regions.gr[!regions.gr$MASK_general]
    regions.gr <- reduce(regions.gr[5:10] + 500,ignore.strand =T)
    # region 1: chr1 898303-899477
    probes.region1 <- c("cg13938959","cg12445832","cg23999112")
    # region 2: chr1 901656-903607
    probes.region2 <- c("cg11527153","cg27573606","cg04195702")
    motif.mapped <- map_motif_probes_to_regions(motifs.probes,
                                                genome,
                                                arrayType,
                                                regions.gr)
    expect_equal(nrow(motif.mapped),2)
    expect_true(all(make_names_from_granges(regions.gr) %in% rownames(motif.mapped)))
    expect_equal(motif.mapped[1,"AHR_HUMAN.H11MO.0.B"],
                 sum(motifs.probes[probes.region1,"AHR_HUMAN.H11MO.0.B"])
    )
    expect_equal(motif.mapped[2,"GSC2_HUMAN.H11MO.0.D"],
                 sum(motifs.probes[probes.region2,"GSC2_HUMAN.H11MO.0.D"])
    )
    expect_equal(motif.mapped[2,"SP1_HUMAN.H11MO.0.A"],
                 sum(motifs.probes[probes.region2,"SP1_HUMAN.H11MO.0.A"])
    )
})


test_that("Mapping motifs only to regions overlapping single/multiple probes", {

    regions.gr <- get_met_probes_info(genome = genome,arrayType = arrayType)
    regions.gr <- regions.gr[!regions.gr$MASK_general]
    regions.gr <- c(reduce(regions.gr[5:10] + 500,ignore.strand = TRUE),regions.gr[1:2])
    # region 1: chr1 898303-899477
    probes.region1 <- c("cg13938959","cg12445832","cg23999112")
    # region 2: chr1 901656-903607
    probes.region2 <- c("cg11527153","cg27573606","cg04195702")
    motif.mapped <- map_motif_probes_to_regions(motifs.probes,
                                                genome,
                                                arrayType,
                                                regions.gr)
    expect_equal(nrow(motif.mapped),4)
    expect_true(all(make_names_from_granges(regions.gr) %in% rownames(motif.mapped)))
    expect_equal(motif.mapped["chr1:69591-69592","AHR_HUMAN.H11MO.0.B"],
                 sum(motifs.probes["cg21870274","AHR_HUMAN.H11MO.0.B"])
    )

    expect_equal(motif.mapped["chr1:864703-864704","GSC2_HUMAN.H11MO.0.D"],
                 sum(motifs.probes["cg08258224","GSC2_HUMAN.H11MO.0.D"])
    )
    expect_equal(motif.mapped["chr1:898303-899477","SP1_HUMAN.H11MO.0.A"],
                 sum(motifs.probes[probes.region1,"SP1_HUMAN.H11MO.0.A"])
    )

    expect_equal(motif.mapped["chr1:901656-903607","SP1_HUMAN.H11MO.0.A"],
                 sum(motifs.probes[probes.region2,"SP1_HUMAN.H11MO.0.A"])
    )
})


test_that("No mapping throws an error", {
    regions.gr <- get_met_probes_info(genome = genome,arrayType = arrayType)
    regions.gr <- regions.gr[!regions.gr$MASK_general]

    expect_error(
        map_motif_probes_to_regions(motifs.probes,
                                    genome,
                                    arrayType,
                                    shift(regions.gr[1],10))
    )

})


