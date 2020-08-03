test_that("get_residuals works", {
    data("gene.exp.chr21")
    data("clinical")
    metadata <- clinical[,c( "gender", "sample_type")]
    cnv <- matrix(
        sample(x = c(-2,-1,0,1,2),
               size = ncol(gene.exp.chr21) * nrow(gene.exp.chr21),
               replace = TRUE
        ),
        nrow = nrow(gene.exp.chr21),
        ncol = ncol(gene.exp.chr21)
    )
    rownames(cnv) <- rownames(gene.exp.chr21)
    colnames(cnv) <- colnames(gene.exp.chr21)
    gene.exp.residuals <- get_residuals(
        data.matrix = gene.exp.chr21[1:3,],
        metadata.samples = metadata,
        metadata.genes = cnv
    )
    gene.exp.residuals <- get_residuals(
        data.matrix = gene.exp.chr21[1:3,],
        metadata.samples = metadata,
        metadata.genes = cnv[1:2,]
    )
    gene.exp.residuals <- get_residuals(
        data.matrix = gene.exp.chr21[1:3,],
        metadata.samples = metadata
    )
})
