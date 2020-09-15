test_that("get_residuals works", {
    data("gene.exp.chr21.log2")
    data("clinical")
    metadata <- clinical[,c( "gender", "sample_type")]
    cnv <- matrix(
        sample(x = c(-2,-1,0,1,2),
               size = ncol(gene.exp.chr21.log2) * nrow(gene.exp.chr21.log2),
               replace = TRUE
        ),
        nrow = nrow(gene.exp.chr21.log2),
        ncol = ncol(gene.exp.chr21.log2)
    )
    rownames(cnv) <- rownames(gene.exp.chr21.log2)
    colnames(cnv) <- colnames(gene.exp.chr21.log2)
    gene.exp.residuals <- get_residuals(
        data = gene.exp.chr21.log2[1:3,],
        metadata.samples = metadata,
        metadata.genes = cnv
    )
    gene.exp.residuals <- get_residuals(
        data = gene.exp.chr21.log2[1:3,],
        metadata.samples = metadata,
        metadata.genes = cnv[1:2,]
    )
    gene.exp.residuals <- get_residuals(
        data = gene.exp.chr21.log2[1:3,],
        metadata.samples = metadata
    )
})
