context("test-extractBoundaries")

test_that("Whether extractBoundaries gives us the same output", {

    data("arrowhead_gm12878_5kb")

    bounds.GR <- extractBoundaries(domains.mat=arrowhead_gm12878_5kb,
                                   filter=FALSE,
                                   CHR=paste0("CHR",c(1:8,10:22)),
                                   resolution=5000)

    expect_equal(length(bounds.GR), 15468)

})
