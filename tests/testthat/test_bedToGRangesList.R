context("test-bedToGRangesList")

test_that("Whether bedToGRangesList gives us the same output", {

    path = system.file("extdata", package = "preciseTAD")

    tfbsList <- bedToGRangesList(filepath=path, pattern = "*.bed", signal=4)

    expect_equal(length(tfbsList[[1]]), 1841)

    expect_equal(length(tfbsList[[2]]), 2077)

})
