ss <- "..(((.....[[...<<<<<___>>>>>.[[.<<<___>>>.]"

test_that("string secondary structure repair works", {
        expect_identical(repair_unmatched_secondary_structure(ss),
                         "..:::.....::...<<<<<___>>>>>.-[.<<<___>>>.]")
    })

test_in <- system.file(file.path("extdata", "fungi_32S.stk"), package = "LSUx")
test_out <- tempfile(pattern = "trunc", fileext = ".stk")

test_that("alignment truncation works", {
    expect_null(truncate_alignment(test_in, test_out, 1, 500))
    expect_null(
        inferrnal::cmbuild("/dev/null", msafile = test_out, force = TRUE)
    )
})

test_that(
    "alignment truncation works with connections",
    {
        alncon <- file(test_out)
        expect_null(truncate_alignment(test_in, alncon, 1, 500))
        test_aln <- file(test_in)
        alncon <- textConnection("test_trunc", open = "w")
        expect_null(truncate_alignment(test_aln, alncon, 1, 500))
        expect_equal(readLines(test_out), test_trunc)
    }
)

test_that(
    "truncate_alignment rejects bad connections",
    {
        expect_error(truncate_alignment(3, test_out, 1, 500))
        expect_error(truncate_alignment(test_in, NA, 1, 500))
    }
)
