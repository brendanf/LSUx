ss <- "..(((.....[[...<<<<<___>>>>>.[[.<<<___>>>.]"

test_that(
    "string secondary structure repair works", {
        expect_identical(repair_unmatched_secondary_structure(ss),
                         "..:::.....::...<<<<<___>>>>>.-[.<<<___>>>.]")
    })

test_out <- tempfile(pattern = "trunc", fileext = ".stk")

test_that(
    "alignment truncation works",
    {
        expect_null(
            truncate_alignment(
                system.file(file.path("extdata", "fungi_32S.stk"),
                            package = "LSUx"),
                test_out,
                1,
                500
            )
        )
        expect_null(
            inferrnal::cmbuild("/dev/null", msafile = test_out, force = TRUE)
        )
    })
