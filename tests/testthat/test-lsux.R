futile.logger::flog.threshold(futile.logger::WARN, "LSUx")

test_file <- inferrnal::sample_rRNA_fasta()
test_that(
    "no regression in lsux result",
    {

        skip_on_cran()
        lsux_file <- lsux(test_file)
        expect_known_value(lsux_file, file = "lsu")
    }
)

test_sread <- ShortRead::readFasta(test_file)[1]
test_file <- tempfile(pattern = "sample_rRNA", fileext = "fasta")
ShortRead::writeFasta(test_sread, test_file)
test_DNAss <- Biostrings::readDNAStringSet(test_file)
test_RNAss <- Biostrings::RNAStringSet(test_DNAss)
test_RNAchar <- as.character(test_RNAss)
test_DNAchar <- as.character(test_DNAss)

test_that(
    "lsux gives same results for file, character, ShortRead, DNAStringSet, and RNAStringSet",
    {
        lsux_file <- lsux(test_file)
        lsux_sread <- lsux(test_sread)
        lsux_DNAss <- lsux(test_DNAss)
        lsux_RNAss <- lsux(test_RNAss)
        lsux_DNAchar <- lsux(test_DNAchar)
        lsux_RNAchar <- lsux(test_RNAchar)

        expect_equal(lsux_file, lsux_sread)
        expect_equal(lsux_file, lsux_DNAss)
        expect_equal(lsux_file, lsux_RNAss)
        expect_equal(lsux_file, lsux_DNAchar)
        expect_equal(lsux_file, lsux_RNAchar)
    }
)

test_that(
    "empty result for query with no 5.8S",
    {
        expect_equal(
            lsux(c("test" = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")),
            tibble::tibble(
                seq_id = character(),
                length = integer(),
                region = character(),
                start = integer(),
                end = integer()
            )
        )
    }
)
