futile.logger::flog.threshold(futile.logger::WARN, "LSUx")
cm <- system.file(file.path("extdata", "fungi_32S_LR5.cm"), package = "LSUx")

test_file <- inferrnal::sample_rRNA_fasta()
test_that(
    "no regression in lsux result",
    {

        skip_on_cran()
        lsux_file <- lsux(test_file, cm_32S = cm)
        expect_known_value(lsux_file, file = "lsu")
    }
)

test_sread <- ShortRead::readFasta(test_file)[1:3]
test_quality <- lapply(Biostrings::width(test_sread), sample.int, n = 40, replace = TRUE)
test_quality <- lapply(test_quality, Biostrings::PhredQuality)
test_quality <- do.call(c, test_quality)
test_sreadq <- ShortRead::ShortReadQ(test_sread@sread, test_quality, test_sread@id)
test_file <- tempfile(pattern = "sample_rRNA", fileext = ".fasta")
ShortRead::writeFasta(test_sread, test_file)
test_fastq <- tempfile(pattern = "sample_rRNA", fileext = ".fastq")
ShortRead::writeFastq(test_sreadq, test_fastq)
test_DNAss <- Biostrings::readDNAStringSet(test_file)
test_RNAss <- Biostrings::RNAStringSet(test_DNAss)
test_RNAchar <- as.character(test_RNAss)
test_DNAchar <- as.character(test_DNAss)
lsux_file <- lsux(test_file, cm_32S = cm)

test_that(
    "lsux gives same results for fasta, fastq, character, ShortRead, DNAStringSet, and RNAStringSet",
    {
        lsux_sread <- lsux(test_sread, cm_32S = cm)
        lsux_DNAss <- lsux(test_DNAss, cm_32S = cm)
        lsux_RNAss <- lsux(test_RNAss, cm_32S = cm)
        lsux_DNAchar <- lsux(test_DNAchar, cm_32S = cm)
        lsux_RNAchar <- lsux(test_RNAchar, cm_32S = cm)
        lsux_fastq <- lsux(test_fastq, cm_32S = cm)

        expect_equal(lsux_file, lsux_sread)
        expect_equal(lsux_file, lsux_DNAss)
        expect_equal(lsux_file, lsux_RNAss)
        expect_equal(lsux_file, lsux_DNAchar)
        expect_equal(lsux_file, lsux_RNAchar)
        expect_equal(lsux_file, lsux_fastq)
    }
)

test_that(
    "lsux gives same results when query is unnamed",
    {
        names(test_DNAss) <- NULL
        lsux_noname <- lsux(test_DNAss)
        expect_true(
            dplyr::all_equal(
                spread_regions(lsux_noname)[-1],
                spread_regions(lsux_file)[-1]
            )
        )
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

test_that(
    "progressive memory allocation works",
    {
        expect_error(lsux(test_sread, cm_32S = cm, mxsize = 1, cpu = 1))
        expect_equal(
            lsux(test_sread, cm_32S = cm, mxsize = c(1, 1024), cpu = 1),
            lsux_file
        )
    }
)
