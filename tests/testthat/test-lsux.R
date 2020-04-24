test_file <- inferrnal::sample_rRNA_fasta()
test_sread <- ShortRead::readFasta(test_file)
test_DNAss <- Biostrings::readDNAStringSet(test_file)
test_RNAss <- Biostrings::RNAStringSet(test_DNAss)
test_RNAchar <- as.character(test_RNAss)
test_DNAchar <- as.character(test_DNAss)

lsux_file <- lsux(test_file)
lsux_sread <- lsux(test_sread)
lsux_DNAss <- lsux(test_DNAss)
lsux_RNAss <- lsux(test_RNAss)
lsux_DNAchar <- lsux(test_DNAchar)
lsux_RNAchar <- lsux(test_RNAchar)

test_that(
    "lsux gives same results for file, character, ShortRead, DNAStringSet, and RNAStringSet",
    {
        expect_equal(lsux_file, lsux_sread)
        expect_equal(lsux_file, lsux_DNAss)
        expect_equal(lsux_file, lsux_RNAss)
        expect_equal(lsux_file, lsux_DNAchar)
        expect_equal(lsux_file, lsux_RNAchar)
    }
)

test_that(
    "no regression in lsux result",
    expect_known_value(lsux_file, file = "lsu")
)
