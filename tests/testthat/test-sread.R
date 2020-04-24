test_file <- inferrnal::sample_rRNA_fasta()
test_sread <- ShortRead::readFasta(test_file)
test_DNAss <- Biostrings::readDNAStringSet(test_file)

sread_file <- methods::as(test_file, "ShortRead")
sread_sread <- methods::as(test_sread, "ShortRead")
sread_DNAss <- methods::as(test_DNAss, "ShortRead")
test_that("loading fasta by Biostrings and ShortRead yield same result from as_sread", {
  expect_equal(sread_file, sread_sread)
  expect_equal(sread_file, sread_DNAss)
  expect_equal(sread_sread, sread_DNAss)
})
