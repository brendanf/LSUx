test_file <- inferrnal::sample_rRNA_fasta()
test_sread <- ShortRead::readFasta(test_file)
test_DNAss <- Biostrings::readDNAStringSet(test_file)

sread_file <- as_sread(test_file)
sread_sread <- as_sread(test_sread)
sread_DNAss <- as_sread(test_DNAss)
test_that("loading fasta by Biostrings and ShortRead yield same result from as_sread", {
  expect_equal(sread_file, sread_sread)
  expect_equal(sread_file, sread_DNAss)
  expect_equal(sread_sread, sread_DNAss)
})
