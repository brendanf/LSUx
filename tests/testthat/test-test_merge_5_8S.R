seqfile <- system.file(
    file.path("extdata", "sample.fasta"),
    package = "inferrnal"
)
seq <- Biostrings::readDNAStringSet(seqfile)[c(1,4,20,32,33,43,46,49)]
cm_5_8S <- system.file(
    file.path("extdata", "RF00002.cm"),
    package = "inferrnal"
)
cm_result <- inferrnal::cmsearch(seq, cm = cm_5_8S, cpu = 1)
test_that("merge_5_8S works", {
  expect_known_value(merge_5_8S(itsx_result$positions, cm_result), "merge.rds")
})
