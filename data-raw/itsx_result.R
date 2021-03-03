## code to prepare `DATASET` dataset goes here
# load sample data from inferrnal
seqfile <- system.file(
    file.path("extdata", "sample.fasta"),
    package = "inferrnal"
)
# ITSx has trouble with some of the reads
seq <- Biostrings::readDNAStringSet(seqfile)[c(1,4,20,32,33,43,46,49)]
itsx_result <- rITSx::itsx(
    in_file = seq,
    positions = TRUE,
    complement = FALSE,
    cpu = 1,
    read_function = Biostrings::readDNAStringSet
)
usethis::use_data(itsx_result, overwrite = TRUE)
