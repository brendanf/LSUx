test_file <- inferrnal::sample_rRNA_fasta()
test_sread <- ShortRead::readFasta(test_file)
test_DNAss <- Biostrings::readDNAStringSet(test_file)

sread_file <- methods::as(test_file, "ShortRead")
sread_sread <- methods::as(test_sread, "ShortRead")
sread_DNAss <- methods::as(test_DNAss, "ShortRead")
test_that(
    "loading fasta by Biostrings and ShortRead yield same result from as_sread",
    {
        expect_equal(sread_file, sread_sread)
        expect_equal(sread_file, sread_DNAss)
        expect_equal(sread_sread, sread_DNAss)
    })

test_that("round-trip converstion between DNAStringSet and ShortRead", {
    expect_equal(test_DNAss, methods::as(sread_DNAss, "DNAStringSet"))
})

test_that("round-trip converstion between character and ShortRead", {
    expect_equal(
        sread_sread,
        methods::as(methods::as(sread_sread, "character"), "ShortRead")
    )
})

test_quality <- lapply(Biostrings::width(test_DNAss), sample.int, n = 40, replace = TRUE)
test_phredquality <- lapply(test_quality, Biostrings::PhredQuality)
test_phredquality <- do.call(c, test_phredquality)
test_qsDNAss <- Biostrings::QualityScaledDNAStringSet(test_DNAss, test_phredquality)
sreadq_qsDNAss <- methods::as(test_qsDNAss, "ShortReadQ")
test_fastq <- tempfile(fileext = ".fastq")
Biostrings::writeQualityScaledXStringSet(test_qsDNAss, test_fastq)
sreadq_sreadq <- ShortRead::readFastq(test_fastq, qualityType = "FastqQuality")
qsDNAss_sreadq <- methods::as(sreadq_sreadq, "QualityScaledDNAStringSet")
test_that(
    "ShortReadQ and QualityScaledDNAStringSet conversion works for base 33", {
    expect_equal(sreadq_qsDNAss@sread, sreadq_sreadq@sread)
    expect_equal(sreadq_qsDNAss@id, sreadq_sreadq@id)
    expect_equal(
        methods::as(sreadq_qsDNAss@quality, "matrix"),
        methods::as(sreadq_sreadq@quality, "matrix")
    )
    expect_equal(
        methods::as(sreadq_qsDNAss@quality@quality, "character"),
        methods::as(sreadq_sreadq@quality@quality, "character")
    )
    expect_equal(
        as.character(test_qsDNAss),
        as.character(qsDNAss_sreadq)
        )
    expect_equal(
        methods::as(test_qsDNAss@quality, "NumericList"),
        methods::as(qsDNAss_sreadq@quality, "NumericList")
    )
    expect_equal(
        methods::as(test_qsDNAss@quality, "character"),
        methods::as(qsDNAss_sreadq@quality, "character")
    )
})


test_solexaquality <- lapply(test_quality, Biostrings::SolexaQuality)
test_solexaquality <- do.call(c, test_solexaquality)
test_qsDNAss <- Biostrings::QualityScaledDNAStringSet(test_DNAss, test_solexaquality)
sreadq_qsDNAss <- methods::as(test_qsDNAss, "ShortReadQ")
test_fastq <- tempfile(fileext = ".fastq")
Biostrings::writeQualityScaledXStringSet(test_qsDNAss, test_fastq)
sreadq_sreadq <- ShortRead::readFastq(test_fastq, qualityType = "SFastqQuality")
qsDNAss_sreadq <- methods::as(sreadq_sreadq, "QualityScaledDNAStringSet")
test_that(
    "ShortReadQ and QualityScaledDNAStringSet conversion works for base 64", {
        expect_equal(sreadq_qsDNAss@sread, sreadq_sreadq@sread)
        expect_equal(sreadq_qsDNAss@id, sreadq_sreadq@id)
        expect_equal(
            methods::as(sreadq_qsDNAss@quality, "matrix"),
            methods::as(sreadq_sreadq@quality, "matrix")
        )
        expect_equal(
            methods::as(sreadq_qsDNAss@quality@quality, "character"),
            methods::as(sreadq_sreadq@quality@quality, "character")
        )
        expect_equal(
            as.character(test_qsDNAss),
            as.character(qsDNAss_sreadq)
        )
        expect_equal(
            methods::as(test_qsDNAss@quality, "NumericList"),
            methods::as(qsDNAss_sreadq@quality, "NumericList")
        )
        expect_equal(
            methods::as(test_qsDNAss@quality, "character"),
            methods::as(qsDNAss_sreadq@quality, "character")
        )
    })

test_that("Conversion from fastq file to ShortReadQ works", {
    expect_equal(sreadq_sreadq, methods::as(test_fastq, "ShortRead"))
})
