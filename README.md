
# LSUx

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/brendanf/LSUx.svg?branch=master)](https://travis-ci.org/brendanf/LSUx)
[![Codecov test
coverage](https://codecov.io/gh/brendanf/LSUx/branch/master/graph/badge.svg)](https://codecov.io/gh/brendanf/LSUx?branch=master)
<!-- badges: end -->

Cut rDNA sequences into domains using covariance models

## Installation

LSUx uses [Infernal](http://eddylab.org/infernal/) (INFERence of RNA
ALignment) to match the 5.8S and LSU regions. This means that Infernal
must be installed for it to function. Follow the instructions on the
Infernal webpage to accomplish this.

Note that Infernal does not work on Windows.

The R interface to Infernal is in package `inferrnal` (two “r”s).

You can install the development versions of `inferrnal` and `LSUx` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("brendanf/inferrnal")
devtools::install_github("brendanf/LSUx")
```

``` r
library(LSUx)
```

## Usage

LSUx extracts alternating variable and conserved domains from the
contiguous rDNA regions which form the eukaryotic ribosomal large
subunit, i.e. 5.8S RNA, ITS2, and 28S RNA. For the purposes of this
document, this region will be referred to as the 32S precursor RNA, as
in humans, although its actual size in Svedberg units varies between
lineages.

Input sequences should contain, at a minimum, a significant fraction of
the 5.8S RNA, which is used to define the 5’ end of 32S. Any base pairs
before the 5’ end of 5.8S will be considered to be ITS1 (`ITS1 = TRUE`)
or discarded (`ITS1 = FALSE`). Input sequences should not extend past
the end of the 32S model at the 3’ end.

LSUx requires two covariance models: one for 5.8S, which is used in to
locate the 5’ end of the target region, and one for 32S. The 5.8S model
can be [RF00002 from Rfam](https://rfam.xfam.org/family/RF00002) (the
default), or an equivalent. If using a custom CM, it must be calibrated
using `cmcalibrate` from Infernal.

The 32S model must include annotations in the reference line (“`#=GC
RF`” in the seed alignment) to distinguish conserved and variable
regions. The annotations should be sequential characters in the range
`1..9A..Z` for conserved domains, “`v`” for variable domains, and “`.`”
for unaligned gaps in the seed alignment. In the output, the conserved
domains will be named “5\_8S”, “LSU1”, “LSU2”, …; the variable domains
will be named “ITS2”, “V1”, “V2”, …

Two example models are included, both based on the [RDP fungal LSU
CM](https://github.com/rdpstaff/fungene_pipeline/blob/master/%%20resources/RRNA_28S/model.cm),
and annotated with variable regions according to Raué (1988). The first,
includes the full LSU region.

``` r
cm_full <- system.file(
    file.path("extdata", "fungal_32S.cm"),
    package = "LSUx"
)
```

The second, is truncated at the binding site of the LR5 primer, and
should be faster for input sequences which do not extend past that
point.

``` r
cm_trunc <- system.file(
    file.path("extdata", "fungal_32S_LR5.cm"),
    package = "LSUx"
)
```

The seed alignments are also provided; they can be accessed using:

``` r
system.file(file.path("extdata", "fungal_32S.stk"), package = "LSUx")
system.file(file.path("extdata", "fungal_32S_LR5.stk"), package = "LSUx")
```

The sample data is from an environmental metabarcoding study focusing on
fungi.

``` r
seq <- system.file("extdata/sample.fasta", package = "inferrnal")
```

The sequences in the sample data were amplified using primers ITS1 and
LR5, so using the truncated CM is appropriate.

``` r
cm_32S_trunc <- system.file(
    file.path("extdata", "fungi_32S_LR5.cm"),
    package = "LSUx"
)
```

LSUx gives a `tibble` with the start and end locations of each region
for each sequence.

``` r
regions <- lsux(seq, cm_32S = cm_32S_trunc, ITS1 = TRUE, cpu = 1)
```

    ## INFO [2020-03-23 18:18:17] Beginning CM search.
    ## INFO [2020-03-23 18:18:20] 48 sequences contained a single 5.8S hit.
    ## INFO [2020-03-23 18:18:24] Beginning CM alignment.
    ## INFO [2020-03-23 18:18:53] Extracting LSU regions.

``` r
regions
```

    ## # A tibble: 48 x 22
    ##    seq_name length ITS1_start ITS1_end `5_8S_start` `5_8S_end` ITS2_start
    ##    <chr>     <int>      <int>    <int>        <int>      <int>      <int>
    ##  1 seq45      1635          1      294          295        451        452
    ##  2 seq3       1571          1      235          236        392        393
    ##  3 seq2       1551          1      192          193        349        350
    ##  4 seq28      1408          1      191          192        348        349
    ##  5 seq23      1447          1      193          194        350        351
    ##  6 seq9       1406          1      169          170        326        327
    ##  7 seq7       1447          1      190          191        347        348
    ##  8 seq48      1535          1      191          192        348        349
    ##  9 seq5       1587          1      256          257        413        414
    ## 10 seq6       1584          1      256          257        413        414
    ## # … with 38 more rows, and 15 more variables: ITS2_end <int>, LSU1_start <int>,
    ## #   LSU1_end <int>, V2_start <int>, V2_end <int>, LSU2_start <int>,
    ## #   LSU2_end <int>, V3_start <int>, V3_end <int>, LSU3_start <int>,
    ## #   LSU3_end <int>, V4_start <int>, V4_end <int>, LSU4_start <int>,
    ## #   LSU4_end <int>
