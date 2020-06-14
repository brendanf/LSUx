
# LSUx

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/brendanf/LSUx.svg?branch=master)](https://travis-ci.com/brendanf/LSUx)
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

    ## INFO [2020-06-14 22:00:26] Beginning CM search.
    ## INFO [2020-06-14 22:00:27] 48 sequences contained a single 5.8S hit.
    ## INFO [2020-06-14 22:00:27] Beginning CM alignment.
    ## INFO [2020-06-14 22:00:45] Extracting LSU regions.

``` r
regions
```

    ## # A tibble: 480 x 5
    ##    seq_id length region start   end
    ##    <chr>   <int> <chr>  <int> <int>
    ##  1 seq1     1404 ITS1       1   114
    ##  2 seq1     1404 5_8S     115   271
    ##  3 seq1     1404 ITS2     272   412
    ##  4 seq1     1404 LSU1     413   517
    ##  5 seq1     1404 V2       518   673
    ##  6 seq1     1404 LSU2     674   832
    ##  7 seq1     1404 V3       833  1106
    ##  8 seq1     1404 LSU3    1107  1147
    ##  9 seq1     1404 V4      1148  1242
    ## 10 seq1     1404 LSU4    1243  1404
    ## # … with 470 more rows
