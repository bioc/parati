parati 0.99.6 (2026-03-16)
-------------------------
* Revised the package interface to better align with Bioconductor standards.
* Updated `parati_run()`:
  - now accepts either a VCF file path or a `VariantAnnotation::VCF` object
  - now returns R objects by default instead of writing files to disk
  - reduced repeated merging inside loops by collecting intermediate results first
* Added internal helpers for:
  - reading family tables from file paths, `data.frame`, or `data.table`
  - converting `VariantAnnotation::VCF` objects into internal `data.table` representation
* Improved Bioconductor integration:
  - added support for `VariantAnnotation`, `SummarizedExperiment`,
    `BiocGenerics`, and `GenomeInfoDb`
  - updated vignette to demonstrate integration with Bioconductor VCF workflows
* Removed incomplete placeholder functionality:
  - removed `vcf_to_plink()`
* Improved documentation:
  - updated function documentation and return value sections
  - cleaned up roxygen2-generated man pages
* Improved tests:
  - replaced minimal structural tests with unit tests covering expected
    transmission inference behavior on toy data
  - ensured core tests are suitable for Bioconductor build checking
* Updated vignette:
  - added abstract and motivation for inclusion in Bioconductor
  - replaced static code blocks with executable R code chunks
  - provided runnable examples using the included toy data
* Minor code cleanup:
  - removed redundant `requireNamespace()` calls for imported packages
  - updated package metadata and `biocViews`
