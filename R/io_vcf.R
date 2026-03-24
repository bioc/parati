#' Read VCF by chromosome
#'
#' Reads a VCF file and subsets variants to a given chromosome.
#'
#' @param vcf_file Character scalar, path to a VCF/VCF.GZ file.
#' @param chr Character or integer chromosome identifier.
#'
#' @return A `data.table` containing VCF rows for the selected chromosome.
#'
#' @examples
#' vcf_file <- system.file("extdata", "Toy_TrioGenotype.vcf.gz", package = "parati")
#' vcf_chr <- read_vcf_by_chr(vcf_file, chr = 1)
#' dim(vcf_chr)
#'
#' @export
read_vcf_by_chr <- function(vcf_file, chr) {
  vcf_all <- data.table::fread(
    file = vcf_file,
    skip = "#CHROM",
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    fill = TRUE
  )

  if (!"#CHROM" %in% names(vcf_all) && "CHROM" %in% names(vcf_all)) {
    data.table::setnames(vcf_all, "CHROM", "#CHROM")
  }

  vcf_chr <- vcf_all[vcf_all[["#CHROM"]] == as.character(chr), ]
  vcf_chr
}

#' Write VCF data.table to file
#'
#' Writes a `data.table` representing VCF rows to a VCF file.
#'
#' @param df A `data.table` containing VCF rows.
#' @param file Character scalar, output file path.
#'
#' @return `NULL`, invisibly. The VCF file is written to `file`.
#'
#' @examples
#' vcf_file <- system.file("extdata", "Toy_TrioGenotype.vcf.gz", package = "parati")
#' vcf_dt <- read_vcf_by_chr(vcf_file, chr = 1)
#' outfile <- tempfile(fileext = ".vcf")
#' write_vcf_dt(vcf_dt, outfile)
#' file.exists(outfile)
#'
#' @export
write_vcf_dt <- function(df, file) {
  if (grepl("\\.gz$", file)) {
    data.table::fwrite(df, file = gzfile(file), sep = "\t", quote = FALSE)
  } else {
    data.table::fwrite(df, file = file, sep = "\t", quote = FALSE)
  }
  invisible(NULL)
}

#' Convert data.table to vcfR object
#'
#' Converts a `data.table` representing VCF rows into a `vcfR` object.
#'
#' @param df A `data.table` containing VCF rows.
#' @param meta Character vector of VCF meta lines.
#'
#' @return A `vcfR` object.
#'
#' @examples
#' vcf_file <- system.file("extdata", "Toy_TrioGenotype.vcf.gz", package = "parati")
#' vcf_dt <- read_vcf_by_chr(vcf_file, chr = 1)
#' vcf_obj <- vcf_dt_to_vcfR(vcf_dt)
#' class(vcf_obj)
#'
#' @export
vcf_dt_to_vcfR <- function(df, meta = character()) {
  fix_mat <- as.matrix(df[, seq_len(8)])
  gt_mat <- as.matrix(df[, 9:ncol(df), with = FALSE])

  cls <- methods::getClass("vcfR", where = asNamespace("vcfR"))
  methods::new(cls, meta = meta, fix = fix_mat, gt = gt_mat)
}

#' Write vcfR object to file
#'
#' Writes a `vcfR` object to a VCF file.
#'
#' @param vcf_obj A `vcfR` object.
#' @param file Character scalar, output file path.
#'
#' @return `NULL`, invisibly. The VCF file is written to `file`.
#'
#' @examples
#' vcf_file <- system.file("extdata", "Toy_TrioGenotype.vcf.gz", package = "parati")
#' vcf_dt <- read_vcf_by_chr(vcf_file, chr = 1)
#' vcf_obj <- vcf_dt_to_vcfR(vcf_dt)
#' outfile <- tempfile(fileext = ".vcf")
#' write_vcf_obj(vcf_obj, outfile)
#' file.exists(outfile)
#'
#' @export
write_vcf_obj <- function(vcf_obj, file) {
  vcfR::write.vcf(vcf_obj, file = file)
  invisible(NULL)
}
