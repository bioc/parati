# Internal helper: read family table
.parati_read_family <- function(fam) {
  if (data.table::is.data.table(fam)) {
    fam_dt <- data.table::copy(fam)
  } else if (is.data.frame(fam)) {
    fam_dt <- data.table::as.data.table(fam)
  } else if (is.character(fam) && length(fam) == 1L) {
    fam_dt <- data.table::as.data.table(openxlsx::read.xlsx(fam))
  } else {
    stop("`fam` must be a file path, data.frame, or data.table.")
  }

  required_cols <- c("FamilyIndex", "IndividualID", "Role")
  if (!all(required_cols %in% names(fam_dt))) {
    stop(
      "Family table must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  fam_dt
}

# Internal helper: convert a VariantAnnotation::VCF object to data.table
.parati_vcf_to_dt <- function(vcf_obj) {
  rr <- SummarizedExperiment::rowRanges(vcf_obj)
  gt <- VariantAnnotation::geno(vcf_obj)$GT

  if (is.null(gt)) {
    stop("The VCF object does not contain GT genotype data.")
  }

  fixed_df <- as.data.frame(VariantAnnotation::fixed(vcf_obj))
  chrom <- as.character(GenomeInfoDb::seqnames(rr))
  pos <- BiocGenerics::start(rr)
  id <- names(rr)
  id[is.na(id)] <- "."

  ref <- as.character(unname(VariantAnnotation::ref(vcf_obj)))
  alt <- vapply(
    VariantAnnotation::alt(vcf_obj),
    function(x) paste(as.character(x), collapse = ","),
    character(1)
  )

  qual <- if ("QUAL" %in% names(fixed_df)) {
    fixed_df$QUAL
  } else {
    rep(".", length(chrom))
  }

  filt <- if ("FILTER" %in% names(fixed_df)) {
    vapply(
      fixed_df$FILTER,
      function(x) paste(as.character(x), collapse = ";"),
      character(1)
    )
  } else {
    rep(".", length(chrom))
  }

  vcf_dt <- data.table::data.table(
    `#CHROM` = chrom,
    POS = pos,
    ID = id,
    REF = ref,
    ALT = alt,
    QUAL = qual,
    FILTER = filt,
    INFO = ".",
    FORMAT = "GT"
  )

  gt_dt <- data.table::as.data.table(gt)
  cbind(vcf_dt, gt_dt)
}

# Internal helper: read family table
.parati_read_family <- function(fam) {
  if (data.table::is.data.table(fam)) {
    fam_dt <- data.table::copy(fam)
  } else if (is.data.frame(fam)) {
    fam_dt <- data.table::as.data.table(fam)
  } else if (is.character(fam) && length(fam) == 1L) {
    fam_dt <- data.table::as.data.table(openxlsx::read.xlsx(fam))
  } else {
    stop("`fam` must be a file path, data.frame, or data.table.")
  }

  required_cols <- c("FamilyIndex", "IndividualID", "Role")
  if (!all(required_cols %in% names(fam_dt))) {
    stop(
      "Family table must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  fam_dt
}

# Internal helper: convert a VariantAnnotation::VCF object to data.table
.parati_vcf_to_dt <- function(vcf_obj) {
  rr <- SummarizedExperiment::rowRanges(vcf_obj)
  gt <- VariantAnnotation::geno(vcf_obj)$GT

  if (is.null(gt)) {
    stop("The VCF object does not contain GT genotype data.")
  }

  fixed_df <- as.data.frame(VariantAnnotation::fixed(vcf_obj))
  chrom <- as.character(GenomeInfoDb::seqnames(rr))
  pos <- BiocGenerics::start(rr)
  id <- names(rr)
  id[is.na(id)] <- "."

  ref <- as.character(unname(VariantAnnotation::ref(vcf_obj)))
  alt <- vapply(
    VariantAnnotation::alt(vcf_obj),
    function(x) paste(as.character(x), collapse = ","),
    character(1)
  )

  qual <- if ("QUAL" %in% names(fixed_df)) fixed_df$QUAL else rep(".", length(chrom))
  filt <- if ("FILTER" %in% names(fixed_df)) {
    vapply(fixed_df$FILTER, function(x) paste(as.character(x), collapse = ";"), character(1))
  } else {
    rep(".", length(chrom))
  }

  vcf_dt <- data.table::data.table(
    `#CHROM` = chrom,
    POS = pos,
    ID = id,
    REF = ref,
    ALT = alt,
    QUAL = qual,
    FILTER = filt,
    INFO = ".",
    FORMAT = "GT"
  )

  gt_dt <- data.table::as.data.table(gt)
  cbind(vcf_dt, gt_dt)
}

# Internal helper: read family table
.parati_read_family <- function(fam) {
  if (data.table::is.data.table(fam)) {
    fam_dt <- data.table::copy(fam)
  } else if (is.data.frame(fam)) {
    fam_dt <- data.table::as.data.table(fam)
  } else if (is.character(fam) && length(fam) == 1L) {
    fam_dt <- data.table::as.data.table(openxlsx::read.xlsx(fam))
  } else {
    stop("`fam` must be a file path, data.frame, or data.table.")
  }

  required_cols <- c("FamilyIndex", "IndividualID", "Role")
  if (!all(required_cols %in% names(fam_dt))) {
    stop(
      "Family table must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  fam_dt
}

# Internal helper: convert a VariantAnnotation::VCF object to data.table
.parati_vcf_to_dt <- function(vcf_obj) {
  rr <- SummarizedExperiment::rowRanges(vcf_obj)
  gt <- VariantAnnotation::geno(vcf_obj)$GT

  if (is.null(gt)) {
    stop("The VCF object does not contain GT genotype data.")
  }

  fixed_df <- as.data.frame(VariantAnnotation::fixed(vcf_obj))
  chrom <- as.character(GenomeInfoDb::seqnames(rr))
  pos <- BiocGenerics::start(rr)
  id <- names(rr)
  id[is.na(id)] <- "."

  ref <- as.character(unname(VariantAnnotation::ref(vcf_obj)))
  alt <- vapply(
    VariantAnnotation::alt(vcf_obj),
    function(x) paste(as.character(x), collapse = ","),
    character(1)
  )

  qual <- if ("QUAL" %in% names(fixed_df)) fixed_df$QUAL else rep(".", length(chrom))
  filt <- if ("FILTER" %in% names(fixed_df)) {
    vapply(fixed_df$FILTER, function(x) paste(as.character(x), collapse = ";"), character(1))
  } else {
    rep(".", length(chrom))
  }

  vcf_dt <- data.table::data.table(
    `#CHROM` = chrom,
    POS = pos,
    ID = id,
    REF = ref,
    ALT = alt,
    QUAL = qual,
    FILTER = filt,
    INFO = ".",
    FORMAT = "GT"
  )

  gt_dt <- data.table::as.data.table(gt)
  cbind(vcf_dt, gt_dt)
}

# Internal helper: read VCF path or VCF object
.parati_read_vcf <- function(vcf, chr = NULL) {
  if (inherits(vcf, "VCF")) {
    vcf_dt <- .parati_vcf_to_dt(vcf)
  } else if (is.character(vcf) && length(vcf) == 1L) {
    vcf_dt <- data.table::fread(
      file = vcf,
      skip = "#CHROM",
      sep = "\t",
      header = TRUE,
      data.table = TRUE,
      fill = TRUE
    )
  } else {
    stop("`vcf` must be a file path or a VariantAnnotation::VCF object.")
  }

  if (!"#CHROM" %in% names(vcf_dt) && "CHROM" %in% names(vcf_dt)) {
    data.table::setnames(vcf_dt, "CHROM", "#CHROM")
  }

  required_vcf_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  if (!all(required_vcf_cols %in% names(vcf_dt))) {
    stop(
      "VCF input is missing required columns: ",
      paste(setdiff(required_vcf_cols, names(vcf_dt)), collapse = ", ")
    )
  }

  if (!is.null(chr)) {
    vcf_dt <- vcf_dt[vcf_dt[["#CHROM"]] == as.character(chr), ]
  }

  vcf_dt
}
