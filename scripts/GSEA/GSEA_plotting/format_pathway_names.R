#' Format Pathway Names with Smart Capitalization
#'
#' Formats pathway names with proper capitalization while preserving common
#' biological abbreviations, roman numerals, and chemical nomenclature.
#'
#' @param text Character vector of pathway names to format
#' @param use_formatting Logical, if TRUE applies smart formatting, if FALSE returns cleaned text only (default: TRUE)
#' @param strip_prefix Logical, if TRUE removes common database prefixes (default: TRUE)
#'
#' @return Character vector of formatted pathway names
#'
#' @examples
#' format_pathway_name("HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING")
#' # Returns: "TNFR1-Induced NF-kappaB Signaling"
#'
#' format_pathway_name("GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY")
#' # Returns: "Type II Interferon Signaling Pathway"
#'
#' format_pathway_name("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING")
#' # Returns: "Class I MHC-Mediated Antigen Processing"
#'
#' @export
format_pathway_name <- function(text, use_formatting = TRUE, strip_prefix = TRUE) {

  if (length(text) == 0 || all(is.na(text))) {
    return(text)
  }

  # Store original for fallback
  original <- text

  # Step 1: Replace underscores with spaces and collapse multiple spaces
  text <- stringr::str_replace_all(text, "_", " ")
  text <- stringr::str_replace_all(text, " +", " ")  # Collapse multiple spaces

  # Step 2: Strip common prefixes if requested
  if (strip_prefix) {
    common_prefixes <- c(
      "HALLMARK ", "KEGG ", "REACTOME ", "BIOCARTA ", "MEDICUS ",
      "GOBP ", "GOCC ", "GOMF ", "PID ", "WIKIPATHWAY ",
      "WP ", "^GO ", "GTRD ", "NABA ", "HP "
    )
    for (prefix in common_prefixes) {
      text <- stringr::str_replace(text, paste0("^", prefix), "")
    }
  }

  # If formatting disabled, return cleaned text in title case
  if (!use_formatting) {
    return(stringr::str_to_title(text))
  }

  # Step 3: Apply smart formatting with exceptions
  text <- apply_smart_capitalization(text)

  return(text)
}

#' Apply Smart Capitalization with Biological Exceptions
#'
#' Internal function that applies intelligent capitalization rules
#' preserving biological abbreviations, roman numerals, and chemical terms.
#'
#' @param text Character vector to format
#' @return Formatted character vector
#' @keywords internal
apply_smart_capitalization <- function(text) {

  # Build exception dictionary
  exceptions <- build_exception_dictionary()

  # Build multi-word pattern dictionary (order matters: longest first)
  multiword_patterns <- build_multiword_patterns()

  # Process each string
  sapply(text, function(s) {
    # Convert to lowercase first for consistent matching
    s_lower <- tolower(s)

    # First pass: Replace multi-word patterns (marked with << >>)
    for (pattern_info in multiword_patterns) {
      s_lower <- gsub(pattern_info$pattern, pattern_info$replacement, s_lower, fixed = TRUE)
    }

    # Split into words (will split around << >> markers too)
    words <- unlist(strsplit(s_lower, " "))

    # Process each word
    formatted_words <- sapply(seq_along(words), function(i) {
      word <- words[i]

      # Skip empty words
      if (nchar(word) == 0) return("")

      # Check if this is a protected multi-word term (starts with <<)
      if (grepl("^<<", word)) {
        # Extract the protected content and return as-is
        protected <- gsub("<<|>>", "", word)
        return(protected)
      }

      # Check if word matches an exception (case-insensitive)
      exception_match <- exceptions[tolower(names(exceptions)) == word]

      if (length(exception_match) > 0) {
        # Use exception formatting
        return(exception_match[[1]])
      }

      # Check for roman numerals (I, II, III, IV, V, VI, VII, VIII, IX, X)
      if (grepl("^(i|ii|iii|iv|v|vi|vii|viii|ix|x)$", word)) {
        return(toupper(word))
      }

      # Check for numbered items (e.g., "il2", "stat5")
      # Format as IL-2, STAT5, etc.
      if (grepl("^([a-z]+)(\\d+[a-z]*)$", word)) {
        prefix <- sub("^([a-z]+)(\\d+[a-z]*)$", "\\1", word)
        number <- sub("^([a-z]+)(\\d+[a-z]*)$", "\\2", word)

        # Check if prefix is a known abbreviation
        prefix_exception <- exceptions[tolower(names(exceptions)) == prefix]
        if (length(prefix_exception) > 0) {
          return(paste0(prefix_exception[[1]], "-", number))
        }
      }

      # Check for Greek letters (but not kappa in "nf kappa b" context)
      greek_letters <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta",
                        "eta", "theta", "iota", "lambda", "mu", "nu",
                        "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon",
                        "phi", "chi", "psi", "omega")

      if (word %in% greek_letters) {
        return(word)  # Keep lowercase for Greek letters
      }

      # Check for chemical prefixes (n-, o-, p-, cis-, trans-, etc.)
      if (grepl("^(n|o|p|m|s|cis|trans|alpha|beta|gamma|delta)-", word)) {
        # Keep prefix lowercase, capitalize rest
        parts <- unlist(strsplit(word, "-"))
        return(paste(parts[1], stringr::str_to_title(parts[-1]), sep = "-", collapse = "-"))
      }

      # Check for ordinal suffixes (1st, 2nd, 3rd, etc.)
      if (grepl("^\\d+(st|nd|rd|th)$", word)) {
        return(word)
      }

      # Default: Title case (first letter capitalized)
      return(stringr::str_to_title(word))
    })

    # Reassemble
    result <- paste(formatted_words, collapse = " ")

    # Clean up any remaining markers
    result <- gsub("<<|>>", "", result)

    return(result)
  }, USE.NAMES = FALSE)
}

#' Build Multi-Word Pattern Dictionary
#'
#' Creates patterns for multi-word biological terms that need special handling.
#' These are applied BEFORE word-by-word processing.
#'
#' @return List of pattern/replacement pairs
#' @keywords internal
build_multiword_patterns <- function() {

  patterns <- list(
    # NF-kappaB family (most important - match all variations)
    list(pattern = "nf kappa b", replacement = "<<NF-kappaB>>"),
    list(pattern = "nf-kappa-b", replacement = "<<NF-kappaB>>"),
    list(pattern = "nfkappab", replacement = "<<NF-kappaB>>"),
    list(pattern = "nfkb", replacement = "<<NF-kappaB>>"),

    # TGF-beta
    list(pattern = "tgf beta", replacement = "<<TGF-beta>>"),
    list(pattern = "tgfbeta", replacement = "<<TGF-beta>>"),

    # IFN variants
    list(pattern = "ifn alpha", replacement = "<<IFN-alpha>>"),
    list(pattern = "ifnalpha", replacement = "<<IFN-alpha>>"),
    list(pattern = "ifn beta", replacement = "<<IFN-beta>>"),
    list(pattern = "ifnbeta", replacement = "<<IFN-beta>>"),
    list(pattern = "ifn gamma", replacement = "<<IFN-gamma>>"),
    list(pattern = "ifngamma", replacement = "<<IFN-gamma>>"),

    # Interferon spelled out
    list(pattern = "interferon alpha", replacement = "<<Interferon-alpha>>"),
    list(pattern = "interferon beta", replacement = "<<Interferon-beta>>"),
    list(pattern = "interferon gamma", replacement = "<<Interferon-gamma>>"),

    # TNF-alpha
    list(pattern = "tnf alpha", replacement = "<<TNF-alpha>>"),
    list(pattern = "tnfa", replacement = "<<TNF-alpha>>"),
    list(pattern = "tumor necrosis factor alpha", replacement = "<<Tumor Necrosis Factor-alpha>>"),

    # HIF-1alpha
    list(pattern = "hif 1 alpha", replacement = "<<HIF-1alpha>>"),
    list(pattern = "hif1alpha", replacement = "<<HIF-1alpha>>"),

    # IL-STAT combinations (common in pathways)
    list(pattern = "il2 stat5", replacement = "<<IL-2/STAT5>>"),
    list(pattern = "il6 jak stat3", replacement = "<<IL-6/JAK/STAT3>>"),

    # Common multi-word process terms
    list(pattern = "cell cycle", replacement = "<<Cell Cycle>>"),
    list(pattern = "cell death", replacement = "<<Cell Death>>"),
    list(pattern = "apoptotic process", replacement = "<<Apoptotic Process>>"),
    list(pattern = "immune response", replacement = "<<Immune Response>>"),
    list(pattern = "inflammatory response", replacement = "<<Inflammatory Response>>"),
    list(pattern = "innate immune", replacement = "<<Innate Immune>>"),
    list(pattern = "adaptive immune", replacement = "<<Adaptive Immune>>"),
    list(pattern = "antigen processing", replacement = "<<Antigen Processing>>"),
    list(pattern = "antigen presentation", replacement = "<<Antigen Presentation>>")
  )

  return(patterns)
}

#' Build Exception Dictionary for Biological Terms
#'
#' Creates a named list of biological abbreviations and their proper formatting.
#' Keys are lowercase for case-insensitive matching.
#'
#' @return Named list with lowercase keys and properly formatted values
#' @keywords internal
build_exception_dictionary <- function() {

  # Common biological abbreviations
  exceptions <- list(
    # NF-kappaB family (multi-word handled in build_multiword_patterns)
    "nf" = "NF",
    "nfkb" = "NF-kappaB",
    "nfkappab" = "NF-kappaB",

    # Interleukins
    "il" = "IL",
    "il1" = "IL-1",
    "il2" = "IL-2",
    "il4" = "IL-4",
    "il6" = "IL-6",
    "il8" = "IL-8",
    "il10" = "IL-10",
    "il12" = "IL-12",
    "il17" = "IL-17",

    # Interferons
    "ifn" = "IFN",
    "ifna" = "IFN-alpha",
    "ifnb" = "IFN-beta",
    "ifng" = "IFN-gamma",
    "interferon" = "Interferon",

    # Tumor necrosis factors
    "tnf" = "TNF",
    "tnfa" = "TNF-alpha",
    "tnfr" = "TNFR",
    "tnfr1" = "TNFR1",
    "tnfr2" = "TNFR2",

    # Immunology
    "mhc" = "MHC",
    "hla" = "HLA",
    "tcr" = "TCR",
    "bcr" = "BCR",
    "tlr" = "TLR",
    "tlr4" = "TLR4",
    "cd" = "CD",
    "cd4" = "CD4",
    "cd8" = "CD8",
    "cd28" = "CD28",
    "cd80" = "CD80",
    "cd86" = "CD86",

    # Signaling molecules
    "stat" = "STAT",
    "stat1" = "STAT1",
    "stat3" = "STAT3",
    "stat5" = "STAT5",
    "jak" = "JAK",
    "mapk" = "MAPK",
    "erk" = "ERK",
    "jnk" = "JNK",
    "pi3k" = "PI3K",
    "akt" = "AKT",
    "mtor" = "mTOR",
    "ampk" = "AMPK",

    # Transcription factors
    "ap" = "AP",
    "ap1" = "AP-1",
    "creb" = "CREB",
    "irf" = "IRF",
    "irf1" = "IRF1",
    "irf3" = "IRF3",
    "irf7" = "IRF7",

    # Cell cycle
    "cdk" = "CDK",
    "apc" = "APC",
    "dna" = "DNA",
    "rna" = "RNA",
    "mrna" = "mRNA",
    "rrna" = "rRNA",
    "trna" = "tRNA",
    "mirna" = "miRNA",

    # Metabolism
    "atp" = "ATP",
    "adp" = "ADP",
    "amp" = "AMP",
    "nadh" = "NADH",
    "nadph" = "NADPH",
    "fadh2" = "FADH2",
    "coa" = "CoA",
    "acetyl" = "Acetyl",
    "tca" = "TCA",
    "oxphos" = "OXPHOS",

    # Lipids
    "ldl" = "LDL",
    "hdl" = "HDL",
    "vldl" = "VLDL",
    "lps" = "LPS",
    "pge2" = "PGE2",

    # Reactive species
    "ros" = "ROS",
    "no" = "NO",
    "nos" = "NOS",

    # Enzymes
    "cox" = "COX",
    "cox2" = "COX-2",
    "lox" = "LOX",
    "lt" = "LT",  # Leukotrienes
    "ex" = "EX",  # Eoxins
    "pla2" = "PLA2",

    # Amino acids (3-letter codes)
    "gln" = "Gln",
    "glu" = "Glu",
    "gls" = "GLS",
    "glul" = "GLUL",
    "asn" = "Asn",
    "asp" = "Asp",
    "arg" = "Arg",

    # Common terms
    "egf" = "EGF",
    "fgf" = "FGF",
    "pdgf" = "PDGF",
    "vegf" = "VEGF",
    "tgf" = "TGF",
    "tgfb" = "TGF-beta",
    "igf" = "IGF",
    "ngf" = "NGF",
    "gdnf" = "GDNF",
    "bdnf" = "BDNF",

    # Receptor types
    "gpcr" = "GPCR",
    "rtk" = "RTK",

    # Organelles
    "er" = "ER",
    "ecm" = "ECM",

    # Misc
    "uv" = "UV",
    "hif" = "HIF",
    "hif1a" = "HIF-1alpha",
    "p53" = "p53",
    "rb" = "Rb",
    "brca1" = "BRCA1",
    "brca2" = "BRCA2",
    "parp" = "PARP",
    "sumo" = "SUMO",
    "ubiquitin" = "Ubiquitin",

    # Processes
    "autophagy" = "Autophagy",
    "apoptosis" = "Apoptosis",
    "necrosis" = "Necrosis",
    "ferroptosis" = "Ferroptosis",
    "pyroptosis" = "Pyroptosis",

    # Direction/location terms
    "via" = "via",
    "and" = "and",
    "or" = "or",
    "of" = "of",
    "in" = "in",
    "to" = "to",
    "by" = "by",
    "from" = "from",
    "the" = "the"
  )

  return(exceptions)
}
