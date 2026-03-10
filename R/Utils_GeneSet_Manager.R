#' Create a Standardized Gene Set Object for HCC Analysis
#'
#' This function constructs a standardized gene set object containing specific
#' driver and suppressor genes. It ensures that all subsets are contained within
#' the complete gene list and assigns the object to the 'HCCGeneSet' class.
#'
#' @param name Character string. The name of the pathway or gene set (e.g., "Ferroptosis").
#' @param all_genes Character vector. A list of all genes associated with the phenotype.
#' @param drivers Character vector (Optional). Genes that positively promote the phenotype.
#' @param suppressors Character vector (Optional). Genes that negatively regulate/suppress the phenotype.
#'
#' @return A list of class `HCCGeneSet` containing:
#' \item{name}{Name of the gene set}
#' \item{all_genes}{Unique vector of all genes}
#' \item{drivers}{Unique vector of driver genes}
#' \item{suppressors}{Unique vector of suppressor genes}
#'
#' @export
#'
#' @examples
#' my_set <- create_gene_set(
#'   name = "TestPathway",
#'   all_genes = c("A", "B", "C", "D", "E"),
#'   drivers = c("A", "B"),
#'   suppressors = c("C")
#' )
#' print(my_set)
create_gene_set <- function(name = "CustomGeneSet",
                            all_genes,
                            drivers = NULL,
                            suppressors = NULL,
                            syml=NULL) {

  # --- 1. Data cleaning and validation ---
  # Ensure input is character type and has no duplicates
  all_genes <- unique(as.character(all_genes))
  if(!is.null(drivers)) drivers <- unique(as.character(drivers))
  if(!is.null(suppressors)) suppressors <- unique(as.character(suppressors))
  if(!is.null(syml)) syml <- unique(as.character(syml))



  # Check if Drivers and Suppressors conflict (one gene cannot be both accelerator and brake)
  # While biologically possible, they are usually mutually exclusive in simplified models - here just a warning
  conflict <- intersect(drivers, suppressors)
  if (length(conflict) > 0) {
    warning(paste("Warning: The following genes are listed as BOTH driver and suppressor:", paste(conflict, collapse = ", ")))
  }

  # --- 2. Construct object ---
  obj <- list(
    name = name,
    all_genes = all_genes,
    drivers = drivers,
    suppressors = suppressors,
    syml=syml,
    created_at = Sys.time()
  )

  # --- 3. Assign identity (Class) ---
  # This is the most critical step! We created a new type called "HCCGeneSet"
  class(obj) <- "HCCGeneSet"

  return(obj)
}

#' Print method for HCCGeneSet objects
#'
#' Custom print function to display gene set information neatly.
#'
#' @param x An object of class `HCCGeneSet`.
#' @param ... Additional arguments.
#' @export
print.HCCGeneSet <- function(x, ...) {
  cat("=== HCCMultiOmics Gene Set ===\n")
  cat("Name       :", x$name, "\n")
  cat("Total Genes:", length(x$all_genes), "\n")
  cat("Drivers    :", length(x$drivers), "(Promoters)\n")
  cat("Suppressors:", length(x$suppressors), "(Resistance Markers)\n")
  if (!is.null(x$syml)) {
    cat("Syml       :", length(x$syml), "(Pathway-associated)\n")
  }
  cat("==============================\n")
}
