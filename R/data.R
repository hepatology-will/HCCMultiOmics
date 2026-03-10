#' Internal Datasets for HCCMultiOmics Package
#'
#' These datasets are used internally by the package functions for
#' demonstration and analysis purposes.
#'
#' @name HCCMultiOmics-internal
#' @keywords datasets
NULL

#' TCGA Liver Cancer Clinical Data
#'
#' Clinical information for TCGA-LIHC samples including survival data.
#'
#' @format A data frame with columns: SampleID, Group, OS, OS.time, ...
#' @source TCGA Research Network: https://www.cancer.gov/tcga
"TCGA_Clin_Data"

#' TCGA Liver Cancer Expression Matrix
#'
#' Gene expression matrix for TCGA-LIHC samples (normalized).
#'
#' @format A matrix with genes as rows and samples as columns.
#' @source TCGA Research Network: https://www.cancer.gov/tcga
"TCGA_Expr_Mat"

#' TCGA-LIHC Differentially Expressed Genes
#'
#' List of differentially expressed genes identified in TCGA-LIHC.
#'
#' @format A character vector of gene symbols.
#' @source TCGA Research Network: https://www.cancer.gov/tcga
"TCGA_LIHC_DEGs"

#' HCC GeneCards Gene List
#'
#' Gene list related to hepatocellular carcinoma from GeneCards.
#'
#' @format A character vector of gene symbols.
#' @source GeneCards: https://www.genecards.org
"HCC_GeneCards_List"

#' CopyKAT Prediction Data
#'
#' CopyKAT analysis results for single-cell classification.
#'
#' @format Data frame with cell IDs and copykat predictions.
#' @source CopyKAT analysis
"copykat_data"

#' Ferroptosis Genes from FerrDb
#'
#' Ferroptosis-related genes including drivers, suppressors and markers.
#'
#' @format A HCCGeneSet object or data frame.
#' @source FerrDb: http://www.nmdon.cn/fed/
"Ferroptosis_FerrDb"

#' Ferroptosis Drivers
#'
#' Genes that promote ferroptosis.
#'
#' @format A character vector of gene symbols.
#' @source FerrDb: http://www.nmdon.cn/fed/
"ADCD"

#' Ferroptosis Suppressors
#'
#' Genes that suppress ferroptosis.
#'
#' @format A character vector of gene symbols.
#' @source FerrDb: http://www.nmdon.cn/fed/
"Autosis"

#' Cuproptosis Genes from FerrDb
#'
#' Cuproptosis-related genes from FerrDb database.
#'
#' @format A data frame or list.
#' @source FerrDb: http://www.nmdon.cn/fed/
"Cuproptosis_FerrDb"

#' Disulfidptosis Genes
#'
#' Genes associated with disulfidptosis cell death.
#'
#' @format A data frame or list.
"Disulfidptosis"

#' Immunogenic Cell Death Genes
#'
#' Genes associated with immunogenic cell death.
#'
#' @format A data frame or list.
"immunogonic_cell_death"

#' Mitotic Death Genes
#'
#' Genes associated with mitotic death.
#'
#' @format A data frame or list.
"mitotic_death"

#' Parthanatos Genes
#'
#' Genes associated with parthanatos cell death.
#'
#' @format A data frame or list.
"parthanatos"
