.onLoad <- function(libname, pkgname) {
  invisible()
}

.check_and_install <- function(pkg, source = c("CRAN", "Bioc")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(TRUE)
  }

  source <- match.arg(source)

  msg <- sprintf("Package '%s' is not installed. Installing...", pkg)
  message(msg)

  if (source == "Bioc") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE)
    }
    BiocManager::install(pkg, force = TRUE, quiet = TRUE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
  }

  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Package '%s' installed successfully.", pkg))
    return(TRUE)
  } else {
    stop(sprintf("Failed to install package '%s'. Please install manually.", pkg))
    return(FALSE)
  }
}

.check_bioc_packages <- function(pkgs) {
  for (pkg in pkgs) {
    .check_and_install(pkg, source = "Bioc")
  }
}

.check_cran_packages <- function(pkgs) {
  for (pkg in pkgs) {
    .check_and_install(pkg, source = "CRAN")
  }
}
