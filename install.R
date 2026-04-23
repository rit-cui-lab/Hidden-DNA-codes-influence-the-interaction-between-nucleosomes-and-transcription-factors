# install.R
# Installs R package dependencies for the dNPS / domain analysis pipeline.
# Run once before using any R scripts in this repository:
#   Rscript install.R

cran_packages <- c("dplyr", "ggplot2", "tidyr", "stringr", "tibble")

to_install <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}
