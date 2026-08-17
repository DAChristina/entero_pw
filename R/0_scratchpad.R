# test
library(tidyverse)
source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}
source_all()

################################################################################
# Extract tools info
# Klebsiella tool versions 
k_kleborate <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-kleborate.csv"))
k_kaptive<- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-kaptive.csv"))
k_mlst<- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-mlst.csv"))
k_speciator <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-speciator.csv"))
k_core<- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-core.csv"))
k_lincodes<- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-klebsiella-lincodes.csv"))

# E. coli tool versions 
e_mlst<- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-mlst.csv"))
e_mlst2<- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-mlst2.csv"))
e_speciator <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-speciator.csv"))
e_metrics<- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-metrics.csv"))

# version table 
tool_versions <- tibble::tribble(
  ~Tool, ~Version, ~Species, ~Source_column,
  
  # Klebsiella
  "Kleborate", as.character(k_kleborate$Kleborate.version[1]), "K. pneumoniae","kleb_ori_pathogenwatch-kleborate.csv",
  "Kleborate wrapper", as.character(k_kleborate$Wrapper.version[1]), "K. pneumoniae","kleb_ori_pathogenwatch-kleborate.csv",
  "Kaptive", as.character(k_kaptive$Kaptive.Version[1]), "K. pneumoniae","kleb_ori_pathogenwatch-kaptive.csv",
  "MLST (Klebsiella)", as.character(k_mlst$Version[1]), "K. pneumoniae","kleb_ori_pathogenwatch-mlst.csv",
  "Speciator (Klebsiella)",as.character(k_speciator$Version[1]), "K. pneumoniae","kleb_ori_pathogenwatch-speciator.csv",
  "Core genome (WHO_K)", as.character(k_core$Version[1]), "K. pneumoniae","kleb_ori_pathogenwatch-core.csv",
  
  # E. coli
  "MLST Achtman (E. coli)", as.character(e_mlst$Version[1]), "E. coli", "ecoli_ori_pathogenwatch-mlst.csv",
  "MLST Pasteur (E. coli)", as.character(e_mlst2$Version[1]), "E. coli", "ecoli_ori_pathogenwatch-mlst2.csv",
  "Speciator (E. coli)", as.character(e_speciator$Version[1]), "E. coli", "ecoli_ori_pathogenwatch-speciator.csv",
) %>%
  dplyr::mutate(Version = if_else(is.na(Version) | Version == "NULL",
                                  "Not recorded in output",
                                  Version))

write.csv(tool_versions,
          "inputs/pathogenwatch_tools.csv",
          row.names = FALSE)
