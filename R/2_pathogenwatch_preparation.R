# Compile results from ABRicate
library(tidyverse)
source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}
source_all()

# Batch list:
# Batch 1 (xx August 2026)
# Batch 2 

################################################################################
# Load all _ori_ files
# Species specific coz' different infromation & tools among species

# 1. Eschericia coli ###########################################################
e_path <- "raw_data/pw_result_test/ecoli_ori/"

e_summary <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-summary.csv"))
e_metrics <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-metrics.csv"))
e_speciator <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-speciator.csv"))
e_mlst <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-mlst.csv"))
e_mlst2 <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-mlst2.csv"))
e_clermont <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-clermontyping.csv"))
e_ectyper <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-ectyper.csv"))
e_hclink <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-hclink.csv"))
e_amr <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-rfplus.csv"))
e_vir <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-virulencefinder.csv"))
e_stec <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-stecfinder.csv"))
e_plasmid <- read.csv(paste0(e_path, "ecoli_ori_pathogenwatch-plasmidfinder.csv"))

# qc
e_interp_qc <- e_summary %>%
  janitor::clean_names() %>% 
  dplyr::select(genome_id, genome_name, qc,
                species_prediction #,
                # insdc_run_accession, country, date
  ) %>%
  dplyr::left_join(
    e_metrics %>%
      janitor::clean_names() %>% 
      dplyr::select(genome_id, genome_length, no_contigs,
                    n50, gc_content, ns_per_100_kbp,
                    largest_contig)
    ,
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    e_speciator %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, species_name, mash_distance)
    ,
    by = "genome_id"
  ) %>%
  glimpse()

write.csv(e_interp_qc,
          paste0(e_path, "ecoli_interp_qc_assembly.csv"),
          row.names = FALSE)

# Typing
e_interp_typing <- e_mlst %>%
  janitor::clean_names() %>% 
  dplyr::select(genome_id, genome_name,
                st_achtman = st,
                adk, fum_c, gyr_b, icd, mdh, pur_a, rec_a
  ) %>%
  dplyr::left_join(
    e_mlst2 %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, st_pasteur = st,
                    din_b, icd_a, pab_b, pol_b, put_p, trp_a, trp_b, uid_a
      )
    ,
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    e_clermont %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, phylogroup)
    ,
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    e_ectyper %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, o_type, h_type, serotype)
    ,
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    e_stec %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, stx_type)
    ,
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    e_hclink %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, closest_st, hc10, hc100, hc200,
                    hc400, hc1100, hc1500, hc2000, hc2350)
    ,
    by = "genome_id"
  ) %>%
  glimpse()

write.csv(e_interp_typing,
          paste0(e_path, "ecoli_interp_typing.csv"),
          row.names = FALSE)

# AMR (binarised)
e_interp_amr <- e_amr %>%
  janitor::clean_names() %>% 
  dplyr::mutate(across(-c(genome_id, genome_name), binarise)) %>%
  glimpse()

write.csv(e_interp_amr,
          paste0(e_path, "ecoli_interp_amr.csv"),
          row.names = FALSE)

# Virulence (binarised, complete/partial
e_interp_vir <- e_vir %>%
  janitor::clean_names() %>% 
  dplyr::mutate(across(-c(genome_id, genome_name),
                       ~ if_else(is.na(.) | str_detect(., "(?i)absent"), 0L, 1L))) %>%
  # flag partial AMR
  dplyr::bind_cols(
    e_vir %>%
      janitor::clean_names() %>% 
      dplyr::select(-genome_id, -genome_name) %>%
      dplyr::rename_with(~ paste0(., "_partial")) %>%
      dplyr::mutate(across(everything(),
                           ~ if_else(!is.na(.) & str_detect(., "(?i)partial"), 1L, 0L)))
  ) %>%
  glimpse()

write.csv(e_interp_vir,
          paste0(e_path, "ecoli_interp_virulence.csv"),
          row.names = FALSE)

# Plasmids (wide)
e_interp_plasmid <- e_plasmid %>%
  janitor::clean_names() %>% 
  dplyr::filter(!is.na(inc_match)) %>%
  dplyr::mutate(
    inc_group = str_remove(inc_match, "_\\d+$"),
    high_conf = as.numeric(percent_identity) >= 90 &
      as.numeric(match_coverage) >= 80
  ) %>%
  dplyr::filter(high_conf) %>%
  dplyr::distinct(genome_id, genome_name, inc_group) %>%
  dplyr::mutate(present = 1L) %>%
  tidyr::pivot_wider(names_from= inc_group,
                     values_from = present,
                     values_fill = 0L) %>%
  janitor::clean_names() %>% 
  glimpse()

write.csv(e_interp_plasmid,
          paste0(e_path, "ecoli_interp_plasmids.csv"),
          row.names = FALSE)

# merged
e_interp_master <- e_interp_qc %>%
  dplyr::left_join(e_interp_typing,by = "genome_id") %>%
  dplyr::left_join(e_interp_amr,by = "genome_id") %>%
  dplyr::left_join(e_interp_vir,by = "genome_id") %>%
  dplyr::left_join(e_interp_plasmid, by = "genome_id")

write.csv(e_interp_master,
          paste0(e_path, "ecoli_interp_master.csv"),
          row.names = FALSE)


# 1. Klebsiella pneumoniae #####################################################
k_path <- "raw_data/pw_result_test/kleb_ori/"

k_summary <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-summary.csv"))
k_metrics <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-metrics.csv"))
k_speciator <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-speciator.csv"))
k_mlst <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-mlst.csv"))
k_kaptive <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-kaptive.csv"))
k_kleborate <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-kleborate.csv"))
k_lincodes <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-klebsiella-lincodes.csv"))
k_plasmid <- read.csv(paste0(k_path, "kleb_ori_pathogenwatch-plasmidfinder.csv"))

# qc
k_interp_qc <- k_summary %>%
  janitor::clean_names() %>% 
  dplyr::select(genome_id, genome_name, qc, species_prediction,
                insdc_run_accession, country, date) %>%
  dplyr::left_join(
    k_metrics %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, genome_length, no_contigs,
                    n50, gc_content, ns_per_100_kbp),
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    k_speciator %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, species_name, mash_distance),
    by = "genome_id"
  ) %>%
  glimpse()

write.csv(k_interp_qc,
          paste0(k_path, "kleb_interp_qc_assembly.csv"),
          row.names = FALSE)

# Typing
k_interp_typing <- k_mlst %>%
  janitor::clean_names() %>% 
  dplyr::select(genome_id, st, gap_a, inf_b, mdh, pgi, pho_e, rpo_b, ton_b) %>%
  dplyr::left_join(
    k_kaptive %>% 
      janitor::clean_names() %>% 
      dplyr::select(
        genome_id,
        K_locus      = k_locus_best_match_locus,
        K_type       = k_locus_best_match_type,
        K_confidence = k_locus_match_confidence,
        K_identity   = k_locus_identity,
        K_problems   = k_locus_problems,
        K_missing_genes = k_locus_missing_expected_genes,
        
        O_locus      = o_locus_best_match_locus,
        O_type       = o_locus_best_match_type,
        O_confidence = o_locus_match_confidence,
        O_identity   = o_locus_identity,
        O_problems   = o_locus_problems
      ),
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    k_lincodes %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, cg_st, lin_code,
                    sublineage, clonal_group),
    by = "genome_id"
  ) %>%
  dplyr::left_join(
    k_kleborate %>% 
      janitor::clean_names() %>% 
      dplyr::select(genome_id, clonal_complex,
                    virulence_score, resistance_score,
                    yersiniabactin, colibactin, aerobactin,
                    salmochelin, rmp_adc, rmp_a2),
    by = "genome_id"
  ) %>%
  glimpse()

write.csv(k_interp_typing,
          paste0(k_path, "kleb_interp_typing.csv"),
          row.names = FALSE)

# AMR (kleborate)
amr_cols_kleb <- c("a_gly_acquired",
                   "col_acquired",
                   
                   "fcyn_acquired",
                   "flq_acquired",
                   "mls_acquired",
                   "phe_acquired",
                   "rif_acquired",
                   "sul_acquired",
                   "tet_acquired",
                   "tmt_acquired",
                   
                   "bla_acquired",
                   "bla_esbl_acquired",
                   "bla_carb_acquired",
                   "shv_mutations",
                   "flq_mutations",
                   "omp_mutations")

k_interp_amr <- k_kleborate %>%
  janitor::clean_names() %>%
  dplyr::select(genome_id, genome_name, resistance_score,
                resistance_class_count, resistance_gene_count,
                all_of(amr_cols_kleb)) %>%
  dplyr::mutate(across(all_of(amr_cols_kleb), binarise))

write.csv(k_interp_amr,
          paste0(k_path, "kleb_interp_amr.csv"),
          row.names = FALSE)

# Plasmids
k_interp_plasmid <- k_plasmid %>%
  janitor::clean_names() %>% 
  dplyr::filter(!is.na(inc_match)) %>%
  dplyr::mutate(inc_group = str_remove(inc_match, "_\\d+$")) %>%
  dplyr::filter(as.numeric(percent_identity) >= 90,
                as.numeric(match_coverage)>= 80) %>%
  dplyr::distinct(genome_id, genome_name, inc_group) %>%
  dplyr::mutate(present = 1L) %>%
  tidyr::pivot_wider(names_from= inc_group,
                     values_from = present,
                     values_fill = 0L) %>%
  janitor::clean_names() %>% 
  glimpse()

write.csv(k_interp_plasmid,
          paste0(k_path, "kleb_interp_plasmids.csv"),
          row.names = FALSE)

# Combine!
k_interp_master <- k_interp_qc %>%
  dplyr::left_join(k_interp_typing,by = "genome_id") %>%
  dplyr::left_join(k_interp_amr,by = "genome_id") %>%
  dplyr::left_join(k_interp_plasmid, by = "genome_id")

write.csv(k_interp_master,
          paste0(k_path, "kleb_interp_master.csv"),
          row.names = FALSE)
