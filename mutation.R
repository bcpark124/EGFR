library(dplyr)
library(stringr)
library(scales)

clinical_sample <- fread(
  r"(/data_clinical_sample.txt)",
  sep = "\t",
  header = TRUE,
  skip = 4
)

maf <- fread(
  r"(/data_mutations_extended.txt)",
  sep = "\t",
  header = TRUE,
)

nsclc_samples <- clinical_sample %>%
  filter(CANCER_TYPE == "Non-Small Cell Lung Cancer")

total_nsclc_patients <- nsclc_samples %>%
  pull(PATIENT_ID) %>%
  unique() %>%
  length()

sample_map <- nsclc_samples %>% select(PATIENT_ID, SAMPLE_ID)

egfr_maf_all <- maf %>%
  filter(Hugo_Symbol == "EGFR") %>%
  inner_join(sample_map, by = c("Tumor_Sample_Barcode" = "SAMPLE_ID"))

total_egfr_all_mutations <- nrow(egfr_maf_all)

egfr_nsclc_patients <- egfr_maf_all %>%
  pull(PATIENT_ID) %>%
  unique() %>%
  length()

cds_variants <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Silent"
)

egfr_maf <- egfr_maf_all %>%
  filter(Variant_Classification %in% cds_variants)

egfr_cds_nsclc_patients <- egfr_maf %>%
  pull(PATIENT_ID) %>%
  unique() %>%
  length()

egfr_maf_cleaned <- egfr_maf %>%
  mutate(
    clean_label = str_remove(HGVSp_Short, "^p\\.")
  ) %>%
  mutate(
    all_nums = str_extract_all(clean_label, "\\d+")
  ) %>%
  filter(sapply(all_nums, length) > 0) %>%
  mutate(
    start_pos = as.numeric(sapply(all_nums, function(x) x[1])),
    end_pos   = as.numeric(sapply(all_nums, function(x) if(length(x) > 1) x[2] else x[1])),
    end_pos   = if_else(grepl("del", clean_label) & end_pos < start_pos, start_pos, end_pos)
  ) %>%
  mutate(
    is_any_ex19del =
      Variant_Classification == "In_Frame_Del" &
      start_pos >= 729 & start_pos <= 761,
    
    is_classical_ex19del =
      is_any_ex19del &
      start_pos >= 745 & start_pos <= 754 &
      end_pos >= 745 & end_pos <= 754,
    
    is_atypical_ex19del =
      is_any_ex19del & !is_classical_ex19del,
    
    is_ex20ins =
      Variant_Classification == "In_Frame_Ins" &
      start_pos >= 762 & start_pos <= 823
  ) %>%
  select(-all_nums) %>%
  distinct(PATIENT_ID, clean_label, .keep_all = TRUE)

total_unique_patient_mutations <- nrow(egfr_maf_cleaned)

n_range_mutations_with_multiples <- egfr_maf_cleaned %>%
  filter(start_pos >= 688 & start_pos <= 982) %>%
  nrow()

pct_range_mutations_with_multiples <- if_else(
  total_mutations_with_multiples > 0, 
  n_range_mutations_with_multiples / total_mutations_with_multiples, 
  0
)

n_total_ex19del     <- sum(egfr_maf_cleaned$is_any_ex19del)
n_classical_ex19del <- sum(egfr_maf_cleaned$is_classical_ex19del)
n_atypical_ex19del  <- sum(egfr_maf_cleaned$is_atypical_ex19del)

ex19del_subtypes <- egfr_maf_cleaned %>%
  filter(n_classical_ex19del) %>%
  group_by(clean_label) %>%
  summarise(
    Count = n(),
    Frequency = Count / total_unique_patient_mutations,
    representative_pos = min(start_pos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

ex19del_top20_print <- ex19del_subtypes %>%
  head(20) %>%
  mutate(Rank = row_number(), Frequency_Pct = scales::percent(Frequency, accuracy = 0.1)) %>%
  select(Rank, clean_label, Count, Frequency_Pct)

print(as.data.frame(ex19del_top20_print), row.names = FALSE)

ex20ins_subtypes <- egfr_maf_cleaned %>%
  filter(is_ex20ins) %>%
  group_by(clean_label) %>%
  summarise(
    Count = n(),
    Frequency = Count / total_unique_patient_mutations,
    representative_pos = min(start_pos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

ex20ins_top20_print <- ex20ins_subtypes %>%
  head(20) %>%
  mutate(Rank = row_number(), Frequency_Pct = scales::percent(Frequency, accuracy = 0.1)) %>%
  select(Rank, clean_label, Count, Frequency_Pct)

print(as.data.frame(ex20ins_top20_print), row.names = FALSE)
ex19del_summary <- egfr_maf_cleaned %>%
  filter(n_classical_ex19del) %>%
  summarise(
    AA_pos_adjusted = 746, 
    mut_label = "Ex19del",
    display_label = "Ex19del",
    n_mutations = n(),
    .groups = "drop"
  )

top_1_ex19 <- head(ex19del_subtypes, 1)

ex19del_top_variant <- tibble(
  AA_pos_adjusted = if_else(nrow(top_1_ex19) > 0, top_1_ex19$representative_pos, 746),
  mut_label = "Ex19del",
  display_label = if_else(nrow(top_1_ex19) > 0, top_1_ex19$clean_label, NA_character_),
  n_mutations = if_else(nrow(top_1_ex19) > 0, top_1_ex19$Count, 0L)
)

ex20ins_summary <- egfr_maf_cleaned %>%
  filter(is_ex20ins) %>%
  summarise(
    AA_pos_adjusted = 770, 
    mut_label = "Ex20ins",
    display_label = "Ex20ins",
    n_mutations = n(),
    .groups = "drop"
  )

top_1_ex20 <- head(ex20ins_subtypes, 1)

ex20ins_top_variant <- tibble(
  AA_pos_adjusted = if_else(nrow(top_1_ex20) > 0, top_1_ex20$representative_pos, 770),
  mut_label = "Ex20ins",
  display_label = if_else(nrow(top_1_ex20) > 0, top_1_ex20$clean_label, NA_character_),
  n_mutations = if_else(nrow(top_1_ex20) > 0, top_1_ex20$Count, 0L)
)

target_top_labels <- c(ex19del_top_variant$display_label, ex20ins_top_variant$display_label)

pure_other_summary <- egfr_maf_cleaned %>%
  filter(!n_classical_ex19del & !is_ex20ins & !(clean_label %in% target_top_labels)) %>%
  group_by(clean_label) %>%
  summarise(
    AA_pos_adjusted = min(start_pos, na.rm = TRUE),
    mut_label = "Other",
    display_label = first(clean_label),
    n_mutations = n(),
    .groups = "drop"
  )
