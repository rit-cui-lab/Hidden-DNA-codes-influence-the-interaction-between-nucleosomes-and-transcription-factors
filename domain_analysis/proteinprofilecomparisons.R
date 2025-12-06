# proteinprofilecomparisons.R
# Author: Chris Carson
# Description: Analyzes protein domains across pioneer vs non-pioneer TFs.
# Compares domain distributions, performs chi-square and Fisher's exact tests.

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# ============================================
# 1. LOAD DATA
# ============================================

# Configuration: Set input files
PIONEER_DOMAINS_FILE <- "Simplified_Domainrecord_PF_pfam.csv"
NONPIONEER_DOMAINS_FILE <- "Simplified_Domainrecord_TF_pfam.csv"
PROTEIN_PROFILES_FILE <- "mammalconsensus.csv"

# Load domain data
pf <- read.csv(PIONEER_DOMAINS_FILE)
nonpf <- read.csv(NONPIONEER_DOMAINS_FILE)

# Combine pioneer and non-pioneer data
all_domains <- bind_rows(
  pf %>% mutate(TF_Category = "Pioneer"),
  nonpf %>% mutate(TF_Category = "NonPioneer")
)

# Load protein-to-profile mapping
protein_profiles <- read.csv(PROTEIN_PROFILES_FILE)

print("Original Protein-Profile Mapping:")
print(head(protein_profiles))
print(paste("Total proteins with profiles:", nrow(protein_profiles)))

# ============================================
# 2. EXPAND PROTEINS WITH MULTIPLE PROFILES
# ============================================

print("\n========================================")
print("EXPANDING MULTI-PROFILE PROTEINS")
print("========================================")

expand_profiles <- function(df) {
  df_expanded <- df %>%
    mutate(
      Profile_List = str_split(Consensus, ",\\s*")
    ) %>%
    unnest(Profile_List) %>%
    rename(Individual_Profile = Profile_List) %>%
    select(-Consensus)
  
  return(df_expanded)
}

protein_profiles_expanded <- expand_profiles(protein_profiles)

print("\nExpanded Protein-Profile Mapping:")
print(head(protein_profiles_expanded, 20))
print(paste("Total protein-profile combinations:", nrow(protein_profiles_expanded)))

protein_counts <- protein_profiles_expanded %>%
  count(Transcription_Factor, name = "Num_Profiles") %>%
  arrange(desc(Num_Profiles))

print("\nProteins with multiple profiles:")
print(protein_counts %>% filter(Num_Profiles > 1))

# ============================================
# 3. JOIN DOMAINS WITH EXPANDED PROFILES
# ============================================

domains_with_profiles <- all_domains %>%
  left_join(protein_profiles_expanded, by = c("Gene" = "Transcription_Factor"), relationship = "many-to-many")

proteins_no_profile <- domains_with_profiles %>%
  filter(is.na(Individual_Profile)) %>%
  distinct(Gene)

if(nrow(proteins_no_profile) > 0) {
  print("\nWarning: The following proteins have domains but no profile assignment:")
  print(proteins_no_profile)
}

domains_with_profiles <- domains_with_profiles %>%
  filter(!is.na(Individual_Profile))

print("\nProfile distribution (after expansion):")
profile_summary <- domains_with_profiles %>%
  distinct(Gene, Individual_Profile) %>%
  count(Individual_Profile) %>%
  arrange(desc(n))
print(profile_summary)

# ============================================
# 4. FILTER TO REGULAR DOMAINS ONLY
# ============================================

domains_regular <- domains_with_profiles %>% 
  filter(Hit_Type != "superfamily")

print(paste("\nTotal regular domain entries:", nrow(domains_regular)))
print(paste("Unique proteins:", n_distinct(domains_regular$Gene)))
print(paste("Unique domains:", n_distinct(domains_regular$Domain)))

print("\nTF_Category distribution:")
print(table(domains_regular$TF_Category))

# ============================================
# 4.5. CONSOLIDATE DOMAIN VARIANTS
# ============================================

print("\n========================================")
print("CONSOLIDATING DOMAIN VARIANTS")
print("========================================")

consolidate_domains <- function(df) {
  df <- df %>%
    mutate(
      Domain_Base = str_replace(Domain, "_\\d+$", ""),
      Domain_Consolidated = Domain_Base
    )
  
  manual_mapping <- c(
    "zf-C4" = "zf-C4_consolidated",
    "ZnF-C4" = "zf-C4_consolidated",
    "zf-H2C2_2" = "zf-C2H2_consolidated",
    "zf-C2H2" = "zf-C2H2_consolidated"
  )
  
  for(from in names(manual_mapping)) {
    df$Domain_Consolidated[df$Domain == from] <- manual_mapping[from]
  }
  
  return(df)
}

domains_regular <- consolidate_domains(domains_regular)

print("\nDomain consolidation applied:")
print("- Domains with trailing _number consolidated")
print("- zf-C4 + ZnF-C4 consolidated")
print("- zf-H2C2_2 + zf-C2H2 consolidated")

print(paste("\nUnique domains after consolidation:", n_distinct(domains_regular$Domain_Consolidated)))

# ============================================
# 5. ANALYZE PIONEER VS NON-PIONEER
# ============================================

print("\n========================================")
print("ANALYZING PIONEER VS NON-PIONEER DOMAINS")
print("========================================")

category_results <- list()
category_plots <- list()

TOP_N_PER_CATEGORY <- 15

for(category in c("Pioneer", "NonPioneer")) {
  
  print(paste("\n--- TF Category:", category, "---"))
  
  category_domains <- domains_regular %>%
    filter(TF_Category == category)
  
  domain_counts <- category_domains %>%
    distinct(Gene, Domain_Consolidated) %>%
    group_by(Domain_Consolidated) %>%
    summarise(Gene_Count = n()) %>%
    arrange(desc(Gene_Count))
  
  print(paste("Total unique domains in", category, ":", nrow(domain_counts)))
  print(paste("Total unique proteins in", category, ":", n_distinct(category_domains$Gene)))
  
  top_domains <- domain_counts %>%
    head(TOP_N_PER_CATEGORY)
  
  print(paste("\nTop", TOP_N_PER_CATEGORY, "domains for", category, ":"))
  print(top_domains)
  
  top_domains$TF_Category <- category
  names(top_domains)[1] <- "Domain"
  category_results[[category]] <- top_domains
  
  if(nrow(top_domains) > 0) {
    fill_color <- if(category == "Pioneer") "#E74C3C" else "#3498DB"
    
    p <- ggplot(top_domains, aes(x = reorder(Domain, Gene_Count), y = Gene_Count)) +
      geom_bar(stat = "identity", fill = fill_color) +
      coord_flip() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      ) +
      labs(title = paste("Top Domains in", category, "TFs"),
           subtitle = paste(n_distinct(category_domains$Gene), "proteins"),
           x = "Domain", 
           y = "Number of Unique Genes")
    
    category_plots[[category]] <- p
    print(p)
  }
}

# ============================================
# 6. COMBINED RESULTS TABLE
# ============================================

all_category_results <- bind_rows(category_results)

print("\n========================================")
print("COMBINED TOP DOMAINS PIONEER vs NON-PIONEER")
print("========================================")
print(all_category_results)

# ============================================
# 7. COMPARE CATEGORY-SPECIFIC DOMAINS
# ============================================

print("\n========================================")
print("CATEGORY-SPECIFIC DOMAIN ANALYSIS")
print("========================================")

for(category in c("Pioneer", "NonPioneer")) {
  category_domains <- domains_regular %>%
    filter(TF_Category == category) %>%
    distinct(Domain_Consolidated) %>%
    pull(Domain_Consolidated)
  
  other_category_domains <- domains_regular %>%
    filter(TF_Category != category) %>%
    distinct(Domain_Consolidated) %>%
    pull(Domain_Consolidated)
  
  unique_domains <- setdiff(category_domains, other_category_domains)
  
  if(length(unique_domains) > 0) {
    print(paste("\nDomains unique to", category, ":", length(unique_domains)))
    print(head(unique_domains, 20))
  } else {
    print(paste("\nNo domains unique to", category))
  }
}

# ============================================
# 8. CHI-SQUARE TEST
# ============================================

print("\n========================================")
print("CHI-SQUARE TEST: PIONEER vs NON-PIONEER")
print("========================================")

domain_counts_wide <- domains_regular %>%
  distinct(Gene, Domain_Consolidated, TF_Category) %>%
  group_by(Domain_Consolidated, TF_Category) %>%
  summarise(Gene_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = TF_Category, 
              values_from = Gene_Count, 
              values_fill = 0)

MIN_GENES_FOR_TEST <- 3
domain_counts_wide <- domain_counts_wide %>%
  mutate(Total = Pioneer + NonPioneer) %>%
  filter(Total >= MIN_GENES_FOR_TEST) %>%
  arrange(Domain_Consolidated)

print(paste("Domains used for chi-square test (Total >=", MIN_GENES_FOR_TEST, "):", nrow(domain_counts_wide)))

contingency_matrix <- as.matrix(domain_counts_wide %>% select(Pioneer, NonPioneer))
rownames(contingency_matrix) <- domain_counts_wide$Domain_Consolidated

domain_counts_wide <- arrange(domain_counts_wide, desc(Total))

write.csv(domain_counts_wide, "contingency_table_chi_square.csv", row.names = FALSE)

chi_test <- chisq.test(contingency_matrix, simulate.p.value = TRUE)

print("\n=== CHI-SQUARE TEST RESULTS ===")
print(chi_test)
print(paste("Chi-square statistic:", round(chi_test$statistic, 4)))
print(paste("P-value:", format(chi_test$p.value, scientific = TRUE)))
print(paste("Degrees of freedom:", chi_test$parameter))

if(!is.na(chi_test$p.value) && chi_test$p.value < 0.05) {
  print("\n[OK] SIGNIFICANT difference in domain distribution (p < 0.05)")
} else {
  print("\n[INFO] NO significant difference in domain distribution (p >= 0.05)")
}

# ============================================
# 9. FISHER'S EXACT TEST
# ============================================

print("\n========================================")
print("FISHER'S EXACT TEST: PIONEER vs NON-PIONEER")
print("========================================")

fisher_test <- fisher.test(as.data.frame(domain_counts_wide[,-3]), simulate.p.value = TRUE, B = 10000)

print("\n=== FISHER'S EXACT TEST RESULTS ===")
print(fisher_test)
print(paste("P-value:", format(fisher_test$p.value, scientific = TRUE)))

if(fisher_test$p.value < 0.05) {
  print("\n[OK] SIGNIFICANT difference (p < 0.05)")
} else {
  print("\n[INFO] NO significant difference (p >= 0.05)")
}

# ============================================
# 10. DOMAIN ENRICHMENT ANALYSIS
# ============================================

print("\n========================================")
print("DOMAIN ENRICHMENT ANALYSIS")
print("========================================")

domain_category_proportions <- domains_regular %>%
  distinct(Gene, Domain_Consolidated, TF_Category) %>%
  group_by(Domain_Consolidated, TF_Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Domain_Consolidated) %>%
  mutate(Total = sum(Count),
         Proportion = Count / Total) %>%
  ungroup() %>%
  filter(Total >= 3) %>%
  arrange(desc(Proportion))

for(category in c("Pioneer", "NonPioneer")) {
  print(paste("\n--- Most enriched domains in", category, "---"))
  
  enriched <- domain_category_proportions %>%
    filter(TF_Category == category) %>%
    arrange(desc(Proportion)) %>%
    head(15)
  
  print(enriched)
}

# ============================================
# 11. EXPORT RESULTS
# ============================================

write.csv(domain_counts_wide, "contingency_table_chi_square.csv", row.names = FALSE)
write.csv(all_category_results, "top_domains_pioneer_vs_nonpioneer.csv", row.names = FALSE)
write.csv(protein_profiles_expanded, "protein_profiles_expanded.csv", row.names = FALSE)

domain_counts_export <- domain_counts_wide %>%
  arrange(desc(Total))
write.csv(domain_counts_export, "domain_counts_by_tf_category.csv", row.names = FALSE)

write.csv(domain_category_proportions, "domain_enrichment_by_tf_category.csv", row.names = FALSE)

for(category in names(category_plots)) {
  filename <- paste0("domains_", category, ".png")
  ggsave(filename, category_plots[[category]], width = 10, height = 8, dpi = 300)
  print(paste("Saved:", filename))
}

# ============================================
# 12. CREATE COMBINED COMPARISON PLOT
# ============================================

combined_top_domains <- bind_rows(category_results)

p_faceted <- ggplot(combined_top_domains, aes(x = reorder(Domain, Gene_Count), y = Gene_Count, fill = TF_Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~TF_Category, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Pioneer" = "#E74C3C", "NonPioneer" = "#3498DB")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(title = "Top 15 Domains: Pioneer vs Non-Pioneer TFs",
       x = "Domain", 
       y = "Number of Unique Genes",
       fill = "TF Category")

print(p_faceted)
ggsave("domains_pioneer_vs_nonpioneer_faceted.png", p_faceted, width = 14, height = 10, dpi = 300)

print("\n========================================")
print("ANALYSIS COMPLETE")
print("========================================")
print("\nOutput files:")
print("- contingency_table_chi_square.csv")
print("- top_domains_pioneer_vs_nonpioneer.csv")
print("- protein_profiles_expanded.csv")
print("- domain_counts_by_tf_category.csv")
print("- domain_enrichment_by_tf_category.csv")
print("- domains_Pioneer.png")
print("- domains_NonPioneer.png")
print("- domains_pioneer_vs_nonpioneer_faceted.png")
