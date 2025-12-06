# proteindomaininfocddomain.R
# Author: Chris Carson
# Description: Analyzes domain distributions and performs statistical tests.
# Separates superfamilies from regular domains, creates visualizations.

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

# Configuration: Set input files
PIONEER_DOMAINS_FILE <- "Simplified_Domainrecord_PF_pfam.csv"
NONPIONEER_DOMAINS_FILE <- "Simplified_Domainrecord_TF_pfam.csv"

pf <- read.csv(PIONEER_DOMAINS_FILE)
nonpf <- read.csv(NONPIONEER_DOMAINS_FILE)

# ============================================
# 0. SEPARATE SUPERFAMILIES FROM REGULAR DOMAINS
# ============================================

pf_regular <- pf %>% filter(Hit_Type != "superfamily")
nonpf_regular <- nonpf %>% filter(Hit_Type != "superfamily")

# ============================================
# 1. REGULAR DOMAINS ANALYSIS
# ============================================

print("\n========================================")
print("REGULAR DOMAINS ANALYSIS")
print("========================================")

pf_counts_regular <- pf_regular %>%
  distinct(Gene, Domain) %>%
  group_by(Domain) %>%
  summarise(Pioneer_Genes = n()) %>%
  arrange(desc(Pioneer_Genes))

nonpf_counts_regular <- nonpf_regular %>%
  distinct(Gene, Domain) %>%
  group_by(Domain) %>%
  summarise(NonPioneer_Genes = n()) %>%
  arrange(desc(NonPioneer_Genes))

combined_stats_regular <- full_join(pf_counts_regular, nonpf_counts_regular, by = "Domain") %>%
  replace_na(list(Pioneer_Genes = 0, NonPioneer_Genes = 0)) %>%
  mutate(
    Total_Genes = Pioneer_Genes + NonPioneer_Genes,
    Pioneer_Proportion = Pioneer_Genes / Total_Genes
  )

print("\nRegular Domain Statistics:")
print(combined_stats_regular %>% arrange(desc(Pioneer_Proportion)))

# ============================================
# CHI-SQUARE TEST
# ============================================

MIN_GENES_FOR_TEST <- 3
combined_stats_regular_for_test <- combined_stats_regular %>%
  filter(Total_Genes >= MIN_GENES_FOR_TEST)

print(paste("\nDomains used for chi-square test (Total_Genes >=", MIN_GENES_FOR_TEST, "):", nrow(combined_stats_regular_for_test)))

contingency_table_regular <- as.matrix(combined_stats_regular_for_test[, c("Pioneer_Genes", "NonPioneer_Genes")])
rownames(contingency_table_regular) <- combined_stats_regular_for_test$Domain

chi_test_regular <- chisq.test(contingency_table_regular)

print("\n=== CHI-SQUARE TEST RESULTS (REGULAR DOMAINS) ===")
print(chi_test_regular)
print(paste("Chi-square statistic:", round(chi_test_regular$statistic, 4)))
print(paste("P-value:", format(chi_test_regular$p.value, scientific = TRUE)))

if(!is.na(chi_test_regular$p.value) && chi_test_regular$p.value < 0.05) {
  print("\n[OK] SIGNIFICANT difference (p < 0.05)")
} else {
  print("\n[INFO] NO significant difference (p >= 0.05)")
}

if(any(chi_test_regular$expected < 5)) {
  warning("Some expected frequencies < 5. Using Fisher's exact test instead.")
  print("\n=== FISHER'S EXACT TEST (REGULAR DOMAINS) ===")
  fisher_test_regular <- fisher.test(contingency_table_regular, simulate.p.value = TRUE, B = 10000)
  print(fisher_test_regular)
  print(paste("Fisher's P-value:", format(fisher_test_regular$p.value, scientific = TRUE)))
}

# ============================================
# FILTERING FOR VISUALIZATION
# ============================================

MIN_TOTAL_GENES_VIZ <- 5
MIN_NONPIONEER_GENES <- 1

combined_stats_regular_for_viz <- combined_stats_regular %>%
  filter(Total_Genes >= MIN_TOTAL_GENES_VIZ, 
         NonPioneer_Genes >= MIN_NONPIONEER_GENES)

TOP_N_REGULAR <- 10
combined_stats_regular_filtered <- combined_stats_regular_for_viz %>%
  arrange(desc(Pioneer_Proportion)) %>%
  head(TOP_N_REGULAR)

print(paste("\nTop", TOP_N_REGULAR, "domains by pioneer proportion:"))
print(combined_stats_regular_filtered)

# ============================================
# 3. VISUALIZATION - FILTERED REGULAR DOMAINS
# ============================================

plot_data_regular <- combined_stats_regular_filtered %>%
  pivot_longer(cols = c(Pioneer_Genes, NonPioneer_Genes),
               names_to = "TF_Type",
               values_to = "Gene_Count") %>%
  mutate(TF_Type = gsub("_Genes", "", TF_Type))

# Bar plot
p1_regular <- ggplot(plot_data_regular, aes(x = reorder(Domain, -Pioneer_Proportion), y = Gene_Count, fill = TF_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 17),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 10)
  ) +
  labs(title = "Regular Domain Distribution: Pioneer vs Non-Pioneer TFs",
       subtitle = paste("Top", TOP_N_REGULAR, "domains by pioneer proportion"),
       x = "Domain", y = "Number of Unique Genes", fill = "TF Type") +
  scale_fill_manual(values = c("Pioneer" = "#E74C3C", "NonPioneer" = "#3498DB"))

print(p1_regular)

# Proportional stacked bar plot
p2_regular <- ggplot(plot_data_regular, aes(x = reorder(Domain, -Pioneer_Proportion), y = Gene_Count, fill = TF_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 17),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 10)
  ) +
  labs(title = "Regular Domain Proportions: Pioneer vs Non-Pioneer TFs",
       subtitle = paste("Top", TOP_N_REGULAR, "domains by pioneer proportion"),
       x = "Domain", y = "Proportion", fill = "TF Type") +
  scale_fill_manual(values = c("Pioneer" = "#E74C3C", "NonPioneer" = "#3498DB")) +
  scale_y_continuous(labels = scales::percent)

print(p2_regular)

# ============================================
# 4. TOP 20 DOMAINS ANALYSIS
# ============================================

get_top_domains <- function(data, n = 20) {
  top_domains <- data %>%
    distinct(Gene, Domain) %>%
    mutate(Domain_Base = gsub("_\\d+$", "", Domain),
           Domain_Base = gsub("zf-H2C2", "zf-C2H2", Domain_Base)) %>%
    group_by(Domain_Base) %>%
    summarise(Gene_Count = n(), .groups = 'drop') %>%
    rename(Domain = Domain_Base) %>%
    arrange(desc(Gene_Count)) %>%
    head(n)
  
  return(top_domains)
}

top_20_pioneer <- get_top_domains(pf_regular, n = 20)
top_20_nonpioneer <- get_top_domains(nonpf_regular, n = 20)

print("\n========================================")
print("TOP 20 DOMAINS FOR PIONEER FACTORS")
print("========================================")
print(top_20_pioneer)

print("\n========================================")
print("TOP 20 DOMAINS FOR NON-PIONEER FACTORS")
print("========================================")
print(top_20_nonpioneer)

# ============================================
# 5. VISUALIZATION: TOP 20 DOMAINS
# ============================================

plot_pioneer <- ggplot(top_20_pioneer, aes(x = reorder(Domain, Gene_Count), y = Gene_Count)) +
  geom_bar(stat = "identity", fill = "#E74C3C") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(title = "Top 20 Domains in Pioneer TFs",
       x = "Domain", y = "Number of Unique Genes")

print(plot_pioneer)

plot_nonpioneer <- ggplot(top_20_nonpioneer, aes(x = reorder(Domain, Gene_Count), y = Gene_Count)) +
  geom_bar(stat = "identity", fill = "#3498DB") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(title = "Top 20 Domains in Non-Pioneer Transcription Factors",
       x = "Domain", y = "Number of Unique Genes")

print(plot_nonpioneer)

# Export results
write.csv(top_20_pioneer, "top_20_domains_pioneer.csv", row.names = FALSE)
write.csv(top_20_nonpioneer, "top_20_domains_nonpioneer.csv", row.names = FALSE)
write.csv(top_20_comparison, "top_20_domains_comparison.csv", row.names = FALSE)

# Save plots
ggsave("top_20_pioneer_domains.png", plot_pioneer, width = 10, height = 8, dpi = 300)
ggsave("top_20_nonpioneer_domains.png", plot_nonpioneer, width = 10, height = 8, dpi = 300)
