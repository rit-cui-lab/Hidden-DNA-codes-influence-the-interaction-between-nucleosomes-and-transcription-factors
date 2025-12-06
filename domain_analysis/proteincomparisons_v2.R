# proteincomparisons_v2.R
# Author: Chris Carson
# Description: Comprehensive analysis of protein domains across chromatin profiles.
# Analyzes both stringent (no ties) and non-stringent datasets, creates visualizations.

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

pf <- read.csv(PIONEER_DOMAINS_FILE)
nonpf <- read.csv(NONPIONEER_DOMAINS_FILE)

all_domains <- bind_rows(
  pf %>% mutate(TF_Category = "Pioneer"),
  nonpf %>% mutate(TF_Category = "NonPioneer")
)

protein_profiles <- read.csv(PROTEIN_PROFILES_FILE)

print("========================================")
print("PROTEIN PROFILE FILTERING")
print("========================================")
print(paste("Total proteins with profiles:", nrow(protein_profiles)))

# Create stringent dataset - only proteins without ties (no commas)
protein_profiles_stringent <- protein_profiles %>%
  filter(!grepl(",", Consensus))

print(paste("Proteins with NO ties (stringent):", nrow(protein_profiles_stringent)))
print(paste("Proteins removed (had ties):", nrow(protein_profiles) - nrow(protein_profiles_stringent)))

print("\nStringent profiles (no ties):")
print(table(protein_profiles_stringent$Consensus))

# ============================================
# FUNCTION TO ANALYZE PROFILES
# ============================================

analyze_profiles <- function(protein_profile_data, analysis_label) {
  
  print(paste("\n########################################"))
  print(paste("ANALYSIS:", analysis_label))
  print(paste("########################################"))
  
  # Expand profiles if needed
  if(analysis_label == "ALL (Non-Stringent)") {
    protein_profiles_expanded <- protein_profile_data %>%
      mutate(Profile_List = str_split(Consensus, ",\\s*")) %>%
      unnest(Profile_List) %>%
      rename(Individual_Profile = Profile_List) %>%
      mutate(Consensus = Individual_Profile) %>%
      select(Transcription_Factor, Individual_Profile, Consensus)
  } else {
    protein_profiles_expanded <- protein_profile_data %>%
      rename(Individual_Profile = Consensus) %>%
      mutate(Consensus = Individual_Profile) %>%
      select(Transcription_Factor, Individual_Profile, Consensus)
  }
  
  print(paste("Protein-profile combinations:", nrow(protein_profiles_expanded)))
  
  # Join domains with profiles
  domains_with_profiles <- all_domains %>%
    left_join(protein_profiles_expanded, by = c("Gene" = "Transcription_Factor"),
              relationship = "many-to-many")
  
  # Remove proteins without profile assignments
  domains_with_profiles <- domains_with_profiles %>%
    filter(!is.na(Individual_Profile))
  
  print(paste("Domain records with profiles:", nrow(domains_with_profiles)))
  print(paste("Unique proteins:", n_distinct(domains_with_profiles$Gene)))
  
  # Filter to regular domains only
  domains_regular <- domains_with_profiles %>% 
    filter(Hit_Type != "superfamily")
  
  # Consolidate domain variants
  domains_regular <- domains_regular %>%
    mutate(
      Domain_Clean = gsub("_[0-9]+$", "", Domain),
      Domain_Clean = gsub("H2C2", "C2H2", Domain_Clean)
    )
  
  print(paste("Regular domain entries:", nrow(domains_regular)))
  print(paste("Unique domains (before combining):", n_distinct(domains_regular$Domain)))
  print(paste("Unique domains (after combining variants):", n_distinct(domains_regular$Domain_Clean)))
  
  # Show examples of domain consolidation
  domain_mapping <- domains_regular %>%
    distinct(Domain, Domain_Clean) %>%
    arrange(Domain_Clean, Domain)
  
  print("\nExample domain consolidation (first 30):")
  print(head(domain_mapping, 30))
  
  # Show specifically the zinc finger consolidation
  zf_domains <- domain_mapping %>%
    filter(grepl("zf", Domain_Clean, ignore.case = TRUE))
  
  if(nrow(zf_domains) > 0) {
    print("\nZinc finger domain consolidation:")
    print(zf_domains)
  }
  
  # Get list of all individual profiles
  all_profiles <- unique(domains_regular$Individual_Profile)
  all_profiles <- all_profiles[order(all_profiles)]
  
  print("\nProfiles found:")
  print(all_profiles)
  
  # Store results for each profile
  profile_results <- list()
  profile_plots <- list()
  
  TOP_N_PER_PROFILE <- 10
  
  for(profile in all_profiles) {
    
    print(paste("\n--- Profile:", profile, "---"))
    
    # Filter domains for this profile
    profile_domains <- domains_regular %>%
      filter(Individual_Profile == profile)
    
    # Count unique genes per domain
    domain_counts <- profile_domains %>%
      distinct(Gene, Domain_Clean) %>%
      group_by(Domain_Clean) %>%
      summarise(Gene_Count = n()) %>%
      arrange(desc(Gene_Count))
    
    print(paste("Unique proteins:", n_distinct(profile_domains$Gene)))
    print(paste("Unique domains (after combining):", nrow(domain_counts)))
    
    # Get top N domains
    top_domains <- domain_counts %>%
      head(TOP_N_PER_PROFILE)
    
    # Store results
    top_domains$Profile <- profile
    top_domains$Analysis <- analysis_label
    profile_results[[profile]] <- top_domains
    
    # Create plot for this profile
    if(nrow(top_domains) > 0) {
      plot_title <- ifelse(analysis_label == "STRINGENT",
                           paste("Top Domains in", profile, "Profile (Stringent)"),
                           paste("Top Domains in", profile, "Profile"))
      
      plot_subtitle <- ifelse(analysis_label == "STRINGENT",
                              paste(n_distinct(profile_domains$Gene), "proteins (no consensus ties)"),
                              paste(n_distinct(profile_domains$Gene), "proteins"))
      
      p <- ggplot(top_domains, aes(x = reorder(Domain_Clean, -Gene_Count), y = Gene_Count)) +
        geom_bar(stat = "identity", fill = "#3498DB") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 34),
          axis.text.y = element_text(size = 32),
          axis.title.x = element_text(size = 40),
          axis.title.y = element_text(size = 40),
          plot.title = element_text(size = 40, face = "bold"),
          plot.subtitle = element_text(size = 28)
        ) +
        labs(title = plot_title,
             subtitle = plot_subtitle,
             x = "Domain", 
             y = "Number of Unique Genes")
      
      profile_plots[[profile]] <- p
      print(p)
    }
  }
  
  # Combine all profile results
  all_profile_results <- bind_rows(profile_results)
  
  # Create domain counts wide format
  domain_counts_wide <- domains_regular %>%
    distinct(Gene, Domain_Clean, Individual_Profile) %>%
    group_by(Domain_Clean, Individual_Profile) %>%
    summarise(Gene_Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Individual_Profile, 
                values_from = Gene_Count, 
                values_fill = 0) %>%
    mutate(Total = rowSums(select(., -Domain_Clean))) %>%
    arrange(desc(Total))
  
  # Create enrichment analysis
  domain_profile_proportions <- domains_regular %>%
    distinct(Gene, Domain_Clean, Individual_Profile) %>%
    group_by(Domain_Clean, Individual_Profile) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Domain_Clean) %>%
    mutate(Total = sum(Count),
           Proportion = Count / Total) %>%
    ungroup() %>%
    filter(Total >= 3) %>%
    arrange(desc(Proportion))
  
  # Return results
  return(list(
    profile_results = all_profile_results,
    profile_plots = profile_plots,
    domain_counts_wide = domain_counts_wide,
    domain_enrichment = domain_profile_proportions,
    protein_profiles_used = protein_profiles_expanded,
    all_profiles = all_profiles,
    domain_mapping = domain_mapping
  ))
}

# ============================================
# RUN BOTH ANALYSES
# ============================================

results_all <- analyze_profiles(protein_profiles, "ALL (Non-Stringent)")
results_stringent <- analyze_profiles(protein_profiles_stringent, "STRINGENT")

# ============================================
# EXPORT RESULTS - NON-STRINGENT
# ============================================

print("\n========================================")
print("EXPORTING NON-STRINGENT RESULTS")
print("========================================")

write.csv(results_all$profile_results, "top_domains_per_profile.csv", row.names = FALSE)
write.csv(results_all$protein_profiles_used, "protein_profiles_expanded.csv", row.names = FALSE)
write.csv(results_all$domain_counts_wide, "domain_counts_by_profile_wide.csv", row.names = FALSE)
write.csv(results_all$domain_enrichment, "domain_enrichment_by_profile.csv", row.names = FALSE)
write.csv(results_all$domain_mapping, "domain_name_mapping.csv", row.names = FALSE)

for(profile in names(results_all$profile_plots)) {
  clean_name <- gsub("/", "_", profile)
  clean_name <- gsub(" ", "_", clean_name)
  clean_name <- gsub(",", "", clean_name)
  
  filename <- paste0("domains_", clean_name, ".png")
  ggsave(filename, results_all$profile_plots[[profile]], width = 14, height = 12, dpi = 300)
}

p_faceted_all <- ggplot(results_all$profile_results, 
                        aes(x = reorder(Domain_Clean, -Gene_Count), y = Gene_Count)) +
  geom_bar(stat = "identity", fill = "#3498DB") +
  facet_wrap(~Profile, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
    axis.text.y = element_text(size = 26),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    strip.text = element_text(size = 28, face = "bold"),
    plot.title = element_text(size = 36, face = "bold"),
    plot.subtitle = element_text(size = 26)
  ) +
  labs(title = "Top 10 Domains per Chromatin Profile",
       subtitle = "All proteins (including those with consensus ties)",
       x = "Domain", 
       y = "Number of Unique Genes")

ggsave("domains_all_profiles_faceted.png", p_faceted_all, width = 20, height = 16, dpi = 300)

# ============================================
# EXPORT RESULTS - STRINGENT
# ============================================

print("\n========================================")
print("EXPORTING STRINGENT RESULTS")
print("========================================")

write.csv(results_stringent$profile_results, "top_domains_per_profile_STRINGENT.csv", row.names = FALSE)
write.csv(results_stringent$protein_profiles_used, "protein_profiles_expanded_STRINGENT.csv", row.names = FALSE)
write.csv(results_stringent$domain_counts_wide, "domain_counts_by_profile_wide_STRINGENT.csv", row.names = FALSE)
write.csv(results_stringent$domain_enrichment, "domain_enrichment_by_profile_STRINGENT.csv", row.names = FALSE)
write.csv(results_stringent$domain_mapping, "domain_name_mapping_STRINGENT.csv", row.names = FALSE)

for(profile in names(results_stringent$profile_plots)) {
  clean_name <- gsub("/", "_", profile)
  clean_name <- gsub(" ", "_", clean_name)
  clean_name <- gsub(",", "", clean_name)
  
  filename <- paste0("domains_", clean_name, "_STRINGENT.png")
  ggsave(filename, results_stringent$profile_plots[[profile]], width = 14, height = 12, dpi = 300)
}

p_faceted_stringent <- ggplot(results_stringent$profile_results, 
                              aes(x = reorder(Domain_Clean, -Gene_Count), y = Gene_Count)) +
  geom_bar(stat = "identity", fill = "#E74C3C") +
  facet_wrap(~Profile, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
    axis.text.y = element_text(size = 26),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    strip.text = element_text(size = 28, face = "bold"),
    plot.title = element_text(size = 36, face = "bold"),
    plot.subtitle = element_text(size = 26)
  ) +
  labs(title = "Top 10 Domains per Chromatin Profile (Stringent)",
       subtitle = "Only proteins with no consensus ties",
       x = "Domain", 
       y = "Number of Unique Genes")

ggsave("domains_all_profiles_faceted_STRINGENT.png", p_faceted_stringent, width = 20, height = 16, dpi = 300)

# ============================================
# COMPARISON SUMMARY
# ============================================

print("\n========================================")
print("ANALYSIS COMPARISON SUMMARY")
print("========================================")

comparison <- data.frame(
  Analysis = c("Non-Stringent", "Stringent"),
  Proteins_With_Profiles = c(nrow(protein_profiles), nrow(protein_profiles_stringent)),
  Profiles_Found = c(length(results_all$all_profiles), length(results_stringent$all_profiles)),
  Total_Domain_Records = c(nrow(results_all$domain_counts_wide), nrow(results_stringent$domain_counts_wide))
)

print(comparison)
write.csv(comparison, "analysis_comparison.csv", row.names = FALSE)

print("\n========================================")
print("ANALYSIS COMPLETE")
print("========================================")
print("\nNON-STRINGENT FILES:")
print("- top_domains_per_profile.csv")
print("- protein_profiles_expanded.csv")
print("- domain_counts_by_profile_wide.csv")
print("- domain_enrichment_by_profile.csv")
print("- domain_name_mapping.csv")
print("- domains_[PROFILE].png (individual plots)")
print("- domains_all_profiles_faceted.png")
print("\nSTRINGENT FILES:")
print("- top_domains_per_profile_STRINGENT.csv")
print("- protein_profiles_expanded_STRINGENT.csv")
print("- domain_counts_by_profile_wide_STRINGENT.csv")
print("- domain_enrichment_by_profile_STRINGENT.csv")
print("- domain_name_mapping_STRINGENT.csv")
print("- domains_[PROFILE]_STRINGENT.png (individual plots)")
print("- domains_all_profiles_faceted_STRINGENT.png")
print("- analysis_comparison.csv")
