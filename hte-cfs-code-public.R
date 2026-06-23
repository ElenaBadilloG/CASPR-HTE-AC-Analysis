

library(haven)
library(tidyverse)
library(twang)
library(MatchIt)
library(cobalt)
library(bcf)
library(survival)
library(gtsummary)
library(WeightIt)
library(dplyr)
library(grf)
library(fpc)
library(haven)
library(DescTools)
library(cluster)
library(survminer)



options(scipen=999)
set.seed(427)


theme_tuf2 <- ggthemes::theme_tufte() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.title.x = element_text(vjust = -0.3), 
    axis.title.y = element_text(vjust = 0.8),
    axis.text.x = element_text(angle = 90),
    legend.background = element_blank(), 
    legend.key = element_blank(), 
    legend.title = element_text(face="plain"),
    panel.background = element_blank(), 
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank(),
    
    plot.subtitle = element_text(hjust = 0.75, vjust= -1),
    
  )




dfc <- read_dta("data.dta")

dfc <- dplyr::filter(dfc, include==1)
dfc$tx <- ifelse(dfc$tx_ap_no_ac ==1, 0,1)

dfc$race <- as.factor(dfc$race)
dfc$sex <- as.factor(dfc$sex)
colnames(dfc)
dfc$valvular_lesion <- as.factor(dfc$valvular_lesion)
dfc$insurance <- as.factor(dfc$insurance_primary2)

################
orig_vars <-  c(
  "lv_injury",
  "mrs_pta",             
  "age",                  
  "sex",                  
  "race",                 
  "hispanic",            
  "insurance",
  "hx_htn",              
  "hx_hld",              
  "hx_tob",              
  "hx_dm",               
  "hx_cad",              
  "hx_stroke",           
  "hx_cancer",           
  "lae_severe",       
  "pfo",                 
  "lvef",                
  "wma",
  "mri_dwi_multifocal",  
  "dwi_cortical",        
  "bnihss")

clustering_vars <- orig_vars

variables <- unique(c(
  clustering_vars,
  "lv_injury",
  "mrs_pta",             
  "age",                  
  "sex",                  
  "race",                 
  "hispanic",            
  "insurance",
  "hx_htn",              
  "hx_hld",              
  "hx_tob",              
  "hx_dm",               
  "hx_cad",              
  "hx_stroke",           
  "hx_cancer",           
  "lae_severe",       
  "wma",
  "lvef",                
  
  "mri_dwi_multifocal",  
  "dwi_cortical",        
  "tx_ap_no_ac",         
  "bnihss",              
  "time_to_event", "stroke_bleed_death2", "stroke_bleed_death2"
))

####
# Pull raw values before scaling
dfc_raw <- dfc %>%
  select(all_of(variables)) %>%
  mutate(complete_case = complete.cases(across(all_of(variables)))) %>%
  mutate(complete_case = factor(complete_case,
                                levels = c(TRUE, FALSE),
                                labels = c("Included", "Excluded")))

dfc_raw <- dfc %>%
  select(all_of(variables)) %>%
  mutate(complete_case = complete.cases(across(all_of(variables)))) %>%
  mutate(
    complete_case = factor(complete_case, levels = c(TRUE, FALSE),
                           labels = c("Included", "Excluded")),
    sex = factor(sex, levels = c(0, 1), labels = c("Male", "Female")),
    race = factor(race, levels = c(1, 2, 3, 4),
                  labels = c("White", "Black", "Asian/Pacific Islander", "Other")),
    insurance = factor(insurance, levels = c(0, 1, 2, 3, 4),
                       labels = c("None", "Private/Veterans Affairs", "Medicare",
                                  "Other", "Medicaid/None/Pending"))
  )

tbl_missing <- dfc_raw %>%
  select(all_of(orig_vars), complete_case) %>%
  tbl_summary(
    by = complete_case,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      lv_injury       ~ "LV Injury",
      mrs_pta         ~ "mRS Prior to Admission",
      age             ~ "Age",
      sex             ~ "Sex",
      race            ~ "Race",
      hispanic        ~ "Hispanic Ethnicity",
      insurance       ~ "Insurance Type",
      hx_htn          ~ "History of Hypertension",
      hx_hld          ~ "History of Hyperlipidemia",
      hx_tob          ~ "History of Tobacco Use",
      hx_dm           ~ "History of Diabetes",
      hx_cad          ~ "History of CAD",
      hx_stroke       ~ "History of Stroke",
      hx_cancer       ~ "History of Cancer",
      lae_severe      ~ "Severe Left Atrial Enlargement",
      pfo             ~ "Patent Foramen Ovale",
      lvef            ~ "Left Ventricular Ejection Fraction",
      wma             ~ "Wall Motion Abnormality",
      mri_dwi_multifocal ~ "Multifocal DWI Lesions",
      dwi_cortical    ~ "Cortical DWI Lesions",
      bnihss          ~ "Baseline NIHSS"
    )
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  bold_p(t = 0.05) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Complete Case Status**")

print(tbl_missing)

flextable::save_as_docx(as_flex_table(tbl_missing),
                        path = "t1_missing_comparison.docx")


####

# rescale
# Min-Max normalization function
minmax_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Apply to your variables
dfc$hdl <- minmax_scale(dfc$hdl)
dfc$ldl <- minmax_scale(dfc$ldl)
dfc$trig <- minmax_scale(dfc$trig)
dfc$lvef <- minmax_scale(dfc$lvef)
dfc$age <- minmax_scale(dfc$age)
dfc$bnihss <- minmax_scale(dfc$bnihss)
dfc$hba1c <- minmax_scale(dfc$hba1c)
dfc$mri_fazekas <- minmax_scale(dfc$mri_fazekas)
dfc$mrs_pta <- minmax_scale(dfc$mrs_pta)

dfc$mri_microhemorrhage_count <- ifelse(dfc$mri_microhemorrhage==0, 0, dfc$mri_microhemorrhage_count)






demog_vars <-  c("age", "sex", "race", "hispanic", "insurance")
mri_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "mri_") &  ! str_detect(x, "72") &
    ! str_detect(x, "23") & ! str_detect(x, "_loc")}))]
cta_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "cta_")}))]
tee_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){x=="tee"}))]
tte_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "tte_la_diameter")}))]
hx_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "hx_") &
    ! str_detect(x, "_type") &
    ! str_detect(x, "cancer_met") &
    ! str_detect(x, "insulin") &
    ! str_detect(x, "timing") &
    ! str_detect(x, "carotid_stenosis")}))]
tx_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "tx_")}))]
time2_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "time_")}))]
wma_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "wma_")}))]
#lesion_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "lesion_")}))]
tcd_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "tcd_")}))]
recurr_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "recurr")}))]
lvo_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){x=="lvo"}))]
lv_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "lv_injury")}))]
lae_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "lae_")}))]
pta_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "pta_")}))]
pfo_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){x=="pfo"}))]
t2_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "t2_")}))]
dwi_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "dwi_")}))]
mlan_vars <- c("mlakinesis", "mlhypokinesis", "mldyskinesis", "mlaneurysm", "mlreportwma")

lablev_vars <- c("hba1c", "ldl", "hdl", "trig")

strkchar_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){(x=="mrs_pta" |
                                                                           str_detect(x, "nihss")|
                                                                           str_detect(x, "etiol")| str_detect(x, "aspects")) & (! str_detect(x, "_90d")) &
    (! str_detect(x, "_cat"))}))]

death_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "death")}))]
bleed_vars <- colnames(dfc)[unlist(lapply(colnames(dfc), function(x){str_detect(x, "isth_bleed")}))]

clustering_vars  <-  unique(c(demog_vars, strkchar_vars, mri_vars, cta_vars, tee_vars, tte_vars, hx_vars,
                              wma_vars, tcd_vars,  lvo_vars, lv_vars, lae_vars, pta_vars,
                              pfo_vars, t2_vars, dwi_vars, lablev_vars, "lvef"))
outc_vars  <-  unique(c("time_to_event", bleed_vars, death_vars, "los", "stroke_bleed_death2",
                        "stroke_bleed_death2", "stroke_bleed_death", "recurrent_ischemic_stroke2", "death", "isth_bleed2"))




dfs <- select(dfc, clustering_vars, tx_vars, outc_vars)

dfs$wma <- ifelse(dfs$wma_apical==1 | dfs$wma_basal==1 | dfs$wma_mid ==1, 1,0)


################
orig_vars <-  c(
  "lv_injury",
  "mrs_pta",             
  "age",                  
  "sex",                  
  "race",                 
  "hispanic",            
  "insurance",
  "hx_htn",              
  "hx_hld",              
  "hx_tob",              
  "hx_dm",               
  "hx_cad",              
  "hx_stroke",           
  "hx_cancer",           
  "lae_severe",       
  "pfo",                 
  "lvef",                
  "wma",
  "mri_dwi_multifocal",  
  "dwi_cortical",        
  "bnihss")

clustering_vars <- orig_vars

variables <- unique(c(
  clustering_vars,
  "lv_injury",
  "mrs_pta",             
  "age",                  
  "sex",                  
  "race",                 
  "hispanic",            
  "insurance",
  "hx_htn",              
  "hx_hld",              
  "hx_tob",              
  "hx_dm",               
  "hx_cad",              
  "hx_stroke",           
  "hx_cancer",           
  "lae_severe",       
  "wma",
  "lvef",                
  
  "mri_dwi_multifocal",  
  "dwi_cortical",        
  "tx_ap_no_ac",         
  "bnihss",              
  "time_to_event", "stroke_bleed_death2", "stroke_bleed_death2"
))



# CASPR Causal Survival Forest Analysis
# Adapted from Sverdrup & Wager (2024) methodology

# Select relevant columns
dfs <- na.omit(select(dfs, all_of(variables)))



# Create data subset for clustering
cluster_data <- dfs %>% dplyr::select(all_of(orig_vars))


# Method 1: Convert haven_labelled to appropriate R types
cluster_data <- cluster_data %>%
  # Convert haven_labelled to numeric or factor as appropriate
  mutate(across(where(is.labelled), as_factor)) %>%  # Convert labelled to factors
  mutate(across(where(is.character), as.factor)) %>%  # Convert character to factors
  # For binary variables that should be numeric (0/1), convert back
  mutate(across(where(~is.factor(.) && length(levels(.)) == 2), as.numeric)) %>%
  # Subtract 1 from binary factors to get 0/1 coding
  mutate(across(where(~is.numeric(.) && all(. %in% c(1,2), na.rm = TRUE)), ~. - 1))




# Fixed Cramer's V function with better missing value handling
cramers_v <- function(x, y) {
  # Convert to factors if not already
  x <- as.factor(x)
  y <- as.factor(y)
  
  
  # Check if we have variation in both variables
  if (length(unique(x)) <= 1 || length(unique(y)) <= 1) {
    return(0)
  }
  
  # Create contingency table
  tbl <- table(x, y)
  
  # Calculate Cramer's V
  chi_sq_test <- tryCatch({
    chisq.test(tbl)
  }, error = function(e) {
    return(list(statistic = 0))
  })
  
  chi_sq <- chi_sq_test$statistic
  n <- sum(tbl)
  min_dim <- min(dim(tbl)) - 1
  
  if (min_dim == 0 || n == 0) return(0)
  
  cramers_v <- sqrt(chi_sq / (n * min_dim))
  return(as.numeric(cramers_v))
}

# Calculate correlation matrix for categorical variables
calc_categorical_correlations <- function(data, cat_vars) {
  n_vars <- length(cat_vars)
  corr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  rownames(corr_matrix) <- colnames(corr_matrix) <- cat_vars
  
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) {
        corr_matrix[i, j] <- 1
      } else if (i < j) {
        corr_matrix[i, j] <- cramers_v(data[[cat_vars[i]]], data[[cat_vars[j]]])
        corr_matrix[j, i] <- corr_matrix[i, j]  # Symmetric
      }
    }
  }
  
  return(corr_matrix)
}

# Function to calculate correlation between numerical variable and categorical variables
calculate_mixed_correlations <- function(data, cat_vars, num_var) {
  
  mixed_correlations <- numeric(length(cat_vars))
  names(mixed_correlations) <- cat_vars
  
  cat(sprintf("Calculating mixed correlations between '%s' and categorical variables...\n", num_var))
  
  for (cat_var in cat_vars) {
    tryCatch({
      # Use eta-squared (effect size) for categorical-numerical association
      # This measures how much of the variance in the numerical variable is explained by the categorical variable
      
      # Create model with categorical variable predicting numerical variable
      model <- lm(data[[num_var]] ~ data[[cat_var]], data = data)
      
      # Calculate eta-squared (proportion of variance explained)
      ss_total <- sum((data[[num_var]] - mean(data[[num_var]], na.rm = TRUE))^2, na.rm = TRUE)
      ss_residual <- sum(residuals(model)^2, na.rm = TRUE)
      ss_explained <- ss_total - ss_residual
      
      # Eta-squared = explained variance / total variance
      eta_squared <- ss_explained / ss_total
      
      # Convert to correlation-like measure (0 to 1)
      mixed_correlations[cat_var] <- sqrt(abs(eta_squared))
      
    }, error = function(e) {
      # If model fails, set correlation to 0
      mixed_correlations[cat_var] <<- 0
      cat(sprintf("Warning: Could not calculate correlation for %s: %s\n", cat_var, e$message))
    })
  }
  
  return(mixed_correlations)
}

# Enhanced correlation weights function that handles mixed types
calculate_correlation_weights_mixed <- function(data, cat_vars, num_vars = NULL, 
                                                weight_method = "inverse_max") {
  
  # Check data
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame")
  }
  
  # Initialize weights
  all_vars <- c(cat_vars, num_vars)
  weights <- rep(1, length(all_vars))
  names(weights) <- all_vars
  
  # Calculate categorical correlations (Cramer's V)
  if (length(cat_vars) > 1) {
    cat("Calculating categorical correlations (Cramer's V)...\n")
    
    # Print diagnostics for categorical variables
    cat("Categorical variables:\n")
    for (var in cat_vars) {
      n_levels <- length(unique(data[[var]][!is.na(data[[var]])]))
      n_missing <- sum(is.na(data[[var]]))
      cat(sprintf("  %s: %d levels, %d missing values\n", var, n_levels, n_missing))
    }
    
    cat_corr <- calc_categorical_correlations(data, cat_vars)
    
    cat("Categorical correlation matrix (Cramer's V):\n")
    print(round(cat_corr, 3))
    
    # Calculate weights for categorical variables
    for (var in cat_vars) {
      correlations <- cat_corr[var, cat_corr[var, ] != 1]  # Exclude self-correlation
      
      if (length(correlations) > 0) {
        max_corr <- max(abs(correlations))
        mean_corr <- mean(abs(correlations))
        
        if (weight_method == "inverse_max") {
          weights[var] <- 1 / (1 + max_corr)
        } else if (weight_method == "inverse_mean") {
          weights[var] <- 1 / (1 + mean_corr)
        } else if (weight_method == "correlation_penalty") {
          weights[var] <- max(0.1, 1 - max_corr)
        }
      }
    }
  }
  
  # Handle numerical variables
  if (!is.null(num_vars)) {
    
    cat(sprintf("\nNumerical variables: %s\n", paste(num_vars, collapse = ", ")))
    
    if (length(num_vars) == 1) {
      # Single numerical variable - calculate its association with categorical variables
      num_var <- num_vars[1]
      
      cat(sprintf("Single numerical variable '%s': calculating associations with categorical variables\n", num_var))
      
      # Print diagnostics for numerical variable
      n_missing <- sum(is.na(data[[num_var]]))
      var_range <- range(data[[num_var]], na.rm = TRUE)
      cat(sprintf("  %s: range [%.2f, %.2f], %d missing values\n", 
                  num_var, var_range[1], var_range[2], n_missing))
      
      # Calculate mixed correlations
      mixed_corr <- calculate_mixed_correlations(data, cat_vars, num_var)
      
      cat("\nCategorical-Numerical associations (eta coefficient):\n")
      print(round(mixed_corr, 3))
      
      # Calculate weight for numerical variable based on its associations
      if (length(mixed_corr) > 0 && any(!is.na(mixed_corr))) {
        max_mixed_corr <- max(abs(mixed_corr), na.rm = TRUE)
        mean_mixed_corr <- mean(abs(mixed_corr), na.rm = TRUE)
        
        if (weight_method == "inverse_max") {
          weights[num_var] <- 1 / (1 + max_mixed_corr)
        } else if (weight_method == "inverse_mean") {
          weights[num_var] <- 1 / (1 + mean_mixed_corr)
        } else if (weight_method == "correlation_penalty") {
          weights[num_var] <- max(0.1, 1 - max_mixed_corr)
        }
        
        cat(sprintf("Numerical variable weight based on max association (%.3f): %.3f\n", 
                    max_mixed_corr, weights[num_var]))
      } else {
        weights[num_var] <- 1.0
        cat("No significant associations found, using default weight 1.0\n")
      }
      
    } else if (length(num_vars) > 1) {
      # Multiple numerical variables - calculate correlations between them
      cat("Calculating numerical correlations (Pearson)...\n")
      
      for (var in num_vars) {
        n_missing <- sum(is.na(data[[var]]))
        var_range <- range(data[[var]], na.rm = TRUE)
        cat(sprintf("  %s: range [%.2f, %.2f], %d missing values\n", 
                    var, var_range[1], var_range[2], n_missing))
      }
      
      num_data <- data[num_vars]
      num_corr <- abs(cor(num_data, use = "pairwise.complete.obs"))
      
      cat("Numerical correlation matrix (Pearson):\n")
      print(round(num_corr, 3))
      
      # Calculate weights for numerical variables
      for (var in num_vars) {
        correlations <- num_corr[var, num_corr[var, ] != 1]
        
        if (length(correlations) > 0 && !all(is.na(correlations))) {
          max_corr <- max(abs(correlations), na.rm = TRUE)
          mean_corr <- mean(abs(correlations), na.rm = TRUE)
          
          if (weight_method == "inverse_max") {
            weights[var] <- 1 / (1 + max_corr)
          } else if (weight_method == "inverse_mean") {
            weights[var] <- 1 / (1 + mean_corr)
          } else if (weight_method == "correlation_penalty") {
            weights[var] <- max(0.1, 1 - max_corr)
          }
        }
      }
      
      # Also calculate mixed correlations for each numerical variable
      for (num_var in num_vars) {
        mixed_corr <- calculate_mixed_correlations(data, cat_vars, num_var)
        max_mixed_corr <- max(abs(mixed_corr), na.rm = TRUE)
        
        # Combine numerical and mixed correlations for final weight
        current_weight <- weights[num_var]
        mixed_penalty <- 1 / (1 + max_mixed_corr)
        
        # Take the minimum weight (most restrictive)
        weights[num_var] <- min(current_weight, mixed_penalty)
        
        cat(sprintf("%s: combined weight considering numerical (%.3f) and mixed (%.3f) correlations: %.3f\n",
                    num_var, current_weight, mixed_penalty, weights[num_var]))
      }
    }
  }
  
  return(weights)
}

clean_names <- function(names) {
  # Replace spaces, slashes, and other special characters with underscores
  cleaned <- gsub("[^A-Za-z0-9]", "_", names)
  # Remove multiple consecutive underscores
  cleaned <- gsub("_{2,}", "_", cleaned)
  # Remove leading/trailing underscores
  cleaned <- gsub("^_|_$", "", cleaned)
  return(cleaned)
}
# Updated main function using mixed correlations
calculate_weighted_gower_mixed <- function(data, cat_vars, num_vars = NULL, 
                                           weight_method = "inverse_max") {
  
  # Calculate weights using mixed correlation approach
  weights <- calculate_correlation_weights_mixed(data, cat_vars, num_vars, weight_method)
  
  # Print weights for inspection
  cat("\nFinal feature weights:\n")
  print(round(weights, 3))
  
  # Prepare data for daisy (only include the variables we're clustering on)
  cluster_vars <- c(cat_vars, num_vars)
  clustering_data <- data[cluster_vars]
  
  # Calculate weighted Gower distance
  weighted_dist <- daisy(clustering_data, 
                         metric = "gower", 
                         weights = weights)
  
  return(weighted_dist)
}

num_vars <- c("age", "mrs_pta", "bnihss","lvef")


# Get categorical variables as the complement (everything else)
cat_vars <- setdiff(clustering_vars, num_vars)


# Clean the names
cat_vars <- clean_names(cat_vars)

colnames(cluster_data) <- clean_names(colnames(cluster_data))

# Calculate weighted Gower distance with mixed correlations
wgower_dist <- calculate_weighted_gower_mixed(
  data = cluster_data, 
  cat_vars = cat_vars,
  num_vars = num_vars,
  weight_method = "inverse_max"
)

#####


# Determine optimal k
silhouette_scores <- numeric(9)
for(k in 2:10) {
  pam_fit <- pam(wgower_dist, k = k, diss = TRUE)
  silhouette_scores[k-1] <- pam_fit$silinfo$avg.width
}

# Plot silhouette scores
plot(2:10, silhouette_scores, type = "b",
     xlab = "Number of clusters", ylab = "Average Silhouette Width",
     main = "Optimal Number of Clusters for PAM")

sil_df <- data.frame(k = 2:10, silhouette = silhouette_scores)

s1 <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  # Red rectangle highlighting optimal k=3
  annotate("rect", 
           xmin = 2.7, xmax = 3.3, 
           ymin = min(silhouette_scores) - 0.005, 
           ymax = max(silhouette_scores) + 0.005,
           fill = NA, color = "darkred", linewidth = 1.5) +
  annotate("text", x = 3, y = max(silhouette_scores) + 0.01,
           label = "Optimal (k = 3)", color = "darkred", fontface = "bold", size = 4, vjust = -0.5) +
  # Line and points
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  # Force all integers on x-axis
  scale_x_continuous(breaks = 2:10) +
  labs(x = "Number of Clusters (k)",
       y = "Average Silhouette Width",
       title = "Optimal Number of Clusters for PAM") +
  theme_minimal()

pdf("fS1.pdf", width = 8, height = 5)
print(s1)
dev.off()





# Fit PAM with optimal k
optimal_k <- which.max(silhouette_scores) + 1
pam_fit <- pam(wgower_dist, k = optimal_k, diss = TRUE)

# Add clusters to data
dfs$cluster_pam <- pam_fit$clustering

# View medoids (representative patients for each cluster)
medoid_indices <- pam_fit$medoids
dfs[medoid_indices, ]

# Add cluster assignment to original data
dfs$cluster <- as.factor(dfs$cluster_pam)



# Visualize clusters using MDS
library(Rtsne)
mds <- cmdscale(wgower_dist)

set.seed(123)

perplexity_value <- min(30, floor(nrow(cluster_data)/3) - 1)

# Run multiple t-SNE with the same parameters and compare
num_runs <- 100
tsne_results <- list()

for (i in 1:num_runs) {
  tsne_results[[i]] <- Rtsne(wgower_dist, 
                             is_distance = TRUE, 
                             perplexity = perplexity_value,
                             theta = 0.3,  # Lower values = more accurate but slower
                             dims = 2,
                             max_iter = 1000,
                             verbose = TRUE,
                             Y_init = NULL)  # Random initialization
}


# Function to calculate pairwise correlations between runs
calculate_stability <- function(tsne_list) {
  n_runs <- length(tsne_list)
  stability_scores <- numeric(n_runs)
  
  for (i in 1:n_runs) {
    coords_i <- tsne_list[[i]]$Y
    correlations <- numeric(n_runs - 1)
    counter <- 1
    
    for (j in 1:n_runs) {
      if (i != j) {
        coords_j <- tsne_list[[j]]$Y
        
        # Calculate correlation between coordinate sets
        cor_x <- cor(coords_i[,1], coords_j[,1])
        cor_y <- cor(coords_i[,2], coords_j[,2])
        
        # Average correlation (considering absolute values since t-SNE can flip axes)
        correlations[counter] <- (abs(cor_x) + abs(cor_y)) / 2
        counter <- counter + 1
      }
    }
    
    # Average correlation with all other runs
    stability_scores[i] <- mean(correlations)
  }
  
  # Return the index of the most stable run
  return(which.max(stability_scores))
}

# Find the most stable t-SNE run
most_stable_index <- calculate_stability(tsne_results)
most_stable_tsne <- tsne_results[[most_stable_index]]
clusters <- dfs$cluster

# Create the data frame using the most stable run
tsne_df <- data.frame(
  tSNE1 = most_stable_tsne$Y[,1],
  tSNE2 = most_stable_tsne$Y[,2],
  cluster = as.factor(clusters)
)

dfs$X <- rownames(dfs)
# Continue with your original code
tsne_df$X <- dfs$X


dfviz <- inner_join(dfs, tsne_df, by=c("X", "cluster"))



# Calculate cluster centroids for annotation placement
cluster_centers <- dfviz %>%
  group_by(cluster) %>%
  summarise(
    center_x = mean(tSNE1),
    center_y = mean(tSNE2),
    .groups = 'drop'
  ) %>%
  mutate(label = paste("Cluster", cluster))


a_colors <- c("#2941cc", "#eab30f", "#3da510")
caption_text <- paste("1) Consensus t-SNE visualization derived from 100 independent runs, showing patient clusters obtained through Partitioning Around Medoids (PAM) clustering based on correlation-weighted Gower distance. \n2) Bootstrapped cluster stability analysis indicated highest stability at k=3 clusters.
3) The clustering incorporated 21 baseline variables: LV injury, mRS prior to admission, age, sex, race, insurance type, history of hypertension, history of hyperlipidemia,\n history of tobacco use, history of diabetes mellitus, history of coronary artery disease, history of stroke, history of cancer, severe LAE, patent foramen ovale,\n left ventricular ejection fraction, WMA, multifocal diffusion-weighted imaging lesions, cortical diffusion-weighted imaging lesions, and baseline NIHSS")


gc0 <- ggplot(dfviz, aes(x = tSNE1, y = tSNE2, color = cluster)) +
  
  geom_point(alpha = 0.6, size=3) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21)),
    shape = guide_legend(override.aes = list(size = 5)),
    size = guide_legend()
  ) +
  
  labs(title = "Patient Clusters Visualization (Consensus t-SNE)",
       subtitle = "Baseline Characteristics among N=1869 patients",
       #caption = caption_text,
       x = "ct-SNE Dimension 1",
       y = "ct-SNE Dimension 2") +
  theme_minimal() +
  # Apply custom colors
  scale_color_manual(values = a_colors) +
  scale_fill_manual(values = a_colors) +
  theme(
    text = element_text(size = 12),
    plot.caption = element_text(hjust = 0, size=8),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

ggsave(gc0, filename = "figures/f2.pdf", width=10, height=7)





variables <- unique(c(
  "lv_injury",
  "mrs_pta",             
  "age",                  
  "sex",                  
  "race",                 
  "hispanic",            
  "insurance",
  "hx_htn",              
  "hx_hld",              
  "hx_tob",              
  "hx_dm",               
  "hx_cad",              
  "hx_stroke",           
  "hx_cancer",           
  "lae_severe",       
  
  "lvef",  
  "wma",
  "pfo",
  "mri_dwi_multifocal",  
  "dwi_cortical",        
  "tx_ap_no_ac",         
  "bnihss",              
  "time_to_event", "stroke_bleed_death2", "stroke_bleed_death2",
  "cluster"
))




dff <- select(dfs, variables)


# CASPR Causal Survival Forest Analysis
# Adapted from Sverdrup & Wager (2024) methodology

# Select relevant columns
dff <- na.omit(select(dff, all_of(variables)))

# ============================================================================
# CAUSAL SURVIVAL FOREST ANALYSIS
# ============================================================================

# =========================================
dff$tx <- ifelse(dff$tx_ap_no_ac ==1, 0, 1)


Y <- dff$time_to_event # Time to event
W <- dff$tx             # Treatment (1 = anticoagulation, 0 = antiplatelet, no anticoagulation)

# Create event indicator - you may need to modify this based on your outcome
# This assumes you have an event indicator variable; if not, you'll need to create one
D <- dff$stroke_bleed_death2


# Create formula programmatically
ps_formula <- reformulate(termlabels = clustering_vars, response = "tx")

# Use in weightit
ps.tx <- weightit(ps_formula,
                  data = as.data.frame(dff),
                  estimand = "ATE")

# Create covariates matrix (excluding treatment and outcome variables)
covariate_vars <- clustering_vars

dff$ps <- ps.tx$ps

X <- dff %>% 
  select(all_of(covariate_vars)) %>%
  # Handle missing values and ensure numeric format
  mutate(across(everything(), ~as.numeric(.))) %>%
  as.matrix()


# Filter Y, W, D to match X dimensions
valid_rows <- which(!is.na(Y) & !is.na(W) & !is.na(D))
Y <- Y[valid_rows]
W <- W[valid_rows]
D <- D[valid_rows]

# Check data structure
cat("Sample size:", length(Y), "\n")
cat("Treatment distribution:", table(W), "\n")
cat("Event distribution:", table(D), "\n")
cat("Covariates dimensions:", dim(X), "\n")




#####

# ============================================================================
# PROPENSITY SCORE WEIGHTED REGRESSION WITH INTERACTION
# ============================================================================

cat("\n=== PS-WEIGHTED REGRESSION WITH TX × CLUSTER INTERACTION ===\n\n")

# Extract propensity scores from the weightit object
dff$ps <- ps.tx$ps


# Calculate IPW weights (for ATE)
dff$ipw <- ifelse(dff$tx == 1, 1/dff$ps, 1/(1-dff$ps))

# Optional: Stabilize weights to reduce variance
dff$ipw_stabilized <- ifelse(dff$tx == 1, 
                             mean(dff$tx) / dff$ps,
                             (1 - mean(dff$tx)) / (1 - dff$ps))

# Check weight distribution
cat("IPW Weight Distribution:\n")
cat("  Mean:", round(mean(dff$ipw), 2), "\n")
cat("  SD:", round(sd(dff$ipw), 2), "\n")
cat("  Min:", round(min(dff$ipw), 2), "\n")
cat("  Max:", round(max(dff$ipw), 2), "\n")
cat("  Stabilized Mean:", round(mean(dff$ipw_stabilized), 2), "\n\n")

# Trim extreme weights if needed (optional)
weight_threshold <- quantile(dff$ipw, 0.99)
dff$ipw_trimmed <- pmin(dff$ipw, weight_threshold)

# Fit Cox model with interaction using IPW weights
cox_interaction <- coxph(Surv(time_to_event, stroke_bleed_death2) ~ 
                           tx * cluster,
                         data = dff,
                         weights = ipw)

cat("Cox Model with Treatment × Cluster Interaction (PS-weighted):\n")
print(summary(cox_interaction))










# ============================================================================
# FIT CAUSAL SURVIVAL FOREST
# ============================================================================

# Since this may not be from a randomized trial, we estimate propensity scores
# If it IS from an RCT, use W.hat = mean(W) instead



cs.forest <- causal_survival_forest(
  X, Y, W, D,
  horizon = 800,     # increase from default (5)
  honesty = TRUE,    # ensure honest estimation
  num.trees = 1000,       # more trees can help
  mtry = min(ceiling(sqrt(ncol(X)) + 10), ncol(X)),
  
  seed = 427
)
####### TREE DIAGNOSTICS
# Check actual leaf sizes in individual trees
check_tree_leaf_sizes <- function(forest, X, n_trees = 50) {
  all_leaf_sizes <- c()
  
  for (i in 1:min(n_trees, forest$`_num_trees`)) {
    tree <- get_tree(forest, i)
    leaves <- get_leaf_node(tree, X)
    leaf_sizes <- table(leaves)
    all_leaf_sizes <- c(all_leaf_sizes, as.numeric(leaf_sizes))
  }
  
  cat("Actual leaf sizes across", n_trees, "trees:\n")
  cat("  Min:", min(all_leaf_sizes), "\n")
  cat("  Mean:", round(mean(all_leaf_sizes), 1), "\n")
  cat("  Median:", median(all_leaf_sizes), "\n")
  cat("  Max:", max(all_leaf_sizes), "\n")
  cat("  % leaves >= min.node.size:", 
      round(100 * mean(all_leaf_sizes >= 45), 1), "%\n")
}

check_tree_leaf_sizes(cs.forest, X)


check_tree_structure <- function(forest, X, n_trees = 50) {
  for (i in 1:min(n_trees, 10)) {  # Just look at first 10
    tree <- get_tree(forest, i)
    leaves <- get_leaf_node(tree, X)
    n_leaves <- length(unique(leaves))
    cat("Tree", i, ": ", n_leaves, "leaves\n")
  }
}

check_tree_structure(cs.forest, X)

######

# Get treatment effect estimates
tau_hat <- predict(cs.forest)$predictions

# Distribution of estimated effects
summary(tau_hat)
hist(tau_hat, breaks = 30, main = "CATE Estimates")


#######
ate_result <- average_treatment_effect(cs.forest)
cat("Average Treatment Effect (ATE):\n")
print(ate_result)

# Extract the 95% CI explicitly
ate_estimate <- ate_result[1]
ate_se <- ate_result[2]

# Calculate 95% CI manually if needed
ate_lower <- ate_estimate - 1.96 * ate_se
ate_upper <- ate_estimate + 1.96 * ate_se

cat("\nDetailed ATE Results:\n")
cat("Estimate:", round(ate_estimate, 3), "\n")
cat("Standard Error:", round(ate_se, 3), "\n") 
cat("95% CI: (", round(ate_lower, 3), ", ", round(ate_upper, 3), ")\n")
cat("p-value:", round(2 * (1 - pnorm(abs(ate_estimate/ate_se))), 3), "\n")
# ============================================================================
# ANALYZE TREATMENT HETEROGENEITY
# ============================================================================

# Get predictions with uncertainty
pred <- predict(cs.forest, 
                estimate.variance = TRUE,
                num.threads = 1)

tau.hat <- pred$predictions
tau.se <- sqrt(pred$variance.estimates)
# Create confidence intervals
alpha <- 0.05  # for 95% CI
z_score <- qnorm(1 - alpha/2)
tau.lower <- tau.hat - z_score * tau.se
tau.upper <- tau.hat + z_score * tau.se

# Summary
results <- data.frame(
  tau_hat = tau.hat,
  tau_se = tau.se,
  tau_lower = tau.lower,
  tau_upper = tau.upper
)

# Save individual treatment effects for further analysis
dff$cate_estimate <- tau.hat
dff$cate_se <- tau.se
dff$cate_lower <- tau.lower
dff$cate_upper <- tau.upper



######
# Visualizing CATE

dff$idx <- as.numeric(rownames(dff))
dff$cluster <- factor(dff$cluster , levels=c("1", "2", "3"))
a_colors <- c("gray")

gcf <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=0.5, alpha=0.9, color="gray",
    position=position_dodge(1.5)) +
  geom_point(size=2, 
             position=position_dodge(1.5),
             alpha=0.6, shape=21, color="black", fill="gray") +
  geom_hline(yintercept = 0, linetype="dashed", color="darkred") +
  
  labs(
    title="Predicted Conditional Average Treatment Effects (CATE, 95%CI)",
    caption="1) Predictions obtained via a Causal Survival Forest model with propensity score adjustment. \n2) Lower CATE = Lower Restricted Mean Survival Time (RMST). \n3) Doubly-robust ATE estimate.") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -1, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )


print(gcf)

pdf("figures/f4.pdf", width=11, height=8)
print(gcf)
dev.off()
#####

# Analyzing variable importance (before any cluster related evaluation)






#########


library(fastshap)

pfun <- function(object, newdata) {
  predict(object, newdata)$predictions
}

# Explain subset of patients for computational efficiency
set.seed(427)
explain_idx <- sample(1:nrow(X), 1000)

shap_vals <- explain(
  cs.forest,
  X = as.data.frame(X),
  pred_wrapper = pfun,
  nsim = 50,
  newdata = as.data.frame(X[explain_idx, ])
)


# 1. SHAP Summary Plot (most important)
shap_importance <- data.frame(
  variable = colnames(shap_vals),
  mean_abs_shap = colMeans(abs(shap_vals))
) %>%
  arrange(desc(mean_abs_shap))

ggplot(shap_importance, aes(x = reorder(variable, mean_abs_shap), 
                            y = mean_abs_shap)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance (Mean |SHAP|)",
       x = "Variable",
       y = "Mean Absolute SHAP Value") +
  theme_minimal()

# 2. SHAP Beeswarm Plot (shows direction + distribution)

# Define binary variables that should be treated as categorical
binary_vars <- c("lv_injury", "hx_cancer", "hx_tob", "hx_htn", "hx_stroke", 
                 "hx_hld", "pfo", "hx_dm", "dwi_cortical", "wma", "hx_cad",
                 "mri_dwi_multifocal", "lae_severe", "hispanic")

# True numeric (continuous) variables only
num_vars_plot <- c("lvef", "age", "bnihss", "mrs_pta")

# Categorical = original factors + binary variables
cat_vars_plot <- c("race", "insurance", "sex", binary_vars)

# Keep only variables that exist in your data
num_vars_plot <- num_vars_plot[num_vars_plot %in% colnames(shap_vals)]
cat_vars_plot <- cat_vars_plot[cat_vars_plot %in% colnames(shap_vals)]

var_labels <- c(
  "lvef"              = "LVEF (%)",
  "age"               = "Age",
  "bnihss"            = "Baseline NIHSS",
  "mrs_pta"           = "Pre-stroke mRS",
  "lv_injury"         = "LV Injury",
  "hx_cancer"         = "History of Cancer",
  "hx_tob"            = "Tobacco Use",
  "hx_htn"            = "Hypertension",
  "hx_stroke"         = "Prior Stroke",
  "hx_hld"            = "Hyperlipidemia",
  "pfo"               = "PFO",
  "hx_dm"             = "Diabetes",
  "dwi_cortical"      = "Cortical DWI Lesion",
  "wma"               = "Wall Motion Abnormality",
  "hx_cad"            = "Coronary Artery Disease",
  "mri_dwi_multifocal"= "Multifocal DWI",
  "lae_severe"        = "Severe LA Enlargement",
  "hispanic"          = "Hispanic Ethnicity",
  "race"              = "Race",
  "insurance"         = "Insurance",
  "sex"               = "Sex"
)
# For numeric
shap_long_num <- shap_vals %>%
  as.data.frame() %>%
  select(all_of(num_vars_plot)) %>%
  mutate(patient_id = as.integer(row_number())) %>%
  pivot_longer(-patient_id, names_to = "variable", values_to = "shap_value") %>%
  left_join(
    dff[explain_idx, num_vars_plot] %>%
      mutate(patient_id = as.integer(row_number())) %>%
      pivot_longer(-patient_id, names_to = "variable", values_to = "feature_value"),
    by = c("patient_id", "variable")
  )

# For categorical
shap_long_cat <- shap_vals %>%
  as.data.frame() %>%
  select(all_of(cat_vars_plot)) %>%
  mutate(patient_id = as.integer(row_number())) %>%
  pivot_longer(-patient_id, names_to = "variable", values_to = "shap_value") %>%
  left_join(
    dff[explain_idx, cat_vars_plot] %>%
      mutate(patient_id = as.integer(row_number()),
             across(all_of(cat_vars_plot), ~as.character(.))) %>%
      pivot_longer(-patient_id, names_to = "variable", values_to = "feature_category"),
    by = c("patient_id", "variable")
  )

# Top variables for each plot
top_num_vars <- shap_importance %>%
  filter(variable %in% num_vars_plot) %>%
  slice_head(n = 10) %>%
  pull(variable)

top_cat_vars <- shap_importance %>%
  filter(variable %in% cat_vars_plot) %>%
  slice_head(n = 15) %>%
  pull(variable)

library(patchwork)

# Plot numeric (continuous) variables with color gradient
p1 <- shap_long_num %>%
  filter(variable %in% top_num_vars) %>%
  ggplot(aes(x = shap_value, 
             y = reorder(variable, abs(shap_value), mean),
             color = feature_value)) +
  geom_jitter(alpha = 0.5, height = 0.2, width = 0) +
  scale_color_viridis_c(option = "plasma") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  
  labs(title = "SHAP Values: Continuous Variables",
       x = "SHAP value (impact on CATE)",
       y = NULL,
       color = "Feature\nValue") +
  theme_minimal() +
  scale_y_discrete(labels = function(x) ifelse(!is.na(var_labels[x]), var_labels[x], x))

# Adjust annotation positions based on actual data range
p2 <- shap_long_cat %>%
  filter(variable %in% top_cat_vars) %>%
  ggplot(aes(x = shap_value, 
             y = reorder(variable, abs(shap_value), mean),
             color = feature_category)) +
  geom_jitter(alpha = 0.6, height = 0.2, width = 0) +
  scale_color_brewer(palette = "Set2") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  
  labs(title = "SHAP Values: Categorical & Binary Variables",
       x = "SHAP value (impact on CATE)",
       y = NULL,
       color = "Category") +
  theme_minimal() +
  scale_y_discrete(labels = function(x) ifelse(!is.na(var_labels[x]), var_labels[x], x))

# Display both
p1 / p2

# Save
pdf("shap_values.pdf", width = 12, height = 10)
print(p1 / p2)
dev.off()
# 3. Dependence plots for key variables
plot_shap_dependence <- function(var_name1, var_name2 = NULL) {
  
  if (is.null(var_name2)) {
    # Single variable dependence plot
    shap_long_num %>%
      filter(variable == var_name1) %>%
      ggplot(aes(x = feature_value, y = shap_value)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_smooth(method = "loess", color = "red") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = paste("SHAP Dependence:", var_name1),
           x = var_name1,
           y = "SHAP value (impact on CATE)") +
      theme_minimal()
    
  } else {
    # Two variable interaction plot
    # Get SHAP values for var_name1, colored by var_name2
    shap_var1 <- shap_long_num %>%
      filter(variable == var_name1) %>%
      select(patient_id, shap_value, feature_value)
    
    shap_var2 <- shap_long_num %>%
      filter(variable == var_name2) %>%
      select(patient_id, color_value = feature_value)
    
    shap_var1 %>%
      left_join(shap_var2, by = "patient_id") %>%
      ggplot(aes(x = feature_value, y = shap_value, color = color_value)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", color = "red", se = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_viridis_c(option = "plasma") +
      labs(title = paste("SHAP Dependence:", var_name1, "by", var_name2),
           x = var_name1,
           y = "SHAP value (impact on CATE)",
           color = var_name2) +
      theme_minimal()
  }
}

# Usage examples:
plot_shap_dependence("lvef")
plot_shap_dependence("lvef", "age")
plot_shap_dependence("age", "lvef")

# Plot for top variables

plot_shap_dependence("age", "lv_injury")


##### Now inspect specific clinically interesting CATE-defining covs
# e.g. LVEF, LV



a_colors <- c("#2941cc",  "#eab30f",  "#3da510", "purple", "red")

gcf2 <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate, 
                        color=cluster,
                        fill=cluster)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=0.5, alpha=0.9,
    position=position_dodge(1.5)) +
  geom_point(size=2, 
             position=position_dodge(1.5),
             alpha=0.6, shape=21, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="darkred") +
  
  # Apply custom colors
  scale_color_manual(values = a_colors) +
  scale_fill_manual(values = a_colors) +
  
  labs(
    title="") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -1, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )


print(gcf2)

pdf("figures/f4.pdf", width=11, height=8)
print(gcf2)
dev.off()




a_colors <- c("#2941cc",  "#eab30f","purple", "red")
dff$lv_injury <- as.factor(dff$lv_injury )

gcfLV <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate, 
                         color=lv_injury,
                         fill=lv_injury)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=0.5, alpha=0.9,
    position=position_dodge(1.5)) +
  geom_point(size=2, 
             position=position_dodge(1.5),
             alpha=0.6, shape=21, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="darkred") +
  
  # Apply custom colors
  scale_color_manual(values = a_colors) +
  scale_fill_manual(values = a_colors) +
  
  labs(
    title="Predicted Conditional Average Treatment Effects (CATE) by LV Injury Status",
    caption="1) Predictions obtained via a Causal Survival Forest model with propensity score adjustment. \n2) Lower CATE = Lower Restricted Mean Survival Time (RMST). \n3) Doubly-robust ATE estimate.") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -1, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  annotate("text", x = 1 * length(unique(dff$idx)), y = 1, 
           label = "> 0 Favors Anticoagulation ", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "darkred", 
           fontface = "bold") +
  annotate("text", x = 1 * length(unique(dff$idx)), y = -10, 
           label = "< 0 Favors Antiplatelet", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "darkred", 
           fontface = "bold") +
  
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )

print(gcfLV)

pdf("2_cate_estims_ordered-lvinj-new.pdf", width=11, height=8)
print(gcfLV)
dev.off()

dff$lv_injury <- as.numeric(dff$lv_injury)

dff$bnihss <- as.numeric(dff$bnihss)

gcfAA <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate, 
                         color=bnihss,
                         fill=bnihss)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=0.5,
    position=position_dodge(1.2)) +
  geom_point(size=2, 
             position=position_dodge(1.2),
             alpha=0.8, shape=21, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="white") +
  
  labs(
    title="Predicted Conditional Average Treatment Effects (CATE) by NIHSS",
    caption="1) Predictions obtained via a Causal Survival Forest model with propensity score adjustment. \n2) Lower CATE = Lower Restricted Mean Survival Time (RMST). \n3) Doubly-robust ATE estimate.") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -5, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  annotate("text", x = 0.07 * length(unique(dfc$idx)), y = 3, 
           label = "> 0 Favors Anticoagulation ", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "white", 
           fontface = "bold") +
  annotate("text", x = 0.07 * length(unique(dfc$idx)), y = -3, 
           label = "< 0 Favors Antiplatelet", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "white", 
           fontface = "bold") +
  #ylim(-25, 20) +
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )


print(gcfAA)

pdf("3_cate_estims_ordered-bnihss-new.pdf", width=11, height=8)
print(gcfAA)
dev.off()

############

gcfAA <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate, 
                         color=age,
                         fill=age)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=0.9,
    position=position_dodge(1.2)) +
  geom_point(size=2, 
             position=position_dodge(1.2),
             alpha=0.8, shape=21, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="white") +
  
  labs(
    title="Predicted Conditional Average Treatment Effects (CATE) by Age",
    caption="1) Predictions obtained via a Causal Survival Forest model with propensity score adjustment. \n2) Lower CATE = Lower Restricted Mean Survival Time (RMST). \n3) Doubly-robust ATE estimate.") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -5, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  annotate("text", x = 0.07 * length(unique(dfc$idx)), y = 3, 
           label = "> 0 Favors Anticoagulation ", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "white", 
           fontface = "bold") +
  annotate("text", x = 0.07 * length(unique(dfc$idx)), y = -3, 
           label = "< 0 Favors Antiplatelet", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "white", 
           fontface = "bold") +
  #ylim(-25, 20) +
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )


print(gcfAA)
pdf("4_cate_estims_ordered-age-new.pdf", width=11, height=8)
print(gcfAA)
dev.off()


##########
dff$mri_fazekas <- as.numeric(dff$mri_fazekas)
a_colors <- c("#2941cc",   "#eab30f","purple", "red", "#3da510")

gcfFaz <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate, 
                          color=mri_fazekas,
                          fill=mri_fazekas)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=1.2, alpha=0.9,
    position=position_dodge(1.5)) +
  geom_point(size=2, 
             position=position_dodge(1.5),
             alpha=0.6, shape=21, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="darkred") +
  
  labs(
    title="Predicted Conditional Average Treatment Effects (CATE) by Fazekas MRI Score",
    caption="1) Predictions obtained via a Causal Survival Forest model with propensity score adjustment. \n2) Lower CATE = Lower Restricted Mean Survival Time (RMST). \n3) Doubly-robust ATE estimate.") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -1, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  annotate("text", x = 1 * length(unique(dff$idx)), y = 1, 
           label = "> 0 Favors Anticoagulation ", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "darkred", 
           fontface = "bold") +
  annotate("text", x = 1 * length(unique(dff$idx)), y = -10, 
           label = "< 0 Favors Antiplatelet", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "darkred", 
           fontface = "bold") +
  
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )

print(gcfFaz)

pdf("5_cate_estims_ordered-fazekas-new.pdf", width=11, height=8)
print(gcfFaz)
dev.off()


gcfLVEF <- ggplot(dff, aes(x=reorder(idx, cate_estimate), y=cate_estimate, 
                           color=lvef,
                           fill=lvef)) +
  geom_errorbar(
    aes(ymin=cate_lower, ymax=cate_upper),
    width=1.2, alpha=0.9,
    position=position_dodge(1.5)) +
  geom_point(size=2, 
             position=position_dodge(1.5),
             alpha=0.6, shape=21, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="darkred") +
  
  labs(
    title="Predicted Conditional Average Treatment Effects (CATE) by LVEF",
    caption="1) Predictions obtained via a Causal Survival Forest model with propensity score adjustment. \n2) Lower CATE = Lower Restricted Mean Survival Time (RMST). \n3) Doubly-robust ATE estimate.") +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Overall ATE (95%CI): ", round(ate_estimate, 2),
                          " (", round(ate_lower, 2), ", ", 
                          round(ate_upper, 2), ")"),
           hjust = 1.1, vjust = -1, 
           size = 4, color = "gray30",
  ) +
  
  ylab("CATE (95% CI)") +
  xlab("Observation") +
  annotate("text", x = 1 * length(unique(dff$idx)), y = 1, 
           label = "> 0 Favors Anticoagulation ", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "darkred", 
           fontface = "bold") +
  annotate("text", x = 1 * length(unique(dff$idx)), y = -10, 
           label = "< 0 Favors Antiplatelet", 
           hjust = 0, vjust = 0.5,
           size = 4, color = "darkred", 
           fontface = "bold") +
  
  
  theme_tuf2 +
  theme(
    text = element_text( size = 12),
    axis.text = element_text(),
    axis.title = element_text(),
    plot.title = element_text( hjust = 0.5, size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Additional academic styling
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50", size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  )


pdf("6_cate_estims_ordered-lvef-new.pdf", width=11, height=8)
print(gcfLVEF)
dev.off()




dfsumm <- select(dfs, c("age_orig", "sex", "race",
                        "hispanic", "insurance",
                        colnames(dfs)[2:3], 
                        colnames(dfs)[9:22],
                        "time_to_event", "cluster"))
# Cluster-based summary table for stroke data
cluster_summary <- dfsumm %>%
  tbl_summary(
    type = list(
      mrs_pta ~ "categorical",
      age_orig ~ "continuous",
      sex ~ "categorical", 
      lv_injury ~ "continuous",
      race ~ "categorical",
      insurance ~ "categorical",
      bnihss ~ "categorical",
      time_to_event ~ "continuous"
      
    ),
    by = cluster,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      lv_injury ~ "LV Injury",
      mrs_pta ~ "mRS Prior to Admission",
      age_orig ~ "Age",
      sex ~ "Sex",
      race ~ "Race",
      hispanic ~ "Hispanic Ethnicity",
      insurance ~ "Insurance Type",
      
      hx_htn ~ "History of Hypertension",
      hx_hld ~ "History of Hyperlipidemia",
      hx_tob ~ "History of Tobacco Use",
      hx_dm ~ "History of Diabetes",
      hx_cad ~ "History of CAD",
      hx_stroke ~ "History of Stroke",
      hx_cancer ~ "History of Cancer",
      lae_severe ~ "Severe Left Atrial Enlargement",
      pfo ~ "Patent Foramen Ovale",
      lvef ~ "Left Ventricular Ejection Fraction",
      wma ~ "Wall Motion Abnormality", 
      mri_dwi_multifocal ~ "Multifocal DWI Lesions",
      dwi_cortical ~ "Cortical DWI Lesions",
      bnihss ~ "Baseline NIHSS",
      time_to_event ~ "Time to Event (days)"
      
    )
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  bold_labels() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Cluster Group**")

print(cluster_summary)

# Export to Word document
ft_cluster_summary <- as_flex_table(cluster_summary)
flextable::save_as_docx(ft_cluster_summary, 
                        path = "Cluster_Characteristics_ISC.docx")



########## KAPLAN-MEIER


# Create survival object

surv_obj <- Surv(time = dff$time_to_event, 
                 event = dff$stroke_bleed_death2)



# Fit by both cluster and treatment
km_fit <- survfit(surv_obj ~ cluster, data = dff)

km_plot <- ggsurvplot(
  km_fit,
  data = dff,
  
  pval = FALSE,
  conf.int = FALSE,
  
  # Truncate at 3 years
  xlim = c(0, 1095),
  break.time.by = 365,
  
  # Risk table
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,
  risk.table.fontsize = 3.5,
  risk.table.title = "Number at Risk",
  
  palette = c("#2941cc", "#eab30f", "#3da510"),
  linetype = "solid",
  size = 1.2,
  
  title = "Kaplan-Meier Survival Curves by Patient Cluster",
  subtitle = "Time to Event",
  xlab = "Time (days)",
  ylab = "Event-free Survival Probability",
  legend.title = "Cluster",
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3")
)



# Save as PDF
pdf("kaplan_meier_by_cluster_isc.pdf", width = 10, height = 8)
print(km_plot)
dev.off()



# Or for better labels:
dff$group <- paste0("Cluster ", dff$cluster, " - ", 
                    ifelse(dff$tx_ap_no_ac == 1, "AP", "AC"))

# Fit by combined variable
km_fit <- survfit(surv_obj ~ group, data = dff)

km_plot <- ggsurvplot(
  km_fit,
  data = dff,
  
  # Basic settings
  pval = FALSE,
  conf.int = FALSE,
  
  # Styling - 6 colors for 6 groups (3 clusters × 2 treatments)
  palette = c("#19297c", "#5c6fd4",  # Cluster 1: AP (dark blue), AC (light blue)
              "#efb30f", "#fff2cc",  # Cluster 2: AP (dark gold), AC (light gold)
              "#2c6115", "#7bc65a"), # Cluster 3: AP (dark green), AC (light green)
  size = 1.2,
  
  # Labels
  title = "Kaplan-Meier Survival Curves by Patient Cluster and Treatment",
  subtitle = "Time to Event",
  xlab = "Time (days)",
  ylab = "Event-free Survival Probability",
  legend.title = "Group",
  
  break.time.by = 360
)

# Save as PDF
pdf("kaplan_meier_by_cluster_and_tx_isc.pdf", width = 12, height = 8)
print(km_plot)
dev.off()

# Fit by cluster only
km_fit <- survfit(surv_obj ~ cluster, data = dff)

# Create survival object
surv_obj <- Surv(time = dff$time_to_event, 
                 event = dff$stroke_bleed_death2)

# Create better labels for faceting
dff$Treatment <- ifelse(dff$tx_ap_no_ac == 1, 
                        "AP", 
                        "AC")

# Fit by cluster only
km_fit <- survfit(surv_obj ~ cluster, data = dff)

km_plot <- ggsurvplot_facet(
  km_fit,
  data = dff,
  facet.by = "Treatment",
  
  # Styling
  palette = c("#2941cc", "#eab30f", "#3da510"),
  linetype = "solid",
  size = 1.2,
  
  # Labels
  xlab = "Time (days)",
  ylab = "Event-free Survival Probability",
  legend.title = "Cluster",
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  
  break.time.by = 360,
  
  # Panel labels
  panel.labs = list(Treatment = c("AP", "AC"))
)

# Save as PDF
pdf("kaplan_meier_by_cluster_faceted_isc.pdf", width = 10, height = 6)
print(km_plot)
dev.off()
########

# Now test for HTE for clinically relevant variables:

# ============================================================================
# BEST LINEAR PROJECTION - Identify which covariates predict heterogeneity
# ============================================================================

# ============================================================================
# BEST LINEAR PROJECTION - HETEROGENEITY TESTING TABLE
# ============================================================================

library(flextable)

# Define clinically hypothesized variables
hte_vars <- c("lvef", "lv_injury", "age", "hx_cancer", "bnihss", "hx_htn", 
              "hx_cad", "wma", "pfo", "cluster")

# Function to extract BLP results for a single variable
extract_blp_results <- function(forest, data, var_name) {
  
  tryCatch({
    key_cov <- data[var_name]
    blp <- best_linear_projection(forest, key_cov)
    
    # BLP is a matrix - extract directly
    coef_matrix <- as.matrix(blp)
    var_names <- rownames(coef_matrix)
    
    # Get non-intercept rows
    non_intercept <- var_names != "(Intercept)"
    
    if (sum(non_intercept) > 0) {
      results <- data.frame(
        variable = var_names[non_intercept],
        estimate = coef_matrix[non_intercept, "Estimate"],
        std_error = coef_matrix[non_intercept, "Std. Error"],
        p_value = coef_matrix[non_intercept, "Pr(>|t|)"],
        predictor = var_name,
        row.names = NULL
      )
      return(results)
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    message(paste("Error with variable:", var_name, "-", e$message))
    return(NULL)
  })
}

# Run BLP for each variable
blp_results_list <- lapply(hte_vars, function(v) {
  cat("Processing:", v, "\n")
  extract_blp_results(cs.forest, dff, v)
})

# Combine results
blp_table <- do.call(rbind, blp_results_list[!sapply(blp_results_list, is.null)])

# Format table
blp_table <- blp_table %>%
  mutate(
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error,
    estimate_ci = sprintf("%.2f (%.2f, %.2f)", estimate, ci_lower, ci_upper),
    p_formatted = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

# Create publication table
pub_table <- blp_table %>%
  select(Variable = variable, 
         `Estimate (95% CI)` = estimate_ci,
         `P-value` = p_formatted)

print(pub_table)

# Export to Word
ft_blp <- flextable(pub_table) %>%
  add_header_lines("Table X. Best Linear Projection Analysis of Treatment Effect Heterogeneity") %>%
  add_footer_lines("Positive coefficients indicate greater benefit from anticoagulation.") %>%
  autofit() %>%
  theme_booktabs()

save_as_docx(ft_blp, path = "Table_BLP_HTE_Analysis.docx")





#######################

#1. Who are the people with positive effects?

# Create CATE-based groups
dff$cate_group <- cut(dff$cate_estimate, 
                      breaks = quantile(dff$cate_estimate, c(0, 0.5, 1)),
                      labels = c("AP relative benefit", "AC relative benefit"),
                      include.lowest = TRUE)

table(dff$cate_group)

# Compare characteristics

dfc_raw <- read_dta("data/CASPR_lvef_092424_FINAL.dta")
dfc_raw <- dplyr::filter(dfc_raw, include == 1)

mrs_min <- min(dfc_raw$mrs_pta, na.rm = TRUE)
mrs_max <- max(dfc_raw$mrs_pta, na.rm = TRUE)
nihss_min <- min(dfc_raw$bnihss, na.rm = TRUE)
nihss_max <- max(dfc_raw$bnihss, na.rm = TRUE)
age_min <- min(dfc_raw$age, na.rm = TRUE)
age_max <- max(dfc_raw$age, na.rm = TRUE)

# Reverse the min-max scaling
dff$mrs_pta_orig <- round(dff$mrs_pta * (mrs_max - mrs_min) + mrs_min)
dff$bnihss_orig <- round(dff$bnihss * (nihss_max - nihss_min) + nihss_min)
dff$age_orig <- round(dff$age * (age_max - age_min) + age_min)


dff$mrs_pta_orig <- factor(dff$mrs_pta_orig)


# ---- Prior antithrombotic therapy classification (Table 1) ----
# dfc_raw is the filtered raw data (include == 1); its row order matches dff$idx,
# so the classification can be derived here and mapped back onto the analysis sample.
ac_agents  <- c("pta_vka", "pta_hep_lmwh", "pta_apix", "pta_dabi",
                "pta_riva", "pta_edox", "pta_argatroban")
ap_agents  <- c("pta_asa", "pta_clopi_prasu", "pta_ticag")
alt_ap_ids <- c(99, 205, 371, 614, 1790, 1847, 1921)

# Coerce agent flags to 0/1 and flag any-positive per row
ac_mat <- sapply(dfc_raw[ac_agents], function(x) as.numeric(x) == 1)
ap_mat <- sapply(dfc_raw[ap_agents], function(x) as.numeric(x) == 1)

# Anticoagulation: any AC agent OR pta_ac already coded 1
on_ac <- (rowSums(ac_mat, na.rm = TRUE) > 0) 

# Antiplatelet: any listed AP agent OR a record on an unlisted alternative AP agent
on_ap <- (rowSums(ap_mat, na.rm = TRUE) > 0) |
  (as.numeric(dfc_raw$record_deid) %in% alt_ap_ids)

dfc_raw$prior_antithrombotic <- factor(
  dplyr::case_when(
    on_ac          ~ "Prior anticoagulation",
    on_ap & !on_ac ~ "Antiplatelet without anticoagulation",
    TRUE           ~ "None"
  ),
  levels = c("None", "Antiplatelet without anticoagulation", "Prior anticoagulation")
)

# Map onto the analysis sample by row index
dff$prior_antithrombotic <- dfc_raw$prior_antithrombotic[dff$idx]

dfsumm0 <- select(dff, c(colnames(dff)[1:20], colnames(dff)[22],
                         "time_to_event","tx", "bnihss_orig", "mrs_pta_orig",
                         "age_orig", "prior_antithrombotic"))



# Cluster-based summary table for stroke data
sample_summary <- dfsumm0 %>%
  tbl_summary(
    type = list(
      mrs_pta_orig ~ "categorical",
      age_orig ~ "continuous",
      sex ~ "categorical", 
      lv_injury ~ "categorical",
      race ~ "categorical",
      insurance ~ "categorical",
      bnihss_orig ~ "continuous",
      time_to_event ~ "continuous",
      lvef ~ "continuous",
      prior_antithrombotic ~ "categorical"
      
    ),
    by = tx,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      lv_injury ~ "LV Injury",
      mrs_pta_orig ~ "mRS Prior to Admission",
      age_orig ~ "Age",
      sex ~ "Sex",
      race ~ "Race",
      hispanic ~ "Hispanic Ethnicity",
      insurance ~ "Insurance Type",
      
      hx_htn ~ "History of Hypertension",
      hx_hld ~ "History of Hyperlipidemia",
      hx_tob ~ "History of Tobacco Use",
      hx_dm ~ "History of Diabetes",
      hx_cad ~ "History of CAD",
      hx_stroke ~ "History of Stroke",
      hx_cancer ~ "History of Cancer",
      lae_severe ~ "Severe Left Atrial Enlargement",
      pfo ~ "Patent Foramen Ovale",
      lvef ~ "Left Ventricular Ejection Fraction",
      wma ~ "Wall Motion Abnormality", 
      mri_dwi_multifocal ~ "Multifocal DWI Lesions",
      dwi_cortical ~ "Cortical DWI Lesions",
      bnihss_orig ~ "Baseline NIHSS",
      time_to_event ~ "Time to Event (days)",
      prior_antithrombotic ~ "Prior Antithrombotic Therapy"
      
    )
  ) %>%
  
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  ),
  pvalue_fun = function(x) style_pvalue(x, digits = 2)
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Treatment Group**")

print(sample_summary)

# Export to Word document
ft_sample_summary <- as_flex_table(sample_summary)
flextable::save_as_docx(ft_sample_summary, 
                        path = "t1-rr.docx")


# Shared labels (define once; reuse for tbl_missing and tbl_pfo_ef)
t1_labels <- list(
  lv_injury ~ "LV Injury", mrs_pta ~ "mRS Prior to Admission", age ~ "Age",
  sex ~ "Sex", race ~ "Race", hispanic ~ "Hispanic Ethnicity",
  insurance ~ "Insurance Type", hx_htn ~ "History of Hypertension",
  hx_hld ~ "History of Hyperlipidemia", hx_tob ~ "History of Tobacco Use",
  hx_dm ~ "History of Diabetes", hx_cad ~ "History of CAD",
  hx_stroke ~ "History of Stroke", hx_cancer ~ "History of Cancer",
  lae_severe ~ "Severe Left Atrial Enlargement", pfo ~ "Patent Foramen Ovale",
  lvef ~ "Left Ventricular Ejection Fraction", wma ~ "Wall Motion Abnormality",
  mri_dwi_multifocal ~ "Multifocal DWI Lesions", dwi_cortical ~ "Cortical DWI Lesions",
  bnihss ~ "Baseline NIHSS"
)

# PFO or HFrEF (20% < EF < 30%) vs. neither
dff$pfo_lowef <- factor(
  +((dff$pfo %in% 1) | (dff$lvef > 20 & dff$lvef < 30 & !is.na(dff$lvef))),
  levels = c(0, 1), labels = c("Neither", "PFO or EF 20-30%")
)

tbl_pfo_ef <- dff %>%
  select(all_of(orig_vars), pfo_lowef) %>%
  tbl_summary(
    by = pfo_lowef, missing = "no",
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    label = t1_labels
  ) %>%
  add_overall() %>%
  add_p(test = list(all_continuous() ~ "kruskal.test",
                    all_categorical() ~ "chisq.test")) %>%
  bold_p(t = 0.05) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**PFO or Reduced EF (20-30%)**")

flextable::save_as_docx(as_flex_table(tbl_pfo_ef),
                        path = "t1_pfo_ef_comparison.docx")


ref <- haven::read_dta("data/CASPR_lvef_092424_FINAL.dta") %>% filter(include == 1)
inv <- function(s, v) { r <- range(ref[[v]], na.rm = TRUE); s * (r[2] - r[1]) + r[1] }

round(inv(c(.58, .46, .70), "age"))      # Age   — Neither
round(inv(c(.54, .37, .65), "age"))      # Age   — PFO/low-EF
round(inv(c(.14, .05, .30), "bnihss"))   # NIHSS — Neither
round(inv(c(.08, .03, .19), "bnihss"))   # NIHSS — PFO/low-EF
round(inv(c(.58, .51, .65), "lvef"), 1)  # EF    — both groups

# Individual components of the composite outcome
cat("Composite outcome (stroke_bleed_death2):", sum(dff$stroke_bleed_death2), 
    "(", round(100*mean(dff$stroke_bleed_death2), 1), "%)\n")

# You'll need to check your exact variable names, but likely:
# Recurrent stroke
cat("Recurrent ischemic stroke:", sum(dff$recurr_ischemic, na.rm=TRUE), 
    "(", round(100*mean(dff$recurr_ischemic, na.rm=TRUE), 1), "%)\n")

# Major bleeding  
cat("Major bleeding (ISTH):", sum(dff$isth_bleed, na.rm=TRUE),
    "(", round(100*mean(dff$isth_bleed, na.rm=TRUE), 1), "%)\n")

# Death
cat("All-cause mortality:", sum(dff$death, na.rm=TRUE),
    "(", round(100*mean(dff$death, na.rm=TRUE), 1), "%)\n")

##### EXTRA VIZ:


# Add cluster-level ATE from Cox model to the CATE plot
cluster_ates <- dff %>%
  group_by(cluster) %>%
  summarise(cluster_ate = mean(cate_estimate))  
# Or extract from your Cox model coefficients

gcf_comparison <- ggplot(dff, aes(x = cluster, y = cate_estimate, fill = cluster)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  scale_fill_manual(values = c("#2941cc", "#eab30f", "#3da510")) +
  labs(title = "Individual CATEs by Cluster",
       subtitle = "Boxplots show within-cluster distribution",
       y = "CATE (days)",
       x = "Cluster") +
  theme_minimal()

pdf("fS3.pdf")
print(gcf_comparison)
dev.off()



#############


# Bootstrap resampling for cluster stability
# ---------------------------------------------------------
set.seed(123)

library(fpc)
library(ggplot2)

# Function to run stability analysis for different k values
run_stability_analysis <- function(dist_matrix, k_values, B = 100) {
  stability_results <- list()
  
  for (k in k_values) {
    cat(sprintf("\nRunning bootstrap stability analysis for k = %d...\n", k))
    
    stab_results <- clusterboot(as.matrix(dist_matrix), 
                                B = B,  # Number of bootstrap replicates
                                distances = TRUE,  # Input is a distance matrix
                                bootmethod = "boot",  # Regular bootstrap
                                clustermethod = pamkCBI,  # PAM clustering
                                k = k,  # Number of clusters
                                count = TRUE)  # Count frequency of cluster recovery
    
    # Store results including bootstrap standard deviations
    stability_results[[paste0("k", k)]] <- list(
      k = k,
      bootmean = stab_results$bootmean,
      bootsd = apply(stab_results$bootresult, 2, sd),  # Calculate SDs from bootstrap results
      bootbrd = stab_results$bootbrd,
      bootrecover = stab_results$bootrecover,
      mean_jaccard = mean(stab_results$bootmean),
      min_jaccard = min(stab_results$bootmean)
    )
    
    # Print results for this k
    cat(sprintf("\nBootstrap Stability Results for k = %d:\n", k))
    cat("Mean Jaccard similarities per cluster:\n")
    print(round(stab_results$bootmean, 3))
    cat("\nDissolution counts per cluster:\n")
    print(stab_results$bootbrd)
    cat("\nRecovery counts per cluster:\n")
    print(stab_results$bootrecover)
    cat(sprintf("Overall mean Jaccard: %.3f\n", mean(stab_results$bootmean)))
    cat(sprintf("Minimum Jaccard: %.3f\n", min(stab_results$bootmean)))
  }
  
  return(stability_results)
}

# Run stability analysis for k = 2, 3, 4, 5, 6
k_values <- 2:10
stability_results <- run_stability_analysis(wgower_dist, k_values, B = 100)

# Enhanced summary function with confidence intervals
create_stability_summary_with_ci <- function(stability_results, confidence_level = 0.95) {
  summary_data <- data.frame(
    k = integer(),
    cluster = integer(),
    jaccard = numeric(),
    jaccard_sd = numeric(),
    jaccard_lower = numeric(),
    jaccard_upper = numeric(),
    dissolved = integer(),
    recovered = integer()
  )
  
  overall_summary <- data.frame(
    k = integer(),
    mean_jaccard = numeric(),
    mean_jaccard_sd = numeric(),
    mean_jaccard_lower = numeric(),
    mean_jaccard_upper = numeric(),
    min_jaccard = numeric(),
    n_stable = integer(),
    n_partly_stable = integer(),
    n_unstable = integer()
  )
  
  # Calculate z-score for confidence interval
  alpha <- 1 - confidence_level
  z_score <- qnorm(1 - alpha/2)
  
  for (result_name in names(stability_results)) {
    result <- stability_results[[result_name]]
    k <- result$k
    
    # Detailed data with confidence intervals
    for (i in seq_along(result$bootmean)) {
      jaccard_mean <- result$bootmean[i]
      jaccard_sd <- result$bootsd[i]
      
      # Calculate confidence intervals
      jaccard_lower <- jaccard_mean - z_score * jaccard_sd / sqrt(100)  # B = 100
      jaccard_upper <- jaccard_mean + z_score * jaccard_sd / sqrt(100)
      
      # Ensure bounds are within [0,1]
      jaccard_lower <- max(0, jaccard_lower)
      jaccard_upper <- min(1, jaccard_upper)
      
      summary_data <- rbind(summary_data, data.frame(
        k = k,
        cluster = i,
        jaccard = jaccard_mean,
        jaccard_sd = jaccard_sd,
        jaccard_lower = jaccard_lower,
        jaccard_upper = jaccard_upper,
        dissolved = result$bootbrd[i],
        recovered = result$bootrecover[i]
      ))
    }
    
    # Overall summary with confidence intervals
    jaccards <- result$bootmean
    jaccards_sd <- result$bootsd
    
    mean_jaccard <- result$mean_jaccard
    mean_jaccard_sd <- mean(jaccards_sd, na.rm = TRUE)
    mean_jaccard_lower <- mean_jaccard - z_score * mean_jaccard_sd / sqrt(100)
    mean_jaccard_upper <- mean_jaccard + z_score * mean_jaccard_sd / sqrt(100)
    
    # Ensure bounds are within [0,1]
    mean_jaccard_lower <- max(0, mean_jaccard_lower)
    mean_jaccard_upper <- min(1, mean_jaccard_upper)
    
    n_stable <- sum(jaccards > 0.75)
    n_partly_stable <- sum(jaccards >= 0.6 & jaccards <= 0.75)
    n_unstable <- sum(jaccards < 0.6)
    
    overall_summary <- rbind(overall_summary, data.frame(
      k = k,
      mean_jaccard = mean_jaccard,
      mean_jaccard_sd = mean_jaccard_sd,
      mean_jaccard_lower = mean_jaccard_lower,
      mean_jaccard_upper = mean_jaccard_upper,
      min_jaccard = result$min_jaccard,
      n_stable = n_stable,
      n_partly_stable = n_partly_stable,
      n_unstable = n_unstable
    ))
  }
  
  return(list(detailed = summary_data, overall = overall_summary))
}

# Create summaries with confidence intervals
summaries_ci <- create_stability_summary_with_ci(stability_results, confidence_level = 0.95)
detailed_summary_ci <- summaries_ci$detailed
overall_summary_ci <- summaries_ci$overall

# Print overall summary with confidence intervals
print(overall_summary_ci)

# Enhanced Plot: Individual cluster stability with error bars
plot1_ci <- ggplot(detailed_summary_ci, aes(x = factor(cluster), y = jaccard)) +
  geom_col(color = "black", fill = "gray", alpha = 0.7) +
  geom_errorbar(aes(ymin = jaccard_lower, ymax = jaccard_upper), 
                width = 0.2, color = "black", size = 0.5) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "orange", alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkred", alpha = 0.7) +
  facet_wrap(~k, scales = "free_x", labeller = label_both, nrow = 1) +
  labs(title = "Bootstrap Stability by Cluster and No. of Clusters (k) with 95% Confidence Intervals",
       subtitle = "Red line: 0.75 (stable), Orange: 0.6 (partly stable), Dark red: 0.5 (unstable)",
       x = "Cluster",
       y = "Mean Jaccard Similarity (across 100 iterations)",
       caption = "Error bars represent 95% confidence intervals constructed with SD across 100 iterations") +
  ylim(0, 1) +
  theme_minimal() +
  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8))

# Display the plot
print(plot1_ci)

pdf("cluster_stab_jaccard.pdf", width=12, height=8)
print(plot1_ci)
dev.off()