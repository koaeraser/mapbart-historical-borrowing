rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggh4x)
library(ggExtra)
library(tidyr)
mainDir <- "/Users/oliviazhang/Desktop/"

data_folder <- "data_v3"  # Options: "data" or "data_v3"

p_obs <- 10
iter_rct <- 50
iter_rwd <- 50

# Load data for all scenarios and sub-scenarios
# Try to load both null and alternative data
all_data <- data.frame()

for (sc in 1:3) {
  if (sc == 1 | sc == 2){
    # For scenarios 1 and 2, loop through all sub-scenarios
    subscenarios <- c("E")
    cor <- 1
  } else {
    # For scenario 3, sub-scenarios with correlation
    subscenarios <- c("E")
    cor <- c(-0.5, 0, 0.5)
  }

  if (sc == 1 | sc == 2) {
    for (subsc in subscenarios) {
      for (c in cor) {
        # Try both null and alternative
        for (hypo in c("null", "alternative")) {
          # Construct filename based on subsc (no cor suffix for sc == 1 or 2)
          if (subsc == "NULL") {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter_rwd,".RData")
            subsc_label <- "NULL"
          } else {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter_rwd,".RData")
            subsc_label <- subsc
          }

          # Read RCT data and extract RCT portion (D==1)
          data_rct_full <- tryCatch({
            readRDS(filename_rct)
          }, error = function(e) {
            message(paste0("RCT file not found or error reading: ", filename_rct))
            return(NULL)
          })

          # Read RWD data and extract RWD portion (D==0)
          data_rwd_full <- tryCatch({
            readRDS(filename_rwd)
          }, error = function(e) {
            message(paste0("RWD file not found or error reading: ", filename_rwd))
            return(NULL)
          })

          # Combine RCT and RWD portions
          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            # Extract RCT portion (D==1)
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y))
            # Add NA columns for U (sc == 1, 2 don't have U)
            data_rct$U1 <- NA
            data_rct$U2 <- NA
            data_rct <- data_rct[data_rct$D == 1 & data_rct$Z == 1, ]   # single-arm: RCT treated only

            # Extract RWD portion (D==0)
            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y))
            # Add NA columns for U (sc == 1, 2 don't have U)
            data_rwd$U1 <- NA
            data_rwd$U2 <- NA
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            # Combine RCT and RWD
            data <- rbind(data_rct, data_rwd)

            # Add scenario and sub-scenario labels
            data$Scenario <- sc
            data$SubScenario <- subsc_label
            data$Correlation <- NA
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  } else {
    # Scenario 3: loop through subscenarios and correlations
    for (subsc in subscenarios) {
      for (c in cor) {
        # Try both null and alternative
        for (hypo in c("null", "alternative")) {
          # Construct filename: data_p10_sc3F_cor0.5_alternative_1.RData
          if (subsc == "NULL") {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
            subsc_label <- "NULL"
          } else {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
            subsc_label <- subsc
          }

          # Read RCT data and extract RCT portion (D==1)
          data_rct_full <- tryCatch({
            readRDS(filename_rct)
          }, error = function(e) {
            message(paste0("RCT file not found or error reading: ", filename_rct))
            return(NULL)
          })

          # Read RWD data and extract RWD portion (D==0)
          data_rwd_full <- tryCatch({
            readRDS(filename_rwd)
          }, error = function(e) {
            message(paste0("RWD file not found or error reading: ", filename_rwd))
            return(NULL)
          })

          # Combine RCT and RWD portions
          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            # Extract RCT portion (D==1)
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y))
            # Add U columns if available (sc == 3), otherwise add NA columns
            if (!is.null(data_rct_full$U)) {
              data_rct <- cbind(data_rct, data_rct_full$U)
            } else {
              data_rct$U1 <- NA
              data_rct$U2 <- NA
            }
            data_rct <- data_rct[data_rct$D == 1 & data_rct$Z == 1, ]   # single-arm: RCT treated only

            # Extract RWD portion (D==0)
            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y))
            # Add U columns if available (sc == 3), otherwise add NA columns
            if (!is.null(data_rwd_full$U)) {
              data_rwd <- cbind(data_rwd, data_rwd_full$U)
            } else {
              data_rwd$U1 <- NA
              data_rwd$U2 <- NA
            }
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            # Combine RCT and RWD
            data <- rbind(data_rct, data_rwd)

            # Add scenario label, subscenario, and correlation
            data$Scenario <- sc
            data$SubScenario <- subsc_label
            data$Correlation <- c
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  }
}

# Check if any data was loaded
if (nrow(all_data) == 0) {
  stop("No data files were found. Please check that the data files exist for the specified scenarios.")
}

# Check which hypotheses were loaded
hypotheses_loaded <- unique(all_data$Hypothesis)
cat("Successfully loaded data for hypothesis:", paste(hypotheses_loaded, collapse = ", "), "\n")

# Define new color palette for groups
group_colors <- c("External Control" = "lightblue", "Treatment" = "red")

# Recreate Group factor (single-arm: no RCT control)
all_data$Group <- factor(ifelse(all_data$Z == 1, "Treatment", "External Control"),
                     levels = c("External Control", "Treatment"))

# Create factor for Hypothesis
all_data$Hypothesis <- factor(all_data$Hypothesis, levels = c("null", "alternative"))

# Create combined Group-Hypothesis variable for stratification
all_data$Group_Hypo <- interaction(all_data$Group, all_data$Hypothesis, sep = " - ")

# Create combined scenario label for faceting
all_data$Scenario_Label <- NA

# Get unique combinations that actually exist in the data
existing_sc12 <- unique(all_data[all_data$Scenario %in% c(1, 2), c("Scenario", "SubScenario")])
existing_sc3 <- unique(all_data[all_data$Scenario == 3, c("Scenario", "SubScenario", "Correlation")])

# For scenarios 1 and 2: create labels like "Scenario 1.NULL", "Scenario 1.A", etc.
if (nrow(existing_sc12) > 0) {
  for (i in 1:nrow(existing_sc12)) {
    sc <- existing_sc12$Scenario[i]
    subsc <- existing_sc12$SubScenario[i]
    idx <- which(all_data$Scenario == sc & all_data$SubScenario == subsc)
    all_data$Scenario_Label[idx] <- paste0("Scenario ", sc, ".", subsc)
  }
}

# For scenario 3: create labels like "Scenario 3.E\nρ=0.5"
if (nrow(existing_sc3) > 0) {
  for (i in 1:nrow(existing_sc3)) {
    subsc <- existing_sc3$SubScenario[i]
    c <- existing_sc3$Correlation[i]
    idx <- which(all_data$Scenario == 3 & all_data$SubScenario == subsc & all_data$Correlation == c)
    all_data$Scenario_Label[idx] <- paste0("Scenario 3.", subsc, "\nρ=", c)
  }
}

# Define factor levels in order (only for scenarios that exist)
# For scenario 3, create all combinations of subscenarios and correlations
sc3_labels <- expand.grid(subsc = c("E"),
                          cor = c(-0.5, 0, 0.5))
sc3_labels <- paste0("Scenario 3.", sc3_labels$subsc, "\nρ=", sc3_labels$cor)

all_scenario_labels <- c(paste0("Scenario 1.", c("E")),
                         paste0("Scenario 2.", c("E")),
                         sc3_labels)
# Keep only labels that actually exist in the data
existing_labels <- all_scenario_labels[all_scenario_labels %in% unique(all_data$Scenario_Label)]
all_data$Scenario_Label <- factor(all_data$Scenario_Label, levels = existing_labels)

# Reorder data so External Control is plotted last (on top)
all_data_reordered <- all_data %>% arrange(desc(Group == "External Control"))

# Create factor for Scenario to ensure proper ordering
all_data_reordered$Scenario_F <- factor(paste0("Scenario ", all_data_reordered$Scenario),
                                         levels = c("Scenario 1", "Scenario 2", "Scenario 3"))

# Create factor for SubScenario to ensure proper ordering
all_data_reordered$SubScenario_F <- factor(all_data_reordered$SubScenario,
                                            levels = c("E"))

# Create Correlation label for column faceting
all_data_reordered$Cor_Label <- ifelse(is.na(all_data_reordered$Correlation),
                                        "",
                                        paste0("ρ=", all_data_reordered$Correlation))

# Outcome distribution plots - separate by hypothesis
# Create separate lists for null and alternative plots
outcome_plots_null <- list()
outcome_plots_alt <- list()

n_scenarios <- length(levels(all_data$Scenario_Label))

for (scenario in levels(all_data$Scenario_Label)) {
  df <- all_data[all_data$Scenario_Label == scenario, ]

  # Get available hypotheses for this scenario
  hypos_in_scenario <- unique(df$Hypothesis)
  hypos_in_scenario <- hypos_in_scenario[!is.na(hypos_in_scenario)]

  # Create null plot if available
  if ("null" %in% hypos_in_scenario) {
    df_null <- df[df$Hypothesis == "null", ]

    p_null <- ggplot(df_null, aes(x = y, color = Group, fill = Group)) +
      geom_density(alpha = 0.3, size = 1.2) +
      labs(title = paste0(scenario, "\n(H0: Null)"), x = "Y", y = "Density") +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_color_manual(values = group_colors, name = "Group") +
      theme_classic(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    outcome_plots_null[[length(outcome_plots_null) + 1]] <- p_null
  }

  # Create alternative plot if available
  if ("alternative" %in% hypos_in_scenario) {
    df_alt <- df[df$Hypothesis == "alternative", ]

    p_alt <- ggplot(df_alt, aes(x = y, color = Group, fill = Group)) +
      geom_density(alpha = 0.3, size = 1.2) +
      labs(title = paste0(scenario, "\n(H1: Alternative)"), x = "Y", y = "Density") +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_color_manual(values = group_colors, name = "Group") +
      theme_classic(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    outcome_plots_alt[[length(outcome_plots_alt) + 1]] <- p_alt
  }
}

# Determine number of rows needed for outcome plots
if (length(outcome_plots_null) > 0 && length(outcome_plots_alt) > 0) {
  outcome_nrow <- 2
} else if (length(outcome_plots_null) > 0 || length(outcome_plots_alt) > 0) {
  outcome_nrow <- 1
} else {
  stop("No outcome plots were created")
}

# Calculate outcome_ncol based on actual number of plots
n_outcome_plots <- max(length(outcome_plots_null), length(outcome_plots_alt))
outcome_ncol <- n_outcome_plots

# Use patchwork to combine outcome plots
if (length(outcome_plots_null) > 0 && length(outcome_plots_alt) > 0) {
  p_outcome <- wrap_plots(outcome_plots_null, ncol = length(outcome_plots_null)) /
          wrap_plots(outcome_plots_alt, ncol = length(outcome_plots_alt))
} else if (length(outcome_plots_null) > 0) {
  p_outcome <- wrap_plots(outcome_plots_null, ncol = length(outcome_plots_null))
} else if (length(outcome_plots_alt) > 0) {
  p_outcome <- wrap_plots(outcome_plots_alt, ncol = length(outcome_plots_alt))
} else {
  stop("No outcome plots were created")
}

# Create a dummy plot with all three groups for legend
dummy_data <- data.frame(
  x = rep(1, 2),
  Group = factor(c("External Control", "Treatment"),
                 levels = c("External Control", "Treatment"))
)

legend_plot <- ggplot(dummy_data, aes(x = x, fill = Group)) +
  geom_bar() +
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

# Extract legend
legend <- get_legend(legend_plot)

# Reshape data to long format for combined covariate plot
# For sc == 1 and 2, only X5, X6 are available (confounders)
# For sc == 3, X5, U1, X6, U2 are available

# Create long format data for covariates
# Handle sc == 1, sc == 2, and sc == 3 separately
sc1_data <- all_data_reordered[all_data_reordered$Scenario == 1, ]
sc2_data <- all_data_reordered[all_data_reordered$Scenario == 2, ]
sc3_data <- all_data_reordered[all_data_reordered$Scenario == 3, ]

# For sc == 1: reshape X5, X6 and include outcome y
if (nrow(sc1_data) > 0) {
  sc1_long <- sc1_data %>%
    dplyr::select(Scenario_F, SubScenario_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc1_long$Covariate <- factor(sc1_long$Covariate, levels = c("X5", "X6"))
} else {
  sc1_long <- NULL
}

# For sc == 2: reshape X5, X6 and include outcome y
if (nrow(sc2_data) > 0) {
  sc2_long <- sc2_data %>%
    dplyr::select(Scenario_F, SubScenario_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc2_long$Covariate <- factor(sc2_long$Covariate, levels = c("X5", "X6"))
} else {
  sc2_long <- NULL
}

# For sc == 3: reshape X5, U1, X6, U2 in the specified order and include outcome y
if (nrow(sc3_data) > 0 && "U1" %in% colnames(sc3_data)) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_F, SubScenario_F, Cor_Label, Group, X5, U1, X6, U2, y) %>%
    pivot_longer(cols = c(X5, U1, X6, U2),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(sc3_long$Covariate, levels = c("X5", "U1", "X6", "U2"))
  has_U_plot <- TRUE
} else if (nrow(sc3_data) > 0) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_F, SubScenario_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(sc3_long$Covariate, levels = c("X5", "X6"))
  has_U_plot <- FALSE
} else {
  sc3_long <- NULL
  has_U_plot <- FALSE
}

# Create covariate vs outcome plot for sc == 1 with marginal distributions
if (!is.null(sc1_long) && nrow(sc1_long) > 0) {
  # Filter out rows with NA in Value or y
  sc1_long_clean <- sc1_long[!is.na(sc1_long$Value) & !is.na(sc1_long$y), ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc1_main <- ggplot(sc1_long_clean, aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(SubScenario_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Top marginal: Covariate distribution histogram (one per covariate and subscenario)
  p_cov_sc1_top <- ggplot(sc1_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(SubScenario_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenario 1", x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combine with patchwork: histogram grid above main contour
  p_cov_sc1 <- wrap_elements(p_cov_sc1_top / p_cov_sc1_main + plot_layout(heights = c(1, 1)))
  has_sc1_plot <- TRUE
} else {
  has_sc1_plot <- FALSE
}

# Create covariate vs outcome plot for sc == 2 with marginal distributions
if (!is.null(sc2_long) && nrow(sc2_long) > 0) {
  # Filter out rows with NA in Value or y
  sc2_long_clean <- sc2_long[!is.na(sc2_long$Value) & !is.na(sc2_long$y), ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc2_main <- ggplot(sc2_long_clean, aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(SubScenario_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Top marginal: Covariate distribution histogram (one per covariate and subscenario)
  p_cov_sc2_top <- ggplot(sc2_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(SubScenario_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenario 2", x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combine with patchwork: histogram grid above main contour
  p_cov_sc2 <- wrap_elements(p_cov_sc2_top / p_cov_sc2_main + plot_layout(heights = c(1, 1)))
  has_sc2_plot <- TRUE
} else {
  has_sc2_plot <- FALSE
}

# Combine sc1 and sc2 plots vertically
if (has_sc1_plot && has_sc2_plot) {
  p_cov_sc12 <- p_cov_sc1 / p_cov_sc2
  has_sc12_plot <- TRUE
} else if (has_sc1_plot) {
  p_cov_sc12 <- p_cov_sc1
  has_sc12_plot <- TRUE
} else if (has_sc2_plot) {
  p_cov_sc12 <- p_cov_sc2
  has_sc12_plot <- TRUE
} else {
  has_sc12_plot <- FALSE
}

# Create covariate vs outcome plot for sc == 3 with marginal distributions
if (!is.null(sc3_long) && nrow(sc3_long) > 0) {
  # Filter out rows with NA in Value or y
  sc3_long_clean <- sc3_long[!is.na(sc3_long$Value) & !is.na(sc3_long$y), ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc3_main <- ggplot(sc3_long_clean, aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(SubScenario_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Top marginal: Covariate distribution histogram (one per covariate and subscenario)
  p_cov_sc3_top <- ggplot(sc3_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(SubScenario_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenario 3", x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combine with patchwork: histogram grid above main contour
  p_cov_sc3 <- wrap_elements(p_cov_sc3_top / p_cov_sc3_main + plot_layout(heights = c(2, 3)))
  has_sc3_plot <- TRUE
} else {
  has_sc3_plot <- FALSE
}

# Dynamically set heights based on number of scenarios
outcome_height <- max(5.0, outcome_nrow * 1.0)

# Calculate covariate plot heights based on number of rows (SubScenario + Covariate)
n_sc1_subsc <- length(unique(sc1_data$SubScenario_F[!is.na(sc1_data$SubScenario_F)]))
n_sc2_subsc <- length(unique(sc2_data$SubScenario_F[!is.na(sc2_data$SubScenario_F)]))
n_sc3_subsc <- length(unique(sc3_data$SubScenario_F[!is.na(sc3_data$SubScenario_F)]))
# sc1 and sc2 are side by side, so use max of their heights
cov_height_sc12 <- max(n_sc1_subsc, n_sc2_subsc, n_sc3_subsc) * 10
cov_height_sc3 <- max(n_sc1_subsc, n_sc2_subsc, n_sc3_subsc) * 10

# Create a plot from the legend
legend_plot_final <- ggplot(dummy_data, aes(x = x, fill = Group)) +
  geom_bar() +
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_void() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

# Arrange the plots using patchwork (no KM plots for linear models)
if (has_sc12_plot && has_sc3_plot) {
  plot <- legend_plot_final / p_outcome / p_cov_sc12 / p_cov_sc3 +
    plot_layout(heights = c(0.5, outcome_height, cov_height_sc12, cov_height_sc3))
  total_height <- 0.1 + outcome_height + cov_height_sc12 + cov_height_sc3 + 1
} else if (has_sc12_plot) {
  plot <- legend_plot_final / p_outcome / p_cov_sc12 +
    plot_layout(heights = c(0.5, outcome_height, cov_height_sc12))
  total_height <- 0.1 + outcome_height + cov_height_sc12 + 1
} else if (has_sc3_plot) {
  plot <- legend_plot_final / p_outcome / p_cov_sc3 +
    plot_layout(heights = c(0.5, outcome_height, cov_height_sc3))
  total_height <- 0.1 + outcome_height + cov_height_sc3 + 1
} else {
  plot <- legend_plot_final / p_outcome +
    plot_layout(heights = c(0.5, outcome_height))
  total_height <- 0.1 + outcome_height + 1
}

# Calculate total plot dimensions dynamically
total_width <- outcome_ncol * 3 + 2  # Scale width with number of columns, +2 for margins

# Include data_folder in output filename to distinguish plots from different data sources
output_suffix <- ifelse(data_folder == "data", "", paste0("_", data_folder))
ggsave(paste0(file.path(mainDir),"/mapbart-sim-gaussian-realcov/inserts/balance_all_scenarios", output_suffix, ".jpg"),
       width = total_width,
       height = total_height,
       plot,
       limitsize = FALSE)
