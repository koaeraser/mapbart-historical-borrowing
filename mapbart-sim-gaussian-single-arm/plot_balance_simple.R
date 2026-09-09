rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
library(patchwork)
library(tidyr)

mainDir <- "/Users/oliviazhang/Desktop/"

data_folder <- "data_v3"  # Options: "data" or "data_v3"

p_obs <- 10
iter_rct <- 50
iter_rwd <- 50

# Load data for all scenarios (without sub-scenario distinction)
all_data <- data.frame()

for (sc in 1:3) {
  if (sc == 1 | sc == 2){
    subscenarios <- c("E")
    cor <- 1
  } else {
    subscenarios <- c("E")
    cor <- c(-0.5, 0, 0.5)
  }

  if (sc == 1 | sc == 2) {
    for (subsc in subscenarios) {
      for (c in cor) {
        for (hypo in c("null", "alternative")) {
          if (subsc == "NULL") {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter_rwd,".RData")
          } else {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter_rwd,".RData")
          }

          data_rct_full <- tryCatch({
            readRDS(filename_rct)
          }, error = function(e) {
            return(NULL)
          })

          data_rwd_full <- tryCatch({
            readRDS(filename_rwd)
          }, error = function(e) {
            return(NULL)
          })

          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y))
            data_rct$U1 <- NA
            data_rct$U2 <- NA
            data_rct <- data_rct[data_rct$D == 1, ]

            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y))
            data_rwd$U1 <- NA
            data_rwd$U2 <- NA
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            data <- rbind(data_rct, data_rwd)
            data$Scenario <- sc
            data$Correlation <- NA
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  } else {
    for (subsc in subscenarios) {
      for (c in cor) {
        for (hypo in c("null", "alternative")) {
          if (subsc == "NULL") {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
          } else {
            filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
          }

          data_rct_full <- tryCatch({
            readRDS(filename_rct)
          }, error = function(e) {
            return(NULL)
          })

          data_rwd_full <- tryCatch({
            readRDS(filename_rwd)
          }, error = function(e) {
            return(NULL)
          })

          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y))
            if (!is.null(data_rct_full$U)) {
              data_rct <- cbind(data_rct, data_rct_full$U)
            } else {
              data_rct$U1 <- NA
              data_rct$U2 <- NA
            }
            data_rct <- data_rct[data_rct$D == 1, ]

            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y))
            if (!is.null(data_rwd_full$U)) {
              data_rwd <- cbind(data_rwd, data_rwd_full$U)
            } else {
              data_rwd$U1 <- NA
              data_rwd$U2 <- NA
            }
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            data <- rbind(data_rct, data_rwd)
            data$Scenario <- sc
            data$Correlation <- c
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  }
}

if (nrow(all_data) == 0) {
  stop("No data files were found.")
}

# Define color palette
group_colors <- c("Control" = "darkblue", "External Control" = "lightblue", "Treatment" = "red")

# Create Group factor
all_data$Group <- factor(ifelse(all_data$Z == 1, "Treatment",
                            ifelse(all_data$Z == 0 & all_data$D == 1, "Control", "External Control")),
                     levels = c("Control", "External Control", "Treatment"))

all_data$Hypothesis <- factor(all_data$Hypothesis, levels = c("null", "alternative"))

# Create scenario label (with rho for Scenario 3)
all_data$Scenario_Label <- ifelse(all_data$Scenario == 3,
                                   paste0("Scenario 3\n(ρ=", all_data$Correlation, ")"),
                                   paste0("Scenario ", all_data$Scenario))
all_data$Scenario_Label <- factor(all_data$Scenario_Label,
                                   levels = c("Scenario 1", "Scenario 2",
                                              "Scenario 3\n(ρ=-0.5)", "Scenario 3\n(ρ=0)", "Scenario 3\n(ρ=0.5)"))

# Reorder data so External Control is plotted last (on top)
all_data_reordered <- all_data %>% arrange(desc(Group == "External Control"))

# ============================================================================
# PLOT 1: Outcome Density Plot
# ============================================================================

# Create density plots for each scenario (Sc1, Sc2, and Sc3 with 3 rho values)
outcome_plots <- list()

# Scenarios 1 and 2
for (sc in 1:2) {
  df <- all_data[all_data$Scenario == sc, ]
  scenario_label <- paste0("Scenario ", sc)

  p <- ggplot(df, aes(x = y, color = Group, fill = Group)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +
    labs(title = scenario_label, x = "Y", y = "Density") +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

  outcome_plots[[length(outcome_plots) + 1]] <- p
}

# Scenario 3 with separate plots for each rho
for (rho in c(-0.5, 0, 0.5)) {
  df <- all_data[all_data$Scenario == 3 & all_data$Correlation == rho, ]
  scenario_label <- paste0("Scenario 3 (ρ=", rho, ")")

  p <- ggplot(df, aes(x = y, color = Group, fill = Group)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +
    labs(title = scenario_label, x = "Y", y = "Density") +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

  outcome_plots[[length(outcome_plots) + 1]] <- p
}

# Combine outcome plots (no legend)
p_outcome_final <- wrap_plots(outcome_plots, ncol = 5)

# Save outcome density plot
output_suffix <- ifelse(data_folder == "data", "", paste0("_", data_folder))
ggsave(paste0(file.path(mainDir),"/mapbart-sim-gaussian/manuscript/inserts/balanced_1_real.jpg"),
       width = 16,
       height = 3,
       p_outcome_final,
       limitsize = FALSE)

# ============================================================================
# PLOT 2: Covariate Balance Plot
# ============================================================================

# Reshape data for covariate histograms
# For sc == 1 and 2: X5, X6
# For sc == 3: X5, U1, X6, U2

sc1_data <- all_data_reordered[all_data_reordered$Scenario == 1, ]
sc2_data <- all_data_reordered[all_data_reordered$Scenario == 2, ]
sc3_data <- all_data_reordered[all_data_reordered$Scenario == 3, ]

# Define covariate label mapping
cov_labels <- c("X5" = "X₅ (Measured)", "X6" = "X₆ (Measured)",
                "U1" = "U₁ (Unmeasured)", "U2" = "U₂ (Unmeasured)")

# Reshape to long format
if (nrow(sc1_data) > 0) {
  sc1_long <- sc1_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc1_long$Covariate <- factor(cov_labels[sc1_long$Covariate],
                                levels = c("X₅ (Measured)", "X₆ (Measured)"))
} else {
  sc1_long <- NULL
}

if (nrow(sc2_data) > 0) {
  sc2_long <- sc2_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc2_long$Covariate <- factor(cov_labels[sc2_long$Covariate],
                                levels = c("X₅ (Measured)", "X₆ (Measured)"))
} else {
  sc2_long <- NULL
}

if (nrow(sc3_data) > 0 && "U1" %in% colnames(sc3_data) && !all(is.na(sc3_data$U1))) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_Label, Group, X5, U1, X6, U2) %>%
    pivot_longer(cols = c(X5, U1, X6, U2),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(cov_labels[sc3_long$Covariate],
                                levels = c("X₅ (Measured)", "U₁ (Unmeasured)",
                                           "X₆ (Measured)", "U₂ (Unmeasured)"))
} else if (nrow(sc3_data) > 0) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(cov_labels[sc3_long$Covariate],
                                levels = c("X₅ (Measured)", "X₆ (Measured)"))
} else {
  sc3_long <- NULL
}

# Combine all scenarios
all_long <- rbind(sc1_long, sc2_long, sc3_long)
all_long <- all_long[!is.na(all_long$Value), ]

# Create covariate balance histogram plot
p_balance <- ggplot(all_long, aes(x = Value, fill = Group)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity") +
  facet_grid(Covariate ~ Scenario_Label, scales = "free") +
  scale_fill_manual(values = group_colors, name = "Group") +
  labs(x = "Covariate Value", y = "Density") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

# Determine height based on number of covariates (no legend)
n_covariates <- length(unique(all_long$Covariate))
plot_height <- n_covariates * 2

# Save covariate balance plot (5 columns: Sc1, Sc2, Sc3 with 3 rho values)
ggsave(paste0(file.path(mainDir),"/mapbart-sim-gaussian/manuscript/inserts/balanced_2_real.jpg"),
       width = 14,
       height = plot_height,
       p_balance,
       limitsize = FALSE)

