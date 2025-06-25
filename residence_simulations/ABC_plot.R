library(dplyr)
library(ggplot2)
library(stringr)

# Read observed data
observed_data <- read.table("Observed.txt", header = TRUE, sep = "\t")

# Function to read simulation data
read_simulation_data <- function(gen_number, family_type = "") {
  sim_dirs <- c("../Simulations1/", "../Simulations2/")
  
  # For each simulation folder, for each replicate, for each j value
  sim_list <- lapply(sim_dirs, function(sim_dir) {
    lapply(1:10, function(i) {
      lapply(seq(0, 100, by = 10), function(j) {
        file_path <- paste0(sim_dir, i, "/matrilocal.", j, 
                            ".patrilocal.", 100 - j, ".", gen_number, ".50", family_type, ".tsv")
        read.delim(file_path, header = TRUE)
      })
    })
  })
  # Flatten the list of lists and bind the rows
  sim_data <- bind_rows(unlist(sim_list, recursive = FALSE))
  
  # Remove unwanted columns
  sim_data <- sim_data %>% select(-Generation, -Simulation)
  return(sim_data)
}

# Function to process simulation data and calculate distances
process_simulation_data <- function(sim_data, reference_data) {
  full_data <- bind_rows(reference_data, sim_data)
  
  numeric_columns <- c("AuXThetaDiffWB", "YHaploHomozygosityWB", "mtHaploHomozygosityWB", "MatrilocalProp")
  full_data[numeric_columns] <- lapply(full_data[numeric_columns], as.numeric)
  
  # Define a normalization function
  normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  normalized_data <- full_data %>%
    mutate(
      AuXThetaDiffWB = normalize(AuXThetaDiffWB),
      YHaploHomozygosityWB = normalize(YHaploHomozygosityWB),
      mtHaploHomozygosityWB = normalize(mtHaploHomozygosityWB)
    )
  
  # Separate the reference data and simulation data
  normalized_reference_data <- normalized_data[1, ]
  normalized_simulation_data <- normalized_data[-1, ]
  
  ref_values <- normalized_reference_data[numeric_columns]
  
  normalized_simulation_data <- normalized_simulation_data %>%
    mutate(
      Distance = sqrt(
        (AuXThetaDiffWB - ref_values$AuXThetaDiffWB)^2 +
          (YHaploHomozygosityWB - ref_values$YHaploHomozygosityWB)^2 +
          (mtHaploHomozygosityWB - ref_values$mtHaploHomozygosityWB)^2
      )
    )
  
  sorted_distances <- sort(normalized_simulation_data$Distance)
  
  return(list(
    normalized_simulation_data = normalized_simulation_data,
    sorted_distances = sorted_distances
  ))
}

# Define the scenarios in a list
scenarios <- list(
  OneGen4 = list(gen_number = 4, family_type = ""),
  OneGen3 = list(gen_number = 3, family_type = ""),
  OneGen2 = list(gen_number = 2, family_type = ""),
  TwoGen4 = list(gen_number = 4, family_type = ".twofamily"),
  TwoGen3 = list(gen_number = 3, family_type = ".twofamily"),
  TwoGen2 = list(gen_number = 2, family_type = ".twofamily")
)

# Process the simulations for each scenario
simulation_results <- list()
for (sim in names(scenarios)) {
  params <- scenarios[[sim]]
  cat("Processing scenario:", sim, "\n")
  sim_data <- read_simulation_data(gen_number = params$gen_number, family_type = params$family_type)
  simulation_results[[sim]] <- process_simulation_data(sim_data, observed_data)
}

# Acceptance level for posterior probabilities (e.g., top 5%)
acceptance_level <- 0.05
posterior_probabilities_list <- list()

# Compute posterior probabilities for each simulation scenario
for (sim in names(simulation_results)) {
  sim_result <- simulation_results[[sim]]
  sim_data <- sim_result$normalized_simulation_data
  sorted_distances <- sim_result$sorted_distances
  
  threshold_index <- ceiling(acceptance_level * length(sorted_distances))
  threshold <- sorted_distances[threshold_index]
  
  accepted_simulations <- sim_data %>% filter(Distance <= threshold)
  
  posterior_probabilities <- accepted_simulations %>%
    group_by(MatrilocalProp) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Probability = Count / sum(Count),
      Acceptance = acceptance_level,
      Threshold = threshold,
      simtype = sim
    )
  
  posterior_probabilities_list[[sim]] <- posterior_probabilities
}

posterior_probabilities_all <- bind_rows(posterior_probabilities_list)

# Ensure MatrilocalProp is treated as a factor with desired levels
posterior_probabilities_all$MatrilocalProp <- factor(
  posterior_probabilities_all$MatrilocalProp,
  levels = seq(0, 100, by = 10)
)

# Create the plot for posterior probabilities
ObservedPlot <- ggplot(
  posterior_probabilities_all %>% filter(Probability > 0),
  aes(x = MatrilocalProp, y = Probability, group = simtype, color = simtype)
) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_x_discrete(drop = FALSE) +
  ylim(0, 1) +
  scale_color_manual(
    values = c(
      'OneGen2' = "pink", 'OneGen3' = "#f56c6c", 'OneGen4' = "#a83232",
      'TwoGen2' = "lightblue", 'TwoGen3' = "#2f74a1", 'TwoGen4' = "darkblue"
    ),
    labels = c(
      'OneGen2' = '1 Family - 2 Generations', 'OneGen3' = '1 Family - 3 Generations',
      'OneGen4' = '1 Family - 4 Generations', 'TwoGen2' = '2 Families - 2 Generations',
      'TwoGen3' = '2 Families - 3 Generations', 'TwoGen4' = '2 Families - 4 Generations'
    )
  ) +
  theme_minimal() +
  labs(
    x = "Matrilocality rate",
    y = "Posterior probability",
    color = "Simulation scenarios",
    title = NULL
  ) +
  guides(color = guide_legend(nrow = 6, byrow = TRUE)) +
  theme(
    legend.position = c(0.22, 0.82),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 12, family = "sans"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12, family = "sans"),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.background = element_rect(fill = "white", linewidth = 0.5, linetype = "solid"),
    plot.background = element_rect(fill = "white")
  )

# Save the plot to a file
ggsave("Results.png", ObservedPlot, width = 20, height = 20, units = "cm")
