#!/usr/env/bin/Rscript

################################################################################
#
# Function to check within simulation convergence by plotting the Gelman-Rubin 
# shrinkage factor (for across chain convergence) and the Geweke Z-scores (for 
# within chain convergence).
#
################################################################################

# Various libraries
library(ggplot2)
library(mdiHelpR)
library(magrittr)
library(coda)
suppressMessages(library(mcclust))
library(data.table)
library(patchwork)
library(stringr)
suppressMessages(library(dplyr))
library(optparse)
suppressMessages(library(stableGR))

# Function to take arguments from the command line using the ``optparse`` library
inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "File path to model output [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-s", "--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save model and predicted clustering to [default= %default]",
      metavar = "character"
    ),

    # Expected chain length
    optparse::make_option(c("-l", "--length"),
      type = "integer",
      default = 1e6,
      help = "Expected chain length; will call an error if any chains are less than this [default= %default]",
      metavar = "integer"
    ),

    # Thinning factor applied within each chain
    optparse::make_option(c("-t", "--thin"),
      type = "integer",
      default = 1e3,
      help = "Thinning factor applied within each chain [default= %default]",
      metavar = "integer"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("-b", "--burn_in"),
      type = "integer",
      default = 1e4,
      help = "Burn-in to apply to each chain (accounts for the already-applied thinning) [default= %default]",
      metavar = "integer"
    ),

    # Columns of interest for convergence
    optparse::make_option(c("-c", "--cols"),
      type = "character",
      default = "1",
      help = "Columns containing variables of interest for convergence tests [default= %default]",
      metavar = "character"
    ),

    # Scenario name
    optparse::make_option(c("--scn"),
      type = "character",
      default = "",
      help = "Scenario currently being tested [default= %default]",
      metavar = "character"
    ),

    # Scenario name
    optparse::make_option(c("--threshold"),
      type = "numeric",
      default = 1.01,
      help = "Threshold for across chain convergence for R hat. [default= %default]",
      metavar = "numeric"
    )



  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# Read in arguments from the command line
args <- inputArguments()

# Directories to load and save to/from
analysis_dir <- args$dir
save_dir <- args$save_dir

# Columns containing continuous random variables
cols_to_keep <- args$cols %>%
  strsplit(" ") %>%
  unlist() %>%
  as.numeric()

n_params <- length(cols_to_keep)

# Parameters regarding chain lengths and samples to keep
expected_length <- args$length %>%
  as.numeric() 

burnin <- args$burn_in %>%
  as.numeric()

thin_by <- args$thin %>%
  as.numeric()

# Scenario being analysed
scn <- args$scn

# Gelman rubin shrinkage factor convergence threshold
gr_thresh <- args$threshold

# effective values for burn-in and length due to thinning
# +1 as 0th iteration is saved
eff_burnin <- floor(burnin/thin_by) + 1
eff_expected_length <- floor(expected_length/thin_by) + 1

# My dewfault ggplot2 theme
setMyTheme()

# Read the files and sort by chain number (recogininsing as a numeric)
my_files <- list.files(analysis_dir, full.names = T) %>%
  stringr::str_sort(numeric = T)

file_order <- my_files %>% stringr::str_order(numeric = T)
n_files <- length(my_files)

# Extract the simulation number (used in saving results)
sim_name <- my_files %>%
  str_extract("simulation_[1-9][0-9]*")

# If any of the files claim to belong to a different simulation, stop
if (length(unique(sim_name)) != 1) {
  stop("Multiple simulations present in file names.")
} else {
  sim_name <- sim_name[1]
}

# Current simulation string
sim_string <- sim_name %>% 
  str_replace("_", " ") %>% 
  str_to_sentence()

samples <- list()
for (i in 1:n_files) {
  samples[[i]] <- fread(my_files[i], select = cols_to_keep)
}

# Check all the chains have the same number of samples; stop analysis if not
fail <- sapply(samples, nrow) != eff_expected_length
if (any(fail)) {
  print(paste("Chain(s)", paste(which(fail), collapse = " & "), "has less samples than expected."))
  stop("Stopping analysis.")
}

# Apply burn in and thinning, and convert to mcmc type (for coda)
convergence_data <- plt_data <- list()
for (i in 1:n_files) {
  if (n_params == 1) {

    convergence_data[[i]] <- plt_data[[i]] <- samples[[i]][seq(eff_burnin, eff_expected_length, 1), ] 
    convergence_data[[i]] <- convergence_data[[i]] %>%
      set_rownames(seq(burnin, expected_length-1, thin_by)) %>%
      set_colnames("Mass Parameter") %>%
      mcmc(start=burnin, end=expected_length-1, thin = thin_by)
  } else {
    convergence_data[[i]] <- plt_data[[i]] <- samples[[i]][seq(burnin, expected_length, 1), ] 
    convergence_data[[i]] <- convergence_data[[i]] %>%
      set_rownames(seq(burnin, expected_length, thin_by)) %>%
      mcmc(start = burnin, end = expected_length)
  }
  plt_data[[i]]$Chain <- paste0("Chain ", i)
}



# First plot the Geweke Z-score 
# Get the geweke data for chains
p_geweke <- lapply(
  convergence_data,
  gewekePlot,
  threshold_line_colour = "red"
)

for (i in 1:n_files) {
  curr_df <- p_geweke[[i]]$data
  curr_df$Chain <- paste0("Chain ", i)
  
  if (i == 1) {
    geweke_df <- curr_df
  } else {
    geweke_df <- rbind(geweke_df, curr_df)
  }
}

# Test within chain convergence
# Set as a factor and sorted numerically to have the facets ordered nicely
geweke_df$Chain <- geweke_df$Chain %>%
  factor(levels = stringr::str_sort(unique(.), numeric = T))

# Test convergence within chain using SW to test for Normalisty
test_data <- geweke_df %>%
  group_by(Chain) %>% 
  summarise(Shapiro_p_value = shapiro.test(Geweke_statistic)$p.value) %>% 
  mutate(Normal = Shapiro_p_value > 0.05)

# Add this to the geweke df which will be saved
geweke_df$Converged <- test_data$Normal[match(geweke_df$Chain, test_data$Chain)]

n_converged <- sum(test_data$Normal)

# Take the converged chains forward for plotting
converged_data <- convergence_data[which(test_data$Normal)]

# Test GR convergence using Vats and Knudson's updated statistic 
across_chain_conv <- stable.GR(convergence_data,
  multivariate = F, 
  autoburnin = F, 
  blather = F,
  size = 991 
)

chains_converged <- across_chain_conv$psrf < gr_thresh

# Data for histograms and whatnot
plt_data <- plt_data[which(test_data$Normal)] 
plt_data <- do.call(rbind, plt_data)

plt_data$Chain <- plt_data$Chain %>%
  factor(levels = stringr::str_sort(unique(.), numeric = T))

print(length(seq(burnin, expected_length-1, thin_by)))
# Gelman-Rubin plot
p1 <- gelmanPlot(convergence_data) +
  labs(
    title = scn,
    subtitle = paste0(sim_string, ": Across chain convergence"))  +
  scale_color_manual(labels = c("97.5%", "Median"),
                     values = c("grey", "black")) +
  scale_linetype_manual(labels = c("97.5%", "Median"),
                       values = c(2, 1))

  # geom_hline(aes(yintercept = 1.1, lty = "Convergence", colour = "Convergence")) +
  # scale_color_manual(labels = c("97.5%", "Convergence", "Median"),
  #                    values = c("grey", "blue", "black")) +
  # scale_linetype_manual(labels = c("97.5%", "Convergence", "Median"),
  #                      values = c(2, 2, 1))

# + 
#  scale_x_continuous(label = function(x) {
#    return(paste(x * thin_by + burnin))
#  })

gelman_data <- p1$data %>% 
  mutate(Simulation = sim_string)

# Define convergence on Vats, Knudson, 2020 (extended from  Brooks and Gelman, 1997)
gelman_data$Converged <- chains_converged # gelman_data$Shrinkage_factor < 1.2
gelman_data$PSRF <-  across_chain_conv$psrf 

# Geweke plot for first chain
p2 <- p_geweke[[1]] +
  labs(
    subtitle = paste0(sim_string, ": Within chain convergence for chain 1"),
    title = scn
  )# +
#  scale_x_continuous(label = function(x) {
#    return(paste(x * thin_by + burnin))
#  }) # ,
# breaks = scales::pretty_breaks(3)            #### Note that 3 might be a bad choice
# )

p3 <- plt_data %>% 
  ggplot(aes(x = MassParameter_1)) + #, colour = Chain)) +
  geom_histogram() + # geom_density(alpha = 0.3) +
  # scale_fill_viridis_d() +
  facet_wrap(~Chain) +
  labs(
    subtitle = paste0(sim_string, ": Sampled density of mass parameter across chains"),
    title = scn, # "Coloured by chain of origin",
    x = "Mass parameter",
    y = "Count" # y = "Density"
  ) +
  theme(legend.position = "none")

# p3

# Trace plot of mass parameter
p4 <- plt_data %>% 
  mutate(Iteration = rep(seq(burnin, expected_length, thin_by), n_converged)) %>% 
  ggplot(aes(x = Iteration, y = MassParameter_1)) +
  geom_line() +
  facet_wrap(~Chain) +
  labs(
    title = scn,
    subtitle =  paste0(sim_string, ": Mass parameter samples across chains"),
    x = "Iteration",
    y = "Mass parameter"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)
  )

# p4

# Drop the first simulation to have a nice 3x3 grid.
plot_geweke_df <- geweke_df %>% dplyr::filter(Chain != "Chain 1")

# The 1.96 threshold for 0.05 significance on a standard normal distribution
c_limit <- stats::qnorm(1 - 0.05 / 2)

# Construct facet wrapped Geweke plots
y_limit <- max(c(c_limit, abs(geweke_df$Geweke_statistic)))
threshold_line_colour <- "red"
plt_title <- paste0(sim_string, ": Within chain convergence for remaining chains")

# Plot the rest of the chains within-chain convergence
p_rest <- ggplot(plot_geweke_df, aes(x = Start_iteration, y = Geweke_statistic)) +
  geom_line() +
  geom_hline(yintercept = c_limit, linetype = "dashed", color = threshold_line_colour) +
  geom_hline(yintercept = -c_limit, linetype = "dashed", color = threshold_line_colour) +
  facet_wrap(~Chain) +
  labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = plt_title,
    subtitle = "Mass parameter"
  ) +
#  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(3),
    limits = c(-y_limit, y_limit)
  ) +
  theme(axis.text.x=element_text(angle = 30, hjust = 1))

# Use patchwork to create a grid of plots
p_patch <- (p1 + p2) / p_rest

# Save the plots
ggsave(paste0(save_dir, scn, sim_name, "Convergence.png"), plot = p_patch)
ggsave(paste0(save_dir, scn, sim_name, "Gelman.png"), plot = p1)
ggsave(paste0(save_dir, scn, sim_name, "Geweke1.png"), plot = p2)
ggsave(paste0(save_dir, scn, sim_name, "Geweke2t10.png"), plot = p_rest)
ggsave(paste0(save_dir, scn, sim_name, "Densities.png"), plot = p3)
ggsave(paste0(save_dir, scn, sim_name, "TracePlots.png"), plot = p4)

# Save the data that generated them
write.csv(geweke_df, paste0(save_dir, scn, sim_name, "GewekeData.csv"), row.names = F)
write.csv(gelman_data, paste0(save_dir, scn, sim_name, "GelmanData.csv"), row.names = F)

