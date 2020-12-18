#!/usr/env/bin/Rscript

################################################################################
#
# Function to check across simulation convergence by plotting the Gelman-Rubin
# shrinkage factor for each simulation and the median.
#
################################################################################

# Various libraries

# For my theme
library(mdiHelpR)

# For plotting
library(ggplot2)

# For pipe, filter and mutate
suppressMessages(library(dplyr))

# For string order (unnecessary)
library(stringr)

# For fread
library(data.table)

# -- Functions -----------------------------------------------------------------

# Function to take arguments from the command line using the ``optparse`` library
inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "File path to read simulation level Gelman-Rubin data in from  [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-s", "--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save plot to [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-f", "--filename"),
      type = "character",
      default = "ConvergenceAcrossChains.png",
      help = "Filename to save plot under [default= %default]",
      metavar = "character"
    ),
    
    # Scenario being analysed
    optparse::make_option(c("--scn"),
                          type = "character",
                          default = "",
                          help = "Scenario being analysed; used in plot name [default= %default]",
                          metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === Main script ==============================================================

setMyTheme()

args <- inputArguments()

# Details for reading-in and saving of files
file_dir <- args$dir
save_dir <- args$save_dir
save_name <- args$filename %>%
  paste0(save_dir, .)
scn <- args$scn %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_to_sentence()

# List files
files <- list.files(file_dir, pattern = "GelmanData.csv", full.names = T) %>%
  str_sort(numeric = T)

# Number of files present
n_files <- length(files)

# Data.frame holding all data
gelman_df <- fread(files[1])
gelman_df$Converged <- sum(gelman_df$Converged) > 0.5*nrow(gelman_df)

print(files[1])
print(str(gelman_df))

# Iterate over all files, reading in and binding to Gelman data.frame
for (f in files[2:n_files]) {
  .d <- fread(f)

  # Replaced this with calculation
  # .d$Converged <- sum(.d$Converged) > 0.5*nrow(.d)
  gelman_df <- rbind(gelman_df, .d)
}

# Find how many failed to converged based upon Brooks & Gelman, 1998
n_failed <- gelman_df %>% 
  group_by(Simulation) %>%
  summarise(Failed = any( ! Converged )) %>% 
  select(Failed) %>% 
  colSums()



# gelman_df$Converged[gelman_df$Simulation == "Simulation 6"] <- F

# Gelman plot
p_base <- gelman_df %>%
  dplyr::filter(Quantity == "median") %>%
  ggplot(aes(x = Last_iter, y = Shrinkage_factor)) +
  geom_line(aes(group = Simulation, colour = "Simulation", lty = "Simulation"))

p <- p_base +
  geom_smooth(aes(color = "Median", lty = "Median"),
    stat = "summary",
    fill = "red",
    alpha = 0.2,
    fun.data = median_hilow,
    fun.args = list(conf.int = 0.5)
  ) +
  labs(
    title = paste0(scn, ": Convergence across chains"),
    subtitle = paste0(n_failed, "% of simulations failed to converge"), #   "Including median and interquartile range",
    x = "Last iteration in chain",
    y = "Gelman-Rubin shrinkage factor",
    colour = "Values"
  ) +
  ggplot2::geom_hline(aes(yintercept = 1L, colour = "Target", linetype = "Target")) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  labs(color = "Quantity", linetype = "Quantity", shape = "Quantity")


ggsave(save_name, plot = p)


