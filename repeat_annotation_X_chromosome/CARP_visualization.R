###########
# TE track
###########
# This is a chunk to generate censor map plot for the TE track generated from censor of CARP analysis
# How to use:
# Run the code chunk "functions to generate TE plots" first to define the required functions
# process_censor_map("<path/to/dir>/ConsensusSequences.fa.map")
# Required input: censor map file (e.g., ConsensusSequences.fa.map)
# Output: A plot showing the TE track with annotations

# Load necessary libraries
library(tibble)
library(Gviz)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
############################################
# Define a function to generate plots
############################################
# This function generates a TE consensus plot using gviz to show the proportion of TE from RepBase that match the consensus.

# PART 1 : Function definition
# Data manipulation
process_censor_map <- function(file_path) {
  # Read the input file into a data frame
  te_map <- read.table(file_path, header = FALSE) %>%
    setNames(c("Consensus", "query_from", "query_to", "class", "match_element_from", "match_element_to",
               "dir", "sim", "pos", "score", "unknown_col1", "unknown_col2"))
  
  # Generate the mapping for Consensus names
  unique_consensus <- unique(te_map$Consensus)
  consensus_mapping <- setNames(
    paste0("Family ", seq_along(unique_consensus) - 1), 
    unique_consensus
  )
  
  # Mutate the data frame to replace the Consensus names and direction
  te_map <- te_map %>%
    mutate(
      Consensus = str_replace_all(Consensus, consensus_mapping),
      dir = str_replace_all(dir, c("d" = "+", "c" = "-"))
    )
  
  return(te_map)
}

# Plot generation function
generate_plots <- function(data, consensus_lengths) {
  # Merge the TE map with the consensus length data
  data <- data %>%
    left_join(consensus_lengths, by = "Consensus") %>% 
    select(Consensus, length, everything())
  
  # Empty list to store all plots
  all_plots <- list()

  # Loop over each unique consensus sequence
  for (consensus in unique(data$Consensus)) {
    subset_data <- filter(data, Consensus == consensus)
    
    # Create GenomeAxisTrack with the correct length for the current consensus sequence
    ax <- GenomeAxisTrack(
      from = 0, 
      to = subset_data$length[1], 
      chromosome = "chrNA", 
      labelPos = "above",        # Position the labels above the axis
      fontsize = 6,             # Increase the font size
      cex = 1.5,                 # Increase text size scaling
      littleTicks = FALSE,       # Disable little ticks to reduce clutter
      trackMargin = 12          # Adjust the track margin
    )
    
    # Create AnnotationTrack for each consensus sequence for the TE
    annTrack <- AnnotationTrack(name = consensus,
                                start = subset_data$query_from, 
                                end = subset_data$query_to,
                                strand = subset_data$dir,
                                feature = subset_data$class,
                                chromosome = "chrNA",
                                id = subset_data$class,
                                stacking = "dense",
                                shape = "fixedArrow",
                                arrowHeadWidth = 8,
                                fontsize = 6,
                                showId = TRUE)
      
    # Capture the plot output as a grob and store it in the list
    plot_grob <- grid.grabExpr(
      plotTracks(
        list(ax, annTrack), 
        size = c(4,6),
        chromosome = "chrNA", 
        showFeatureId = TRUE,
        featureAnnotation = "id", # Use "id" or "feature" to specify what feature to annotate
        fontcolor.feature = "black", # Adjust feature font color if needed
        from = 0, 
        to = subset_data$length[1]
      )
    )
    
    # Append the generated plot to the list
    all_plots <- c(all_plots, list(plot_grob))
  }
  
  # Return the list of grobs
  return(all_plots)
}

##################
# Example run 1: 
##################
# Wagyu chrmosome X centromere region in the CARP analysis with 94% identity and 400bp length

# PART 2 STEP 1:
# Read the TE map file, motify the consensus name from family00000X_consensus into Family X, and change the direction from d/c to +/-

# Use teh process_censor_map function to process the censor map file
X_ctr_carp_te_map_94_400bp <- process_censor_map("ConsensusSequences.fa.map")

# Create consensus length table
consensus_length_data_94_400bp <- tibble(Consensus = paste0("Family ", 0:7), length = c(550, 2179, 17433, 636, 1425, 2360, 5149, 26870))

# PART 2 STEP 2:
# Generate the plots by function generate_plots
carp_94_400bp_plot <- generate_plots(X_ctr_carp_te_map_94_400bp, consensus_length_data_94_400bp)

# Arrange the plots in a single figure (e.g., using grid.arrange)
n_plots <- length(carp_94_400bp_plot)
grid.arrange(grobs = carp_94_400bp_plot, 
             ncol = 1, # Display the plots in a single column
)
