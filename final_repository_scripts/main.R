# Load necessary libraries
library(dplyr)

### Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/new-output/"

# Read in the data frame
df <- read.csv("your_data_frame.csv")

# Extract unique species codes
species_code_list <- unique(df$species_code)

# Set the file paths for the three scripts
script1_path <- "path/to/script1.R"
script2_path <- "path/to/script2.R"
script3_path <- "path/to/script3.R"

# Source the three scripts in the desired order
source(script1_path, local = TRUE)
source(script2_path, local = TRUE)
source(script3_path, local = TRUE)