# Two simulated data sets are available for download in csv format.
# However, loading CSV files into memory is not the most efficient way to go.
# To reduce loading times, convert the downloaded csv files to rds format.
# Next time, use the readRDS() function to load your data into memory.

# Step 1: After downloading the csv file, load the csv file
# Note: Do not forget to change 'YOURUSERNAME' in the path to the csv file
df_sim_250 <- read.csv("C:/Users/YOURUSERNAME/Downloads/df_sim_250.csv")
df_sim     <- read.csv("C:/Users/YOURUSERNAME/Downloads/df_sim.csv")

# Step 2: Create a subdirectory called data (and add it to the .gitignore file)

# Step 3: Save the data object
saveRDS(df_sim_250, "data/df_sim_250.rds")
saveRDS(df_sim, "data/df_sim.rds")
