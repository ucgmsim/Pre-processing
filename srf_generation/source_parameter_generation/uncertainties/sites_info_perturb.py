import pandas as pd
from qcore.uncertainties import distributions
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Process a CSV file.')

# Add an argument for the CSV file path
parser.add_argument('csv_file_path', help='The path to the CSV file')

# Parse the command-line arguments
args = parser.parse_args()

# Load the CSV file into a DataFrame, ensuring station name is the first column
sites_info_df = pd.read_csv(args.csv_file_path, index_col=0)

# Extract the rrup column data
rrup_data = sites_info_df["rrup"]

# Update the rrup column in the DataFrame
sites_info_df["rrup"] = distributions.truncated_log_normal(rrup_data, 0.1, 4)

# Generate a new file name with "__perturbed.csv" suffix
new_file_name = args.csv_file_path.replace(".csv", "_perturbed.csv")

# Save the updated DataFrame to a new file
sites_info_df.to_csv(new_file_name)

