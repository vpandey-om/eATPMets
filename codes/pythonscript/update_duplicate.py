import pandas as pd
import numpy as np

# Read the Excel file
file_path = 'filtered_table_with_reactions_lb.xlsx'
data = pd.read_excel(file_path)

# Define ratio columns
ratio_cols = ['ratio_t0', 'ratio_t5', 'ratio_t30', 'ratio_t60', 'ratio_t120', 'ratio_t180']

# Utility function to handle NaN and 0 values, replacing them with a small value like 1e-5
def handle_missing_and_zero_values(data, ratio_cols, small_value=1e-5):
    """
    This function fills NaN values with a small value and replaces 0s with a small value in the ratio columns to make
    median calculation valid.
    """
    data[ratio_cols] = data[ratio_cols].replace({0: small_value, np.nan: small_value})
    return data

# Check for missing values in ratio columns and provide details
def check_for_missing_values(data, ratio_cols):
    """
    This function checks for NaN values in the ratio columns and provides details on where they occur.
    """
    nan_summary = data[ratio_cols].isna().sum()
    print("\nMissing (NaN) values in each ratio column:")
    print(nan_summary)
    return nan_summary

# Check for invalid operations (e.g., division by zero, log of zero or negative)
def check_for_invalid_operations(data, ratio_cols):
    """
    This function checks for rows where invalid operations (e.g., dividing by zero) could result in NaN.
    """
    invalid_operations = data[(data[ratio_cols] == 0).any(axis=1)]
    print("\nRows with zero values that could cause invalid operations (e.g., division by zero):")
    print(invalid_operations)
    return invalid_operations

# Check for non-numeric data
def check_data_types(data, ratio_cols):
    """
    This function checks if the ratio columns contain only numeric data.
    """
    dtype_summary = data[ratio_cols].dtypes
    print("\nData types for ratio columns:")
    print(dtype_summary)
    non_numeric_rows = data[~data[ratio_cols].applymap(np.isreal).all(axis=1)]
    print("\nRows with non-numeric data in ratio columns:")
    print(non_numeric_rows)
    return non_numeric_rows

# Check if 'Gene_Name' exists in the data
if 'Gene_Name' in data.columns:
    print("Handling duplicates in 'Gene_Name' column and calculating median for ratio columns.")
    
    # Count the number of duplicates in 'Gene_Name'
    duplicate_counts = data['Gene_Name'].value_counts()

    # Identify the duplicates (those that appear more than once)
    duplicates = duplicate_counts[duplicate_counts > 1]
    
    # If there are duplicates, calculate the median for each group
    if not duplicates.empty:
        print(f"Duplicated Gene Names found: {duplicates}")
        
        # Check for missing values in ratio columns
        nan_summary = check_for_missing_values(data, ratio_cols)
        
        # Check for invalid operations (e.g., division by zero)
        invalid_operations = check_for_invalid_operations(data, ratio_cols)
        
        # Check for non-numeric data
        non_numeric_rows = check_data_types(data, ratio_cols)
        
        # Handle NaN and 0 values in ratio columns by replacing them with small_value (1e-5)
        data = handle_missing_and_zero_values(data, ratio_cols)
        
        # Group the data by 'Gene_Name' and apply median to each group for the ratio columns
        def calc_median(group):
            # Calculate the median for each ratio column
            return pd.Series(group[ratio_cols].median(), index=ratio_cols)
        
        # Group by Gene_Name, and apply the median function while keeping Gene_Name
        grouped_data = data.groupby('Gene_Name', as_index=False).apply(lambda group: calc_median(group)).reset_index(drop=True)
        data_reduced = data.drop_duplicates(subset='Gene_Name', keep='first')

        merged_data = pd.merge(grouped_data,data_reduced, how='left', on='Gene_Name')
        # Instead of merging, directly assign the calculated medians to the original data for those Gene_Name duplicates
        # data.update(grouped_data)
     
        # Remove duplicates, keeping only the first occurrence (after median calculation)
        merged_data = merged_data.drop_duplicates(subset='Gene_Name', keep='first')
        print("Median for duplicates has been calculated and updated. Duplicates removed.")
    else:
        print("No duplicates found in 'Gene_Name' column.")
else:
    print("'Gene_Name' column is missing from the data.")

# Display the updated data
print(data.head())

# Save the updated data with median values and no duplicates
output_file = 'filtered_table_with_reactions_lb_median_no_duplicates.xlsx'
merged_data.to_excel(output_file, index=False)

print(f"Updated data with median and no duplicates saved to: {output_file}")
