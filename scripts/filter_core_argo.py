import pandas as pd

# Define chunk size
chunksize = 100_000

# Output file name
output_file = "ArgoFloats_filtered_qc12_pres1000.csv"

# Initialize flag to write header only once
first_chunk = True

# Stream-read and filter
for chunk in pd.read_csv("ArgoFloats_160E_200E_1998_2024.csv", chunksize=chunksize):
    # Apply filters
    filtered_chunk = chunk[
        chunk['PSAL_ADJUSTED_QC'].isin([1, 2]) &
        chunk['TEMP_ADJUSTED_QC'].isin([1, 2]) &
        (chunk['PRES_ADJUSTED'] <= 1000)
    ]
    
    # Write filtered rows to CSV in append mode
    filtered_chunk.to_csv(output_file, mode='a', header=first_chunk, index=False)
    first_chunk = False  # After the first write, don't write header again

print("✅ Filtered data written incrementally to:", output_file)
