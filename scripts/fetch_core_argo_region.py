from argopy import DataFetcher
import pandas as pd

# Define years to loop over
years = list(range(1998, 2025))
data_all = []

for year in years:
    start = f"{year}-01"
    end = f"{year}-12"

    # Define both longitude boxes to span 160E to 200E (across dateline)
    box1 = [160, 180, 14, 28, 0, 2000, start, end]
    box2 = [-180, -160, 14, 28, 0, 2000, start, end]

    try:
        df1 = DataFetcher().region(box1).to_dataframe()
        df2 = DataFetcher().region(box2).to_dataframe()
        data_all.append(pd.concat([df1, df2]))
        print(f"✓ Collected data for {year}")
    except Exception as e:
        print(f"⚠️ Failed to fetch for {year}: {e}")

# Combine all dataframes
df_total = pd.concat(data_all)

# Save to CSV
df_total.to_csv("ArgoFloats_160E_200E_1998_2024.csv", index=False)
print("✅ Saved to 'ArgoFloats_160E_200E_1998_2024.csv'")
