import pandas as pd
import glob
import os

# Ensure output folder exists
os.makedirs("data", exist_ok=True)

files = glob.glob("bgc_float_*.csv")
columns_to_keep = [
    'TIME', 'LATITUDE', 'LONGITUDE', 'PLATFORM_NUMBER',
    'PRES', 'PSAL', 'TEMP', 'CHLA','CHLA_QC', 'CHLA_ADJUSTED',
    'CHLA_ADJUSTED_QC', 'CDOM', 'DOXY', 'DOXY_ADJUSTED',
    'NITRATE', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED'
]

dfs = []
for file in files:
    try:
        df = pd.read_csv(file)
        if 'CHLA_ADJUSTED_QC' not in df.columns:
            print(f"⚠️ Skipped {file} (no CHLA_ADJUSTED_QC column)")
            continue
        # Keep only good QC for CHLA_ADJUSTED
        df = df[df['CHLA_ADJUSTED_QC'].astype(str).isin(['1', '2'])]
        available = [c for c in columns_to_keep if c in df.columns]
        df = df[available]
        # Normalize TIME column name if needed
        if 'time' in df.columns and 'TIME' not in df.columns:
            df['TIME'] = pd.to_datetime(df['time'], errors='coerce')
            df.drop(columns=['time'], inplace=True)
        dfs.append(df)
    except Exception as e:
        print(f"❌ Error processing {file}: {e}")

if dfs:
    merged = pd.concat(dfs, ignore_index=True)
    out = "data/merged_bgc_filtered2.csv"
    merged.to_csv(out, index=False)
    print(f"✅ Saved merged CSV to {out} with {len(merged):,} rows")
else:
    print("No input files matched pattern bgc_float_*.csv")
