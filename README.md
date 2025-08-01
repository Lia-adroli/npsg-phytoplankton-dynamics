# NPSG Phytoplankton Dynamics

This repository contains code, data, and documentation for analyzing long-term physical and biogeochemical changes in the **North Pacific Subtropical Gyre (NPSG)**. The study focuses on the region between **14–28°N and 160–200°E** during **1998–2024**, integrating satellite observations, core Argo, and BGC-Argo profiles.

---

## Repository Contents

### 🔹 Python Scripts (`/scripts`)
- `fetch_bgc_floats.py` — Download BGC-Argo profiles using WMO numbers
- `filter_bgc_data.py` — Filter CHLA, DOXY, NITRATE from BGC-Argo profiles (QC 1–2)
- `merge_bgc_csvs.py` — Merge cleaned BGC float CSVs into a master dataset
- `filter_core_argo.py` — Extract 0–1000 dbar profiles from core Argo float data

### 🔹 MATLAB Analysis (`/Matlab`)
- `npsg_all_analysis.m` — Unified script with:
  1. Mixed Layer Depth (MLD) estimation
  2. Trend analysis with ±3σ outlier filtering (DOY-wise)
  3. Monthly lagged cross-correlation
  4. Granger causality testing using z-score anomalies
  5. Continuous wavelet transform (CWT) on raw Chl-a
  6. Cross-wavelet coherence with ENSO and PDO
  7. SSTa warming expansion & OHC trends (gridded, Core-Argo derived)
  8. BGC-Argo vertical trends (CHLA, NO₃⁻, DOXY)

### Reproducible Figures (`/figures`)

| Filename                             | Description                                                                                     |
|--------------------------------------|-------------------------------------------------------------------------------------------------|
| `fig1_surface_trend_profile_freq.m`  | Surface chlorophyll-a trend map (1998–2024) and Argo float profile density (Fig 1a/1b) |
| `fig2_chl_wavelet.m`                 | Continuous wavelet transform (CWT) and monthly time series of surface Chl-a anomaly (Fig 2a/2b) |
| `fig3_trend_crosscor_granger.m`      | Long-term drivers of chlorophyll variability: regression, lag, and coherence (Fig 3a–d) |
| `fig4_heat_climate.m`                | Spatial/temporal SSTa and climate patterns; cross-wavelet with Chl-a (Fig 4a–b) |
| `fig5_climate_wavelet_spatial.m`     | SSTa trends, warming expansion, SSTa-Chl relationships (Fig 5a–c) |
| `fig6_ohc_dcm_shift_bgc.m`           | Seasonal/vertical trends in SSTa, OHC, and chlorophyll + interpretation (Fig 6a/6b) |
| `fig7_bgc_hovmoller_profiles.m`      | Seasonal Hovmöller plots of CHLA, NO₃⁻, and DOXY (2016/2020 vs 2024) (Fig 7a–c) |

## Key Datasets (`/data`)

- `daily_mean.csv` — Daily merged surface variables (CHL, SSTa, AOD, SLA, etc.)
- `ArgoFloats_filtered_qc12_pres1000.csv` — QC1‑2 Core‑Argo profiles 0–1000 m
- `merged_bgc_filtered2.csv` — QC1‑2 BGC‑Argo profiles (CHLA, NO₃⁻, DOXY)
- `ersst_v6_ssta.nc` — NOAA ERSST v6 SST anomalies (1998–2024)
- `climate_chla_monthly.csv` — Monthly Chl‑a anomalies + ENSO/PDO indices

---

## 📊 Data Processing Methods

- Satellite PSC estimates derived from **ESA OC-CCI** using the Xi et al. (2021) algorithm
- AOD (550–555 nm) merged from **MERRA-2**, **SeaWiFS**, **MISR**, **MODIS-Terra**, **MODIS-Aqua**, and **OMI**
- SSTa from **NOAA ERSST v6** (1998–2024)
- Core Argo and BGC-Argo profiles retrieved via **Argopy** (Maze & Balem, 2020, expert mode)
- QC filtering: only **adjusted values** with flags **1 & 2** retained
- BGC-Argo platforms: 5904655, 5906040, 5906502, 5906516, 5906529, 3902559 → **65,265 BGC obs (0–300 m)**
- Core Argo (2002–2024): **>9.7 million T/S records** (0–1000 m)
- All in situ profiles interpolated into **10 m bins**, linearly filled and smoothed (3-point MA)
- Biogeochemical trends (CHLA, NO₃⁻, DOXY) assessed by comparing early years (2016/2020) to 2024
- OHC from Core-Argo using:  
  	**OHC = ρ × cₚ × T × Δz**  (with ρ = 1025 kg/m³, cₚ = 3985 J/kg/°C, Δz = 10 m)
- SSTa anomaly threshold: +1°C to identify warming extent
- SSTa/OHC trends: OLS regression on 1° grids; smoothed with 3-point moving avg
- Z-score anomalies used in Granger/wavelet; cross-correlation used ±3σ-filtered monthly anomalies

---

## Installation & Requirements

### Python (data preparation)
```bash
pip install argopy pandas xarray matplotlib
```

### MATLAB (analysis + figures)
Toolboxes:
- Signal Processing Toolbox (`cwt`, `wcoherence`)
- Statistics Toolbox (`fitlm`, `corr`, `granger_cause`)
- Optional: `wavelet.m`, `cross_wavelet.m` for classic implementations


MATLAB functionality includes:
- MLD calculation from sigma-T structure
- Outlier-filtered trends by DOY anomaly (3σ)
- Monthly cross-correlation (±12 lags)
- Daily Granger causality (30-day lag)
- Continuous wavelet transform (raw Chl-a)
- Cross-wavelet coherence with ENSO/PDO (monthly)
- SSTa warming analysis from ERSST v6
- Ocean Heat Content from Core Argo
- BGC-Argo vertical trends in CHLA, NO₃⁻, DOXY (2016/2020 vs 2024)

---

## Example Python Usage
```python
from argopy import DataFetcher
ds = DataFetcher(ds='bgc').float(5906502).to_xarray()
print(ds.chla.mean(dim='N_LEVELS'))
```

---

## Running MATLAB Analyses
```matlab
run_mld_analysis('data/ArgoFloats_filtered_qc12_pres1000.csv');
run_trend_analysis('data/daily_mean.csv');
run_crosscorr_analysis('data/daily_mean.csv');
run_granger_test('data/daily_mean.csv');
run_wavelet_analysis('data/daily_mean.csv');
run_crosswavelet_analysis('data/climate_chla_monthly.csv');
run_ssta_analysis('data/ssta.nc', 'data/ohc.nc');
run_bgc_argo_depth_change('data/merged_bgc_filtered2.csv');
```
---

## 🌍 Data Sources
See [`docs/dataset_sources.md`](docs/dataset_sources.md) for full dataset list.
- **BGC-Argo:** https://biogeochemical-argo.org
- **Core Argo:** https://argo.ucsd.edu
- **ESA OC-CCI:** https://esa-oceancolour-cci.org
- **CMEMS (SLA, currents):** https://marine.copernicus.eu and https://duacs.cls.fr/ 
- **MERRA-2 AOD:** https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
- **NOAA ERSST v6:** https://www.ncei.noaa.gov/products/climate-data-records/extended-reconstructed-sst

---

## 📄 License
MIT License — see [LICENSE](LICENSE) file.

---

## Author
**Lia Adroli**  
PhD Candidate, Earth Science  
Remote Sensing of Ocean Biogeochemistry  
[Research School of Earth Sciences](https://earthsciences.anu.edu.au)  
The Australian National University
