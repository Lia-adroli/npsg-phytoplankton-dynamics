# NPSG Phytoplankton Dynamics

This repository provides scripts and processed datasets for investigating long-term trends in phytoplankton biomass (ocean color chlorophyll), nutrients, and hydrographic profiles in the North Pacific Subtropical Gyre (NPSG) within coordinate 14-28 degree north and 160-200 degree east year 1998-2024. It includes Argo and BGC-Argo float data retrieval, filtering, and merging routines.

## 📚 Contents

- `scripts/fetch_bgc_floats.py`: Download BGC-Argo profiles from selected WMO floats
- `scripts/filter_bgc_data.py`: Filter quality-controlled chlorophyll, DOXY, NITRATE
- `scripts/merge_bgc_csvs.py`: Merge individual float CSVs into a single file
- `scripts/filter_core_argo.py`: Filter temperature and salinity profiles to 0–1000 dbar

## 📁 Datasets

- `data/merged_bgc_filtered2.csv`: Merged and quality-filtered BGC Argo data (CHLA_ADJUSTED QC = 1 or 2) 
- `data/ArgoFloats_filtered_qc12_pres1000.csv`: Core Argo temperature/salinity profiles (QC = 1 or 2)

## 🔧 Requirements

Install required packages:

```bash
pip install argopy xarray pandas matplotlib
```

## 🧪 Example Usage

```python
from argopy import DataFetcher
ds = DataFetcher(ds='bgc').float(5906502).to_xarray()
print(ds)
```

## 🌍 Data Sources

- **BGC-Argo:** https://biogeochemical-argo.org
- **Argo Core Data:** https://argo.ucsd.edu

## 📄 License

This project is shared under the [MIT License](LICENSE).

## 👩‍🔬 Author

Lia Adroli — PhD candidate in Earth Science | Remote Sensing of Ocean Biogeochemistry | Research School of Earth Sciences | the Australian National University

