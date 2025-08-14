# 📊 Dataset Sources for NPSG Phytoplankton Dynamics Study

This document lists the datasets used in our study, grouped by surface and depth-profile sources.

---

## 1. Surface Datasets

| No | Variable | Source | Type | Resolution | Temporal Coverage | Link |
|----|----------|--------|------|------------|--------------------|------|
| 1.1 | Surface chlorophyll, PSC (Micro, Nano, Pico) | ESA-CCI L3 (SeaWiFS, MODIS, MERIS, VIIRS, OLCI) | Satellite L-3 | 4 × 4 km | Daily (1998–2024) | [ESA CCI Chlorophyll](https://www.oceancolour.org/thredds/ncss/grid/CCI_ALL-v6.0-1km-DAILY/dataset.html) |
| 1.2 | Euphotic depth (Zeu) & Kd490 | Copernicus Ocean Colour | Satellite L-4 | 4 × 4 km | Daily (1998–2023) | [Zeu - Copernicus](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_NRT_009_102) |
| 1.3 | AOD (500–555 nm) | MERRA-2, SeaWiFS, MODIS, MISR, OMI | Satellite L-3 | 0.5°–1° | Daily/hourly (1998–2024) | [MERRA Giovanni](https://giovanni.gsfc.nasa.gov/) |
| 1.9 | Sea Level Anomaly (SLA) | DUACS Multi-Mission | Satellite L-4 | 0.125° | Daily (1998–2023) | [DUACS SLA](https://duacs.cls.fr/) |
| 1.10 | Surface geostrophic currents (ugos, vgos) | DUACS Multi-Mission | Satellite L-4 | 0.125° | Daily (1998–2023) | [DUACS Currents](https://duacs.cls.fr/) |
| 1.11 | Wind Speed | MERRA-2 | Reanalysis | 0.5° | Hourly (1998–2024) | [MERRA-2 Wind](https://giovanni.gsfc.nasa.gov/) |
| 1.12 | SST | MERRA-2 + ERSST v6 | Reanalysis | 0.5° | Hourly (1998–2024) | [MERRA-2 SST](https://giovanni.gsfc.nasa.gov/) |
| 1.13 | Cloud Fraction | MERRA-2 | Reanalysis | 0.5° | Hourly (1998–2024) | [MERRA-2 Cloud](https://giovanni.gsfc.nasa.gov/) |
| 1.14 | Precipitation & Evapotranspiration | GPM | Satellite L-3 | 0.1° | Hourly (1998–2024) | [GPM Giovanni](https://giovanni.gsfc.nasa.gov/) |
| 1.15 | PDO | ERSST v5 (1991–2020 baseline) | Satellite/In-situ | 2° | Monthly (1998–2023) | [NOAA ERSST PDO](https://www.ncei.noaa.gov/products/extended-reconstructed-sst) |

---

## 2. Depth-Profile Datasets

| No | Variable | Source | Type | Depth Range | Years | Link |
|----|----------|--------|------|-------------|-------|------|
| 2.1 | Temperature, MLD | Core Argo | In-situ | 0–1000 m | 2002–2024 | [Argo Program](https://biogeochemical-argo.org) |
| 2.2 | Chlorophyll, Nitrate, DOXY, DCM | BGC-Argo | In-situ | 0–1000 m | 2016–2024 | [BGC Argo](https://argo.ucsd.edu) |
