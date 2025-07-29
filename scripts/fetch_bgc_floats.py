from argopy import DataFetcher, set_options
import xarray as xr

# Enable expert mode to access full BGC variable set
set_options(mode='expert')

# List of BGC Argo float WMO numbers
wmo_list = [2902822, 2903860, 3902559, 5904655, 5906040, 5906502, 5906516, 5906529]

def main():
    fetcher = DataFetcher(ds='bgc')  # use BGC dataset
    for wmo in wmo_list:
        print(f"\n=== Fetching BGC float {wmo} ===")
        try:
            ds: xr.Dataset = fetcher.float(wmo).to_xarray()
            print(ds)  # quick overview
            print(f"✅ Retrieved {wmo} with {len(ds['N_PROF']) if 'N_PROF' in ds.dims else 'unknown'} profiles")
        except Exception as e:
            print(f"❌ Failed to fetch {wmo}: {e}")

if __name__ == "__main__":
    main()
