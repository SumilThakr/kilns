import netCDF4 as nc
import numpy as np
import pandas as pd
import geopandas as gpd

# Define file paths
netcdf_file = "$PROJECT/GC/diff.nc4"
output_netcdf_csv = "../outputs/popwtd/netcdf_averages.csv"
output_shapefile_csv = "../outputs/popwtd/shapefile_averages.csv"

# Dictionary of shapefile paths and their scaling factors
shapefile_info = {
    "../inputs/PNO3_regridded.shp":   0.1143236952329479 / 1.1,
    "../inputs/PrimPM_regridded.shp": 0.12778081510958195,
    "../inputs/PSO4_regridded.shp":   3.319161231867947 / 1.1,
    "../inputs/SOA_regridded.shp":    1.0 / 1.05,
}

# Total Population Shapefile
pop_shapefile = "../inputs/TotalPop_regridded.shp"

### Step 1: Load Population Data
try:
    pop_gdf = gpd.read_file(pop_shapefile)
    if "MortalityRa" in pop_gdf.columns:
        population = pop_gdf["MortalityRa"].to_numpy()  # Convert to NumPy array
    else:
        raise ValueError(f"Column 'MortalityRa' not found in {pop_shapefile}")
except Exception as e:
    print(f"Error loading population shapefile: {e}")
    population = None  # Handle case where population data is missing

### Step 2: Process NetCDF File (Surface Level Only)
with nc.Dataset(netcdf_file, 'r') as dataset:
    netcdf_averages = []
    
    surface_level_idx = 0  # Assuming surface is the first level

    for var_name in dataset.variables:
        var_data = dataset.variables[var_name][:]
        
        # Ensure variable has "lev" dimension (e.g., 4D: time, lev, lat, lon)
        if var_data.ndim == 4:  
            var_data = var_data[:, surface_level_idx, :, :]  # Select surface level
        
        # Skip non-numeric and coordinate variables
        if var_data.ndim > 0 and np.issubdtype(var_data.dtype, np.number):
            simple_avg = np.nanmean(var_data)
            
            # Compute population-weighted average
            if population is not None and len(population) == var_data.size:
                var_data_flat = var_data.flatten()
                pop_weighted_avg = np.nansum(var_data_flat * population) / np.nansum(population)
            else:
                pop_weighted_avg = np.nan  # Skip if population data is unavailable
            
            netcdf_averages.append([var_name, simple_avg, pop_weighted_avg])

# Save NetCDF averages to CSV
df_netcdf = pd.DataFrame(netcdf_averages, columns=["Variable", "Simple_Average", "Population_Weighted_Average"])
df_netcdf.to_csv(output_netcdf_csv, index=False)
print(f"Saved NetCDF averages to {output_netcdf_csv}")

### Step 3: Process Shapefiles
shapefile_averages = []

for shapefile, scale_factor in shapefile_info.items():
    try:
        # Load shapefile as a GeoDataFrame
        gdf = gpd.read_file(shapefile)
        
        if "MortalityRa" in gdf.columns:
            # Divide by scaling factor
            adjusted_values = gdf["MortalityRa"] / scale_factor
            
            # Convert to NumPy array for consistency
            adjusted_values_np = adjusted_values.to_numpy()
            
            # Compute simple mean
            simple_avg = np.nanmean(adjusted_values_np)
            
            # Compute population-weighted mean
            if population is not None and len(population) == len(adjusted_values_np):
                pop_weighted_avg = np.nansum(adjusted_values_np * population) / np.nansum(population)
            else:
                pop_weighted_avg = np.nan  # Skip if population data is unavailable
            
            shapefile_averages.append([shapefile, simple_avg, pop_weighted_avg])
        else:
            print(f"Warning: 'MortalityRa' field not found in {shapefile}")

    except Exception as e:
        print(f"Error processing {shapefile}: {e}")

# Save shapefile averages to CSV
df_shapefiles = pd.DataFrame(shapefile_averages, columns=["Shapefile", "Simple_Average", "Population_Weighted_Average"])
df_shapefiles.to_csv(output_shapefile_csv, index=False)
print(f"Saved shapefile averages to {output_shapefile_csv}")
