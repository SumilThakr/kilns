import numpy as np
import pandas as pd
import geopandas as gpd

# Define file paths
output_shapefile_csv = "../outputs/popwtd/shapefile_averages_inmapgrid.csv"

# Dictionary of shapefile paths and their scaling factors
shapefile_info = {
    "../inputs/PNO3_native.shp":   0.1143236952329479 / 1.1,
    "../inputs/PrimPM_native.shp": 0.12778081510958195,
    "../inputs/PSO4_native.shp":   3.319161231867947 / 1.1,
    "../inputs/SOA_native.shp":    1.0 / 1.05,
}

# Total Population Shapefile
pop_shapefile = "$PROJECT/concentrations/TotalPop_regridded_native.shp"

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

### Step 3: Process Shapefiles
shapefile_averages = []

for shapefile, scale_factor in shapefile_info.items():
    try:
        # Load shapefile as a GeoDataFrame
        gdf = gpd.read_file(shapefile)
        
        if "MortalityR" in gdf.columns:
            # Divide by scaling factor
            adjusted_values = gdf["MortalityR"] / scale_factor
            
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
            print(f"Warning: 'MortalityR' field not found in {shapefile}")

    except Exception as e:
        print(f"Error processing {shapefile}: {e}")

# Save shapefile averages to CSV
df_shapefiles = pd.DataFrame(shapefile_averages, columns=["Shapefile", "Simple_Average", "Population_Weighted_Average"])
df_shapefiles.to_csv(output_shapefile_csv, index=False)
print(f"Saved shapefile averages to {output_shapefile_csv}")
