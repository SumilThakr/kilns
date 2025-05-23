import geopandas as gpd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

# Load data from shapefile
shapefile_path = "../intermediate/total_PM25.shp"
gdf = gpd.read_file(shapefile_path)
mortality_ra = gdf["Total_PM25"].values

# Load data from NetCDF file
netcdf_path = "$PROJECT/GC/diff.nc4"
dataset = nc.Dataset(netcdf_path)
aermass_bc = dataset.variables["PM25"][0, 0, :, :]  # Select time=0, level=0

# Flatten the NetCDF data to match the shapefile structure
aermass_bc_flat = aermass_bc.flatten()

# Ensure the data lengths match
if len(mortality_ra) != len(aermass_bc_flat):
    raise ValueError("Mismatch in grid cell lengths between shapefile and NetCDF data.")

# Scatterplot
plt.figure(figsize=(8, 8))
plt.scatter(aermass_bc_flat, mortality_ra, label="Data Points", alpha=0.6)

# Compute slope manually (forcing intercept = 0)
slope = np.sum(aermass_bc_flat * mortality_ra) / np.sum(aermass_bc_flat ** 2)

# Generate predicted values
predicted = slope * aermass_bc_flat

# Plot the regression line
plt.plot(aermass_bc_flat, predicted, color="red", label="Best Fit Line (Intercept=0)")

# Plot y = x line
x_vals = np.linspace(min(aermass_bc_flat), max(aermass_bc_flat), 100)
plt.plot(x_vals, x_vals, color="blue", linestyle="--", label="y = x")

# Calculate statistics
SSE = np.sum((mortality_ra - predicted) ** 2)
TSS = np.sum(mortality_ra ** 2)
r2 = 1 - (SSE / TSS)
mean_error = np.mean(np.abs(mortality_ra - predicted))
mean_bias = np.mean(mortality_ra - predicted)
normalized_mean_error = (np.sum(np.abs(mortality_ra - aermass_bc_flat))/np.sum(aermass_bc_flat)) *100.0
normalized_mean_bias = (np.sum(mortality_ra - aermass_bc_flat)/np.sum(aermass_bc_flat)) *100.0

# Annotate the statistics on the plot
stats_text = (f"$R^2$: {r2:.2f}\n"
              f"Slope: {slope:.2f}\n"
              f"Normalized Mean Error: {normalized_mean_error:.1f}%\n"
              f"Normalized Mean Bias: {normalized_mean_bias:.1f}%")
plt.annotate(stats_text, xy=(0.03, 0.87), xycoords='axes fraction',
             verticalalignment='top', bbox=dict(boxstyle="round", alpha=0.5),
             fontsize=10)

# Final plot adjustments
plt.xlabel("GEOS-Chem Predictions (ug/m3)")
plt.ylabel("InMAP Predictions (ug/m3)")
plt.title("Total PM2.5 Model Comparison (not pop-wtd.)")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show the plot
plt.savefig("../outputs/totalpm_scatter_scaled.png")
