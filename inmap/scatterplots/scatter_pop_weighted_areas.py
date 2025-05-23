import geopandas as gpd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

# Load data from shapefile
shapefile_path = "../intermediate/total_PM25.shp"
gdf = gpd.read_file(shapefile_path)
mortality_ra = gdf["Total_PM25"].values

# Load population data
pop_shapefile_path = "../inputs/TotalPop_regridded.shp"
pop_gdf = gpd.read_file(pop_shapefile_path)
population = pop_gdf["MortalityRa"].values  # Adjust column name if needed

# Load data from NetCDF file
netcdf_path = "$PROJECT/GC/diff.nc4"
dataset = nc.Dataset(netcdf_path)
aermass_bc = dataset.variables["PM25"][0, 0, :, :]  # Select time=0, level=0

# Flatten the NetCDF data to match the shapefile structure
aermass_bc_flat = aermass_bc.flatten()

# Ensure data lengths match
if len(mortality_ra) != len(aermass_bc_flat) or len(population) != len(mortality_ra):
    raise ValueError("Mismatch in grid cell lengths between datasets.")

# Normalize population for circle sizes
min_pop, max_pop = population.min(), population.max()
norm_pop = (population - min_pop) / (max_pop - min_pop)
marker_sizes = 10 + (norm_pop * 190)  # Scale sizes from 10 to 200

# Scatterplot with circle size representing population weight
plt.figure(figsize=(8, 8))
plt.scatter(aermass_bc_flat, mortality_ra, s=marker_sizes, alpha=0.6, edgecolors='black')

# Add y=x line
plt.plot([0.0, 1.5], [0.0, 1.5], linestyle="--", color="black", label="y = x")

# Compute population-weighted statistics
def weighted_mean(x, weights):
    return np.sum(x * weights) / np.sum(weights)

def weighted_variance(x, weights):
    mean_x = weighted_mean(x, weights)
    return np.sum(weights * (x - mean_x) ** 2) / np.sum(weights)

def weighted_r2(y_true, y_pred, weights):
    ss_total = np.sum(weights * (y_true - weighted_mean(y_true, weights)) ** 2)
    ss_residual = np.sum(weights * (y_true - y_pred) ** 2)
    return 1 - (ss_residual / ss_total)

# Compute regression slope manually (forcing intercept = 0)
slope = np.sum(population * aermass_bc_flat * mortality_ra) / np.sum(population * aermass_bc_flat ** 2)
predictions = slope * aermass_bc_flat

# Compute weighted R^2
r2_weighted = weighted_r2(mortality_ra, predictions, population)

# Compute Normalized Mean Bias (NMB) and Normalized Mean Error (NME)
nmb = np.sum(population * (mortality_ra - predictions)) / np.sum(population * mortality_ra) * 100.0
nme = np.sum(population * np.abs(mortality_ra - predictions)) / np.sum(population * mortality_ra) *100.0

# Add line of best fit
#plt.plot(aermass_bc_flat, predictions, color="red", label=f"Best Fit (Slope: {slope:.2f})")
plt.plot(aermass_bc_flat, predictions, color="red", label=f"Best Fit")
# Display statistics
plt.xlabel("GEOS-Chem Predictions (ug/m3)")
plt.ylabel("InMAP Predictions (ug/m3)")
plt.title("Total PM2.5 Model Comparison (Population-Weighted)")
plt.text(0.03, 0.80, f"Weighted R^2: {r2_weighted:.2f}\nNMB: {nmb:.1f}%\nNME: {nme:.1f}%\nSlope: {slope:.2f}", transform=plt.gca().transAxes)
plt.legend()
plt.xlim(0.0, 1.2)
plt.ylim(0.0, 1.5)
plt.savefig("../outputs/pop_weighted_scatter.pdf")
