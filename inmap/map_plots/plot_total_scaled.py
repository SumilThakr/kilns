import geopandas as gpd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the input shapefiles and their scaling factors

# Linear bias correction:

# Quantile scaling:
shapefile_info = {
    "../inputs/PNO3_regridded.shp":   0.1143236952329479 / (1.1 * (1.0 + 18.04/62.0049)), # The 1.1 refers to the scaling in the GC PM2.5 eqn, and the other factor accounts for pNH4
    "../inputs/PrimPM_regridded.shp": 0.12778081510958195,
    "../inputs/PSO4_regridded.shp":   3.319161231867947 / (1.1 * (1.0 + 2 * 18.04/96.06)),
    "../inputs/SOA_regridded.shp":    1.0 / 1.05,
}

# Read and bias-correct each shapefile
corrected_values = []
gdf_template = None  # To store the first GeoDataFrame for metadata

for shapefile, scaling_factor in shapefile_info.items():
    gdf = gpd.read_file(shapefile)
    
    if gdf_template is None:
        gdf_template = gdf.copy()  # Store a copy to use for output metadata
    
    # Identify the first numeric column (excluding geometry)
    numeric_columns = gdf.select_dtypes(include=[np.number]).columns
    if len(numeric_columns) == 0:
        raise ValueError(f"No numeric columns found in {shapefile}")

    column_name = numeric_columns[0]  # Get the first numerical column
    corrected_values.append(gdf[column_name] / scaling_factor)  # Apply bias correction

# Sum the corrected values
gdf_template["Total_PM25_Corrected"] = sum(corrected_values)

# Save the new shapefile
output_shapefile = "../intermediate/total_PM25_with_bias_correction.shp"
gdf_template.to_file(output_shapefile)

# Load the NetCDF file
netcdf_path = "$PROJECT/GC/diff.nc4"
dataset = nc.Dataset(netcdf_path)

# Extract "PM25" (taking the first time and level)
pm25_netCDF = dataset.variables["PM25"][0, 0, :, :]

# Extract lat/lon data
lat = dataset.variables["lat"][:]
lon = dataset.variables["lon"][:]

# Create a meshgrid for NetCDF plotting
lon_mesh, lat_mesh = np.meshgrid(lon, lat)

# Set the common color scale for comparison
vmin = min(gdf_template["Total_PM25_Corrected"].min(), pm25_netCDF.min())
vmax = max(gdf_template["Total_PM25_Corrected"].max(), pm25_netCDF.max())

# Create the figure with equal aspect ratio
fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

# Plot the bias-corrected total PM2.5 from the shapefiles
gdf_template.boundary.plot(ax=axes[0], color="none")  # Removes black border
gdf_template.plot(column="Total_PM25_Corrected", cmap="RdYlGn", vmin=vmin, vmax=vmax, legend=True,
                  ax=axes[0], legend_kwds={"label": "Bias-Corrected PM2.5 (Shapefile)"})
axes[0].set_title("Total PM2.5 (Bias-Corrected Shapefile)")
axes[0].set_xticks([])
axes[0].set_yticks([])
axes[0].set_frame_on(False)

# Plot the NetCDF PM2.5 data
im = axes[1].pcolormesh(lon_mesh, lat_mesh, pm25_netCDF, cmap="RdYlGn", shading="auto",
                        vmin=vmin, vmax=vmax)
fig.colorbar(im, ax=axes[1], label="PM2.5 (NetCDF, ug/m³)")
axes[1].set_title("NetCDF PM2.5")
axes[1].set_xticks([])
axes[1].set_yticks([])
axes[1].set_aspect((lon.max() - lon.min()) / (lat.max() - lat.min()))  # Maintain aspect ratio

# Load Bangladesh outline
bangladesh_outline_path = "$PROJECT/GC/mask/country_mask.shp"
bangladesh_outline = gpd.read_file(bangladesh_outline_path)

# Plot Bangladesh outline
bangladesh_outline.boundary.plot(ax=axes[0], color="black", linewidth=1)
bangladesh_outline.boundary.plot(ax=axes[1], color="black", linewidth=1)

# Ensure both plots have the same aspect ratio
axes[0].set_aspect(1)
axes[1].set_aspect(1)

# Show the plot
plt.savefig("../outputs/maps/total_pm_scaled.png")
