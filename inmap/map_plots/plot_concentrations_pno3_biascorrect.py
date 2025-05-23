import geopandas as gpd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Load the shapefile
shapefile_path = "../inputs/PNO3_regridded.shp"
gdf = gpd.read_file(shapefile_path)

# Extract "MortalityRa" values and apply bias correction
bias_correction_factor = 0.1143236952329479
gdf["MortalityRa_corrected"] = gdf["MortalityRa"] / bias_correction_factor

# Load the NetCDF file
netcdf_path = "$PROJECT/GC/diff.nc4"
dataset = nc.Dataset(netcdf_path)

# Extract "AerMassBC" (taking the first time and level)
aermass_bc = dataset.variables["AerMassNIT"][0, 0, :, :]

# Extract lat/lon data
lat = dataset.variables["lat"][:]
lon = dataset.variables["lon"][:]

# Create a meshgrid for NetCDF plotting
lon_mesh, lat_mesh = np.meshgrid(lon, lat)

# Set the common color scale
vmin = min(gdf["MortalityRa_corrected"].min(), aermass_bc.min())
vmax = max(gdf["MortalityRa_corrected"].max(), aermass_bc.max())

# Create the figure with equal aspect ratio
fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

# Plot the shapefile data
gdf.boundary.plot(ax=axes[0], color="none")  # Removes black border
gdf.plot(column="MortalityRa_corrected", cmap="RdYlGn", vmin=vmin, vmax=vmax, legend=True,
         ax=axes[0], legend_kwds={"label": "pNO3 concentrations (ug/m3)"})
axes[0].set_title("InMAP (regridded)")
axes[0].set_xticks([])
axes[0].set_yticks([])
axes[0].set_frame_on(False)

# Plot the NetCDF data ensuring equal aspect ratio
im = axes[1].pcolormesh(lon_mesh, lat_mesh, aermass_bc, cmap="RdYlGn", shading="auto",
                        vmin=vmin, vmax=vmax)
fig.colorbar(im, ax=axes[1], label="pNO3 concentrations (ug/m3)")
axes[1].set_title("GEOS-Chem")
axes[1].set_xticks([])
axes[1].set_yticks([])
axes[1].set_aspect((lon.max() - lon.min()) / (lat.max() - lat.min()))  # Ensure aspect ratio

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
plt.savefig("../outputs/maps/plot_pno3.png")
