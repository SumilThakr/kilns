import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Load the shapefile
shapefile_path = "../inputs/clipped_max.shp"
gdf = gpd.read_file(shapefile_path)

# Set custom vmin and vmax
vmin = gdf["TotalPM25"].min()
vmax = gdf["TotalPM25"].max()

# Load Bangladesh outline
bangladesh_outline_path = "$PROJECT/GC/mask/country_mask.shp"
bangladesh_outline = gpd.read_file(bangladesh_outline_path)

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

# Plot data
gdf.boundary.plot(ax=ax, color="none")
gdf.plot(
    column="TotalPM25",
    cmap="RdYlGn",
    vmin=vmin,
    vmax=vmax,
    legend=True,
    ax=ax,
    legend_kwds={"label": "Difference in PM2.5 concentrations (µg/m³)"}
)

# Plot country boundary
bangladesh_outline.boundary.plot(ax=ax, color="black", linewidth=1)

# Clean layout
ax.set_aspect('equal')
plt.axis('off')
plt.tight_layout()

# Save output
plt.savefig("../outputs/maps/plot_max_emis_concentrations.png", dpi=300, bbox_inches='tight')
plt.show()
