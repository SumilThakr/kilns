import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Load the shapefile
shapefile_path = "../inputs/clipped_difference_2019_to_mean.shp"
gdf = gpd.read_file(shapefile_path)

# Set custom vmin and vmax
vmin = -13.05
vmax = 0.27

# Load Bangladesh outline
bangladesh_outline_path = "$PROJECT/GC/mask/country_mask.shp"
bangladesh_outline = gpd.read_file(bangladesh_outline_path)

# Define custom colormap
colors = [
    (0.0, "#000000"),   # Very dark blue/black at vmin (-13.05)
    ((-2 - vmin) / (vmax - vmin), "blue"),  # Blue at -2
    ((0 - vmin) / (vmax - vmin), "white"),  # White at 0
    ((vmax - vmin) / (vmax - vmin), "red"),  # Red at vmax (0.27)
]

custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

# Plot data
gdf.boundary.plot(ax=ax, color="none")
gdf.plot(
    column="TotalPM25",
    cmap=custom_cmap,
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
plt.savefig("../outputs/maps/plot_difference_custom_colormap.png", dpi=300, bbox_inches='tight')
plt.show()
