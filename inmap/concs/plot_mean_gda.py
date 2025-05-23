import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.patches import FancyArrow
from matplotlib_scalebar.scalebar import ScaleBar

# Load the shapefile
gdf = gpd.read_file("../inputs/clipped_mean.shp")

# Get value range
#vmin = gdf["TotalPM25"].min()
#vmax = gdf["TotalPM25"].max()
vmin = 0.0
vmax = 18.6

# Projection
proj = ccrs.PlateCarree()

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': proj})

# Plot PM2.5 data
gdf.plot(column="TotalPM25", cmap="RdYlGn", vmin=vmin, vmax=vmax,
         ax=ax, transform=proj)

# Set extent to Greater Dhaka Area (adjusted manually to tightly frame it)
ax.set_extent([90.0, 90.99, 23.35, 24.35], crs=proj)

# Gridlines
gl = ax.gridlines(draw_labels=True, linewidth=0.0, color='grey', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# Add Dhaka label
ax.plot(90.4125, 23.8103, marker='o', color='black', markersize=3, transform=proj)
#ax.text(90.4125 + 0.03, 23.8103, 'Dhaka', transform=proj, fontsize=9, verticalalignment='center')

# North arrow (top right)
arrow = FancyArrow(x=0.95, y=0.9, dx=0, dy=0.05, width=0.005,
                   transform=ax.transAxes, color='black')
ax.add_patch(arrow)
ax.text(0.95, 0.97, 'N', transform=ax.transAxes, ha='center', va='bottom', fontsize=10, weight='bold')

# Scalebar (~1 km = 0.009 degrees)
scalebar = ScaleBar(dx=111000, units='m', length_fraction=0.1, location='lower left',
                    box_alpha=0, scale_loc='bottom')
ax.add_artist(scalebar)

# Load Greater Dhaka Area district boundaries
gda = gpd.read_file("../inputs/GDA.shp")

# Plot Greater Dhaka districts
gda.boundary.plot(ax=ax, color="black", linewidth=1)

# Colorbar
sm = plt.cm.ScalarMappable(cmap="RdYlGn", norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label("PM2.5 concentrations (µg/m³)")

plt.tight_layout()
plt.savefig("../outputs/maps/plot_mean_gda_clean.png", dpi=300, bbox_inches='tight')
plt.show()
