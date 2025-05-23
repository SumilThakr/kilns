import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = "../outputs/popwtd/for_stacked_bar_nosoa.csv"
df = pd.read_csv(file_path)

# Rename columns for clarity
column_mapping = {
    "GC_avg": "GEOS-Chem",
    "GC_popwtd": "GEOS-Chem (pop-wtd.)",
    "InMAP_avg": "InMAP",
    "InMAP_popwtd": "InMAP (pop-wtd.)"
}
df = df.rename(columns=column_mapping)

# Set first column as index (pollutant types)
df.set_index(df.columns[0], inplace=True)

# Transpose so bars are models
df = df.T

# Plot
fig, ax = plt.subplots(figsize=(9, 6))
df.plot(kind='bar', stacked=True, ax=ax, edgecolor='none')

# Axes labels & title
ax.set_ylabel("Avg. pollutant concentrations (μg/m³)", fontsize=12)
ax.set_xlabel("Model", fontsize=12)
ax.set_title("Pollutant Concentrations by Model", fontsize=14, pad=15)

# Move legend to the right
ax.legend(title="Pollutant Type", bbox_to_anchor=(1.05, 1), loc='upper left')

# Remove gap between bars and x-axis
ax.margins(y=0)

# Add horizontal dashed grid lines behind bars
ax.yaxis.grid(True, linestyle='--', alpha=0.7)
ax.set_axisbelow(True)

# Set y-axis limit to 0.6
ax.set_ylim(0, 0.6)

# Rotate x-axis labels nicely
plt.xticks(rotation=0)

# Clean layout
plt.tight_layout()

# Save figure
plt.savefig("../outputs/popwtd/pop_wtd_concentrations.pdf", dpi=300)

plt.show()
