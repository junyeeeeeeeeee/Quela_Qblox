import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

# Data dictionary
data_dict = {
    "FQv1_TaAl_1217OS, AS16": [
        {"T2/T1": 0.19, "T2/T2*": 1.14, "eff_T": "nan"},
        {"T2/T1": 0.33, "T2/T2*": 0.94, "eff_T": "nan"},
        {"T2/T1": 0.82, "T2/T2*": 1.08, "eff_T": "nan"},
        {"T2/T1": 0.7, "T2/T2*": 1.24, "eff_T": "nan"}
    ],
    "RBFQ_1209OS, AS16": [{"T2/T1": 1.97, "T2/T2*": 1.08, "eff_T": 42}, {"T2/T1": 1.62, "T2/T2*": 1.03, "eff_T": 63}, {"T2/T1": 1.64, "T2/T2*": 2.16, "eff_T": 58}, {"T2/T1": 1.39, "T2/T2*": 2.17, "eff_T": 53}],
    "FQv1_WJv7_01, AS16v3": [{"T2/T1": 0.66, "T2/T2*": 1.56, "eff_T": 47}, {"T2/T1": 0.57, "T2/T2*": 1.36, "eff_T": 43}, {"T2/T1": 1.06, "T2/T2*": 1.9, "eff_T": 40}, {"T2/T1": 0.74, "T2/T2*": 2.75, "eff_T": 44}],
    "6FQ_AS217OS_050122_4, AS16v3": [{"T2/T1": 0.86, "T2/T2*": 0.96, "eff_T": "nan"}, {"T2/T1": 0.49, "T2/T2*": 1.00, "eff_T": 96}, {"T2/T1": 0.37, "T2/T2*": 0.93, "eff_T": 87}]
}

# Extract unique sample names and assign color indices
sample_names = list(data_dict.keys())
sample_indices = {name: i for i, name in enumerate(sample_names)}

# Extract data points
x_vals, y_vals, z_vals, sample_labels = [], [], [], []
for sample, measurements in data_dict.items():
    for entry in measurements:
        x_vals.append(entry["T2/T1"])
        y_vals.append(entry["T2/T2*"])
        z_vals.append(entry["eff_T"])
        sample_labels.append(sample)

# Define grid bins
x_bins = np.linspace(0, 3, 20)  # 10 bins in x
y_bins = np.linspace(0, 3, 20)  # 10 bins in y
sample_map = np.full((len(x_bins) - 1, len(y_bins) - 1), np.nan)

# Assign data to grid
for i, (x, y, z, label) in enumerate(zip(x_vals, y_vals, z_vals, sample_labels)):
    x_idx = np.searchsorted(x_bins, x, side="right") - 1
    y_idx = np.searchsorted(y_bins, y, side="right") - 1
    sample_map[x_idx, y_idx] = z  # Store eff_T values

# Create colormap with fixed colors per sample
color_map = {name: plt.cm.tab10(i / len(sample_names)) for i, name in enumerate(sample_names)}
indexed_map = np.full_like(sample_map, np.nan)

for i, (x, y, z, label) in enumerate(zip(x_vals, y_vals, z_vals, sample_labels)):
    x_idx = np.searchsorted(x_bins, x, side="right") - 1
    y_idx = np.searchsorted(y_bins, y, side="right") - 1
    indexed_map[x_idx, y_idx] = sample_indices[label]

# Use a ListedColormap to enforce correct colors
cmap = ListedColormap([color_map[name] for name in sample_names])

# Plot colormesh
fig, ax = plt.subplots(figsize=(8, 6))

mesh = ax.pcolormesh(x_bins, y_bins, indexed_map.T, cmap=cmap, edgecolors="black", shading="auto", vmin=0, vmax=len(sample_names) - 1)

# Add text labels for eff_T values
for i in range(len(x_bins) - 1):
    for j in range(len(y_bins) - 1):
        if not np.isnan(sample_map[i, j]):  
            ax.text((x_bins[i] + x_bins[i+1]) / 2, (y_bins[j] + y_bins[j+1]) / 2, 
                    f"{int(sample_map[i, j])}", ha="center", va="center", fontsize=10, 
                    color="white", fontweight="bold") # ,bbox=dict(facecolor="black", alpha=0.4, edgecolor="none")


# ax.scatter(2,1,marker="*",s=240,edgecolors='red',label="ee")
ax.vlines(x=2, ymin=0, ymax=3, linestyles="--")
ax.hlines(y=1, xmin=0, xmax=3, linestyles="--")
# Create legend with correctly matched colors
legend_patches = [mpatches.Patch(color=color_map[name], label=name) for name in sample_names]
ax.legend(handles=legend_patches, title="Sample Names", loc="upper left", bbox_to_anchor=(1, 1))

# Labels and title
ax.set_xlabel("T2/T1")
ax.set_ylabel("T2/T2*")
ax.set_title("Decoherence Relationship")

plt.tight_layout()
plt.savefig("./fff.png")
plt.close()

