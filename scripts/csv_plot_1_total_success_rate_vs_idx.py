import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Read the CSV file into a DataFrame
directory = '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_single_line_09-09_02-03-44_big.csv'

file_path = directory
df = pd.read_csv(file_path)

# Group by 'map_seed' and calculate the mean of the 'success' column
grouped_df = df.groupby(['map_seed', 'density_index', 'clutter_index', 'structure_index'])['success'].mean().reset_index()

# Convert the mean success rate to percentage
grouped_df['success_rate_percentage'] = grouped_df['success'] * 100

# To use a log scale, ensure that all values are greater than zero by adding a small constant (epsilon)
epsilon = 1e-6

# Apply log transformation to the original columns
grouped_df['density_index_log'] = np.log(grouped_df['density_index'] + epsilon)
grouped_df['clutter_index_log'] = np.log(grouped_df['clutter_index'] + epsilon)
grouped_df['structure_index_log'] = np.log(grouped_df['structure_index'] + epsilon)

# Extract the log-transformed data for plotting
density_index_log = grouped_df['density_index_log']
clutter_index_log = grouped_df['clutter_index_log']
structure_index_log = grouped_df['structure_index_log']
success_rate_percentage = grouped_df['success_rate_percentage']

# Create the 3D scatter plot with log-transformed axes
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(density_index_log, clutter_index_log, structure_index_log, c=success_rate_percentage, cmap='viridis', s=10)

# Add color bar and labels
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Success Rate (%)')
ax.set_xlabel('Density Index (Log Scale)')
ax.set_ylabel('Clutter Index (Log Scale)')
ax.set_zlabel('Structure Index (Log Scale)')
ax.set_title('3D Scatter Plot of Success Rate vs Log-transformed Density, Clutter, and Structure Indices')

plt.savefig('fig1.png')
