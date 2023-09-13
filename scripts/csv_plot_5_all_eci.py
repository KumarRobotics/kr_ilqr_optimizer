import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
# File paths
file_paths = [
    '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_Real_single_line_09-09_18-14-35.csv',
    '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_Maze_single_line_09-09_14-50-20.csv',
    '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_Structure_single_line_09-09_02-03-44.csv'
]

# Initialize 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Colors for each file
colors = ['r', 'g', 'b']

# Loop through each file and plot
for i, file_path in enumerate(file_paths):
    df = pd.read_csv(file_path)
    ax.scatter(np.log(df['density_index']), np.log(df['clutter_index']), np.log(df['structure_index']), c=colors[i], label=file_path.split('/')[-1], alpha=0.2)

# Labels and title
ax.set_xlabel('Density Index')
ax.set_ylabel('Clutter Index')
ax.set_zlabel('Structure Index')
ax.set_title('3D Scatter Plot of Environment Indices')
ax.legend()

plt.tight_layout()
plt.savefig('fig5.png')