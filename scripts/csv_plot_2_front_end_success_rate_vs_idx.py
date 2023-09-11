
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Load the CSV file into a DataFrame
directory = '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_Structure_single_line_09-09_02-03-44.csv'

file_path = directory
df = pd.read_csv(file_path)

# Group the data by 'planner_frontend' and 'planner_backend'
grouped_data = df.groupby(['planner_frontend', 'planner_backend'])

# Initialize a new figure to hold the 3D scatter plots with natural log scale and success rates
fig = plt.figure(figsize=(8, 20))

# Iterate through each group and plot the 3D scatter plot with log scale and success rates
for i, ((frontend, backend), group) in enumerate(grouped_data):
    ax = fig.add_subplot(5, 2, i+1, projection='3d')
    ax.set_title(f'Frontend: {frontend}, Backend: {backend}')
    ax.set_xlabel('ln(Density Index)')
    ax.set_ylabel('ln(Clutter Index)')
    ax.set_zlabel('ln(Structure Index)')
    
    # Calculate the success rate for this group
    total_trials = len(group)
    successful_trials = len(group[group['success'] == True])
    success_rate = (successful_trials / total_trials) * 100 if total_trials > 0 else 0
    
    # Add success rate as text in the plot
    ax.text2D(0.05, 0.95, f'Success Rate: {success_rate:.2f}%', transform=ax.transAxes, color='r')
    
    # Apply natural log transformation to the indices, adding a small constant to avoid log(0)
    small_constant = 1e-10
    
    # Plot points where success is True
    success_group = group[group['success'] == True]
    ax.scatter(np.log(success_group['density_index'] + small_constant), 
               np.log(success_group['clutter_index'] + small_constant), 
               np.log(success_group['structure_index'] + small_constant), c='g', label='Success', s=5)
    
    # Plot points where success is False
    failure_group = group[group['success'] == False]
    ax.scatter(np.log(failure_group['density_index'] + small_constant), 
               np.log(failure_group['clutter_index'] + small_constant), 
               np.log(failure_group['structure_index'] + small_constant), c='r', label='Failure', s=5)
    
    ax.legend()
    ax.set_xlim([-6, 0])

plt.tight_layout()
plt.savefig('fig2.png')