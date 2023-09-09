# Importing necessary libraries for plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
import numpy as np
directory = '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_single_line_09-09_02-03-44_big.csv'

file_path = directory
df = pd.read_csv(file_path)
epsilon = 1e-6
df['density_index'] = np.log(df['density_index'] + epsilon)
df['clutter_index'] = np.log(df['clutter_index'] + epsilon)
df['structure_index'] = np.log(df['structure_index'] + epsilon)
fig = plt.figure(figsize=(18, 12))

#
# Filter the rows where 'planner_backend' is 'gcopter'
filtered_rows = df[df['planner_backend'] == 'gcopter']
filtered_rows2 = df[df['planner_backend'] == 'iLQR(Altro)']
# Get unique front-end planners
unique_frontend_planners = filtered_rows['planner_frontend'].unique()
unique_frontend_planners2 = filtered_rows['planner_frontend'].unique()

grouped_df = df.groupby(['planner_frontend', 'planner_backend'])
success_rate_list = []
for (planner_frontend, planner_backend), group in grouped_df:
    total_count = len(group)
    success_count = group['success'].sum()
    collision_count = group['collision_status'].sum()
    count_minus_one = (group['success_detail'] == -1).sum()
    
    success_rate = (success_count / total_count) * 100
    success_rate_list.append(success_rate)
# Initialize the plot

# Loop over each unique front-end planner to create a separate subplot
for i, planner in enumerate(unique_frontend_planners):
    ax = fig.add_subplot(2, 5, i+1, projection='3d')
    
    # Filter data for the current front-end planner
    planner_data = filtered_rows[filtered_rows['planner_frontend'] == planner]
    
    # Scatter plot
    scatter = ax.scatter(planner_data['density_index'], 
                         planner_data['clutter_index'], 
                         planner_data['structure_index'], 
                         c=planner_data['success'].apply(lambda x: 'g' if x else 'r'), 
                         marker='o',
                         label=['Success', 'Failure'],  s=1)
    
    # Labels and title
    ax.set_xlabel('Density Index')
    ax.set_ylabel('Clutter Index')
    ax.set_zlabel('Structure Index')
    ax.set_title(f'Planner: GCOPTER + {planner}')

    # add text in red showing the success rate
    ax.text2D(0.05, 0.95, f"Success Rate: {success_rate_list[i]:.1f}%", transform=ax.transAxes, color='red')
    

#do the same for altro
for i, planner in enumerate(unique_frontend_planners2):
    ax = fig.add_subplot(2, 5, i+6, projection='3d')
    
    # Filter data for the current front-end planner
    planner_data = filtered_rows2[filtered_rows2['planner_frontend'] == planner]
    
    # Scatter plot
    scatter = ax.scatter(planner_data['density_index'], 
                         planner_data['clutter_index'], 
                         planner_data['structure_index'], 
                         c=planner_data['success'].apply(lambda x: 'g' if x else 'r'), 
                         marker='o',
                         label=['Success', 'Failure'], s=1)
    
    # Labels and title
    ax.set_xlabel('Density Index')
    ax.set_ylabel('Clutter Index')
    ax.set_zlabel('Structure Index')
    ax.set_title(f'Planner: iLQR + {planner}')
    
    ax.text2D(0.05, 0.95, f"Success Rate: {success_rate_list[i+5]:.1f}%", transform=ax.transAxes, color='red')


plt.tight_layout()
plt.savefig('fig2.png')




