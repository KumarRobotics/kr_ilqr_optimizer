import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

# Load the entire CSV file into a DataFrame
directory = '/home/yifei/ws/src/kr_autonomous_flight/autonomy_core/map_plan/action_planner/scripts/res/ECI_single_line_09-08_14-58-11.csv'

file_path = directory
df = pd.read_csv(file_path)

# Filter the DataFrame to keep only the rows where 'success' is True
df_success = df[df['success'] == True]

# Count the number of unique values in the 'success_detail' column
unique_success_detail_values = df_success['success_detail'].nunique()
print(f"Unique values in success_detail: {unique_success_detail_values}")

# Create a new DataFrame to hold the rows with the shortest 'traj_time(s)' for each combination of
# 'density_index', 'clutter_index', and 'structure_index' and 'planner_frontend'
grouped = df_success.groupby(['density_index', 'clutter_index', 'structure_index', 'planner_frontend'])
df_min_time = grouped.apply(lambda x: x[x['traj_time(s)'] == x['traj_time(s)'].min()])

# Create the 3D scatter plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Color map for different 'planner_frontend'
colors = sns.color_palette('husl', n_colors=df_min_time['planner_frontend'].nunique())
color_map = dict(zip(df_min_time['planner_frontend'].unique(), colors))

for planner, color in color_map.items():
    subset = df_min_time[df_min_time['planner_frontend'] == planner]
    print(subset.count())
    ax.scatter(subset['density_index'], subset['clutter_index'], subset['structure_index'],
               c=[color]*len(subset), label=planner, s=20, alpha=1)

# Labels and legend
ax.set_xlabel('Density Index')
ax.set_ylabel('Clutter Index')
ax.set_zlabel('Structure Index')
ax.legend(title='Planner Frontend', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_title('Shortest Trajectory Time by Planner Frontend')

plt.tight_layout()
plt.show()