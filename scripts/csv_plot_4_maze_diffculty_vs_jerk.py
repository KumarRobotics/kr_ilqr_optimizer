# Import the necessary plotting libraries again
import matplotlib.pyplot as plt
import seaborn as sns

# Import the necessary libraries again
import pandas as pd

directory = '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_Maze_single_line_09-09_14-50-20.csv'

# Define the file path again
file_path = directory

# Re-load the original CSV file
df_original = pd.read_csv(file_path)

# Extract the part of the string that appears after 'maze' and before the first '_' in the 'map_filename' column
df_original['maze_description'] = df_original['map_filename'].str.extract('maze(\d+)_')[0]

# Display the first few rows to verify the new column
df_original.head()

# Define the columns to describe again
columns_to_describe = ['compute_time_frontend(ms)', 'compute_time_backend(ms)', 'traj_time(s)', 'traj_length(m)', 'traj_jerk']
df_original = df_original.sort_values(by=['maze_description'])
# Group by 'maze_description', 'planner_frontend', and 'planner_backend' and compute descriptive statistics
grouped_descriptive_stats = df_original.groupby(['maze_description', 'planner_frontend', 'planner_backend'])[columns_to_describe].describe()

# Show the first few rows to verify
grouped_descriptive_stats.head()

# Initialize a larger matplotlib figure for the single boxplot with counts
plt.figure(figsize=(20, 10))

# Create the boxplot
# df_melted = pd.melt(df_original, id_vars=['planner_frontend', 'planner_backend'], value_vars=['planner_frontend', 'planner_backend'])

ax = sns.boxplot(data=df_original, x='maze_description', y='traj_jerk', hue=df_original[['planner_frontend', 'planner_backend']].apply(tuple, axis=1), showfliers=False)

# Add the counts above the boxplots
for i, maze_desc in enumerate(df_original['maze_description'].unique()):
    # Filter data for the current maze_description
    filtered_data = df_original[df_original['maze_description'] == maze_desc]
    
    # Calculate counts for each planner_frontend
    counts = filtered_data['planner_frontend'].value_counts().sort_index()
    
    # Add counts as text above the boxplots
    for j, (planner_frontend, count) in enumerate(counts.items()):
        ax.text(i + j * 0.2 - 0.3, ax.get_ylim()[1] * 0.8, f'N={count}', ha='center')

# Set title and legend
plt.title('Boxplot of Trajectory Jerk Grouped by Maze Description and Planner Frontend (Without Outliers)')
plt.legend(title='Planner Frontend', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()