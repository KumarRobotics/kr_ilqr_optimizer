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

#convert success collumn to bool
df_original['success'] = df_original['success'].astype(bool)
df_original['collision_status'] = df_original['collision_status'].astype(bool)

#filter success = true and collision = false
df_original = df_original[(df_original['success'] == True) & (df_original['collision_status'] == False)]

# Extract the part of the string that appears after 'maze' and before the first '_' in the 'map_filename' column
df_original['maze_description'] = df_original['map_filename'].str.extract('maze(\d+)_')[0].astype(int)

# Display the first few rows to verify the new column
df_original.head()

# Define the columns to describe again
columns_to_describe =['compute_time_frontend(ms)', 'compute_time_backend(ms)', 'traj_time(s)', 'traj_length(m)', 'traj_jerk']
df_original = df_original.sort_values(by=['maze_description', 'planner_frontend', 'planner_backend'])
# Group by 'maze_description', 'planner_frontend', and 'planner_backend' and compute descriptive statistics
grouped_descriptive_stats = df_original.groupby(['maze_description', 'planner_frontend', 'planner_backend'])[columns_to_describe].describe().to_csv("des.csv")

# Show the first few rows to verify
print(grouped_descriptive_stats)

# Initialize a larger matplotlib figure for the single boxplot with counts
plt.figure(figsize=(20, 10))

# Create the boxplot
# df_melted = pd.melt(df_original, id_vars=['planner_frontend', 'planner_backend'], value_vars=['planner_frontend', 'planner_backend'])

tupled_col = df_original[['planner_frontend', 'planner_backend']].apply(tuple, axis=1)


ax = sns.boxplot(data=df_original, x='maze_description', y='traj_time(s)', hue=tupled_col, showfliers=False)



# Set title and legend
plt.title('Boxplot of Trajectory Jerk Grouped by Maze Description and Planner Frontend (Without Outliers)')
plt.legend(title='Planner Frontend', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()