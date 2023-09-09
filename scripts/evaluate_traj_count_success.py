import pandas as pd
pd.set_option('display.max_columns', 50)

# Step 1: Load the data from the provided CSV file
directory = '/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_single_line_09-08_14-58-11_big.csv'

file_path = directory
df = pd.read_csv(file_path)

# Step 2: Group the data by 'planner_frontend' and 'planner_backend'
grouped_df = df.groupby(['planner_frontend', 'planner_backend'])

# Step 3: Initialize an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['planner_frontend', 'planner_backend', 'success_rate', 'collision_rate'])

# Step 4: Iterate through each group to calculate the success rate and collision rate
for (planner_frontend, planner_backend), group in grouped_df:
    total_count = len(group)
    success_count = group['success'].sum()
    collision_count = group['collision_status'].sum()
    
    success_rate = (success_count / total_count) * 100  # in percentage
    collision_rate = (collision_count / total_count) * 100  # in percentage
    
    result_df = result_df.append({
        'planner_frontend': planner_frontend,
        'planner_backend': planner_backend,
        'success_rate': success_rate,
        'collision_rate': collision_rate
    }, ignore_index=True)

# Step 5: Sort the result by 'planner_frontend' in alphabetical order
result_df_sorted = result_df.sort_values(by=['planner_frontend', 'planner_backend'])

# Show the sorted result
# print(result_df_sorted)

# Initialize additional columns for average front-end and back-end planning time for successful trials
result_df['avg_frontend_time_success'] = 0.0
result_df['avg_backend_time_success'] = 0.0

# Iterate through each group to calculate the average front-end and back-end planning time for successful trials
for (planner_frontend, planner_backend), group in grouped_df:
    # Filter out only the successful trials
    successful_group = group[group['success'] == True]
    
    # Calculate the average front-end and back-end planning time for successful trials
    avg_frontend_time_success = successful_group['compute_time_frontend(ms)'].mean()
    avg_backend_time_success = successful_group['compute_time_backend'].mean()
    
    # Update the result DataFrame with the calculated values
    result_df.loc[(result_df['planner_frontend'] == planner_frontend) & 
                  (result_df['planner_backend'] == planner_backend), 
                  ['avg_frontend_time_success', 'avg_backend_time_success']] = avg_frontend_time_success, avg_backend_time_success

# Sort the result by 'planner_frontend' in alphabetical order again
result_df_sorted_with_time = result_df.sort_values(by=['planner_frontend', 'planner_backend'])

# Initialize additional columns for average front-end and back-end planning time for unsuccessful trials
result_df['avg_frontend_time_not_success'] = 0.0
result_df['avg_backend_time_not_success'] = 0.0

# Iterate through each group to calculate the average front-end and back-end planning time for unsuccessful trials
for (planner_frontend, planner_backend), group in grouped_df:
    # Filter out only the unsuccessful trials
    unsuccessful_group = group[group['success'] == False]
    
    # Calculate the average front-end and back-end planning time for unsuccessful trials
    avg_frontend_time_not_success = unsuccessful_group['compute_time_frontend(ms)'].mean()
    avg_backend_time_not_success = unsuccessful_group['compute_time_backend'].mean()
    
    # Update the result DataFrame with the calculated values
    result_df.loc[(result_df['planner_frontend'] == planner_frontend) & 
                  (result_df['planner_backend'] == planner_backend), 
                  ['avg_frontend_time_not_success', 'avg_backend_time_not_success']] = avg_frontend_time_not_success, avg_backend_time_not_success


print(result_df)