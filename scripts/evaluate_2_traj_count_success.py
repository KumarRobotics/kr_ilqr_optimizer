import pandas as pd
pd.set_option('display.max_columns', 50)

##Save prompt here!

# Step 1: Load the data from the provided CSV file
c# Re-load the data from the re-uploaded CSV file
file_path = directory
df = pd.read_csv(file_path)

# Re-group the data by 'planner_frontend' and 'planner_backend'
grouped_df = df.groupby(['planner_frontend', 'planner_backend'])

# Initialize an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['count','planner_frontend', 'planner_backend', 
                                  'success_rate', 'collision_rate', 'server failure rate', 
                                  'avg_frontend_time_success', 'avg_poly_time_success','avg_backend_time_success',
                                  'avg_frontend_time_not_success', 'avg_poly_time_not_success', 'avg_backend_time_not_success',
                                    'dist to goal(m)'])

# Iterate through each group to calculate various metrics
for (planner_frontend, planner_backend), group in grouped_df:
    total_count = len(group)
    success_count = group['success'].sum()
    collision_count = group['collision_status'].sum()
    count_minus_one = (group['success_detail'] == -1).sum()
    
    success_rate = (success_count / total_count) * 100
    collision_rate = (collision_count / total_count) * 100
    server_failure_rate = (count_minus_one / total_count) * 100

    # Average time for successful and unsuccessful trials
    avg_frontend_time_success = group[group['success'] == True]['compute_time_frontend(ms)'].mean()
    avg_backend_time_success = group[group['success'] == True]['compute_time_backend(ms)'].mean()
    avg_poly_time_success = group[group['success'] == True]['compute_time_poly(ms)'].mean()
    avg_poly_time_success = group[group['success'] == True]['compute_time_poly(ms)'].mean()
    dist_to_goal = group[group['success'] == True]['dist_to_goal(m)'].mean()

    avg_frontend_time_not_success = group[group['success'] == False]['compute_time_frontend(ms)'].mean()
    avg_backend_time_not_success = group[group['success'] == False]['compute_time_backend(ms)'].mean()
    avg_poly_time_not_success = group[group['success'] == False]['compute_time_poly(ms)'].mean()


    
    result_df = result_df.append({
        'planner_frontend': planner_frontend,
        'planner_backend': planner_backend,
        'success_rate': success_rate,
        'collision_rate': collision_rate,
        'avg_frontend_time_success': avg_frontend_time_success,
        'avg_poly_time_success': avg_poly_time_success,
        'avg_backend_time_success': avg_backend_time_success,
        'avg_frontend_time_not_success': avg_frontend_time_not_success,
        'avg_poly_time_not_success': avg_poly_time_not_success,
        'avg_backend_time_not_success': avg_backend_time_not_success,
        'server failure rate': server_failure_rate,
        'count': total_count,
        'dist to goal(m)': dist_to_goal
    }, ignore_index=True)

# Sort the result DataFrame
result_df_sorted_with_all_data = result_df.sort_values(by=['planner_frontend', 'planner_backend'])
print(result_df_sorted_with_all_data)

#save the result_df_sorted_with_all_data to csv
result_df_sorted_with_all_data.to_csv('/home/yifei/ws/src/kr_autonomous_flight/kr_ilqr_optimizer/scripts/res/ECI_single_line_09-09_02-03-44_big.csv.csv', index=False)
