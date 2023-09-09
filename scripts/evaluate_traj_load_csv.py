import pickle
# import yaml
import csv

import numpy as np
import os
import pandas as pd

#increase number of pandas prints width
pd.set_option('display.max_columns', 50)


directory = '/home/yifei/ws/src/kr_autonomous_flight/autonomy_core/map_plan/action_planner/scripts/res'
#what do we need
#we have 10 methods and 3 environemnt index

#hanlde files from a the folder if there is no result associated with that file
# step 1: iterate through the folder
for filename in os.listdir(directory):
    if "csv" in filename[-3:]:
        print(filename)
        full_path = os.path.join(directory, filename)
        # step 2: check if there is a result file associated with that file
        # if os.path.isfile(full_path + '.txt'):
        #     print("file exists")
        #     continue
        # step 3: if there is no result file, then load the file and run the evaluation
        df = pd.read_csv(full_path)
        #divide environemnt index in column 'density_index' into bins and based on percentile, divide into three dfs
        print(df.describe())
        name_of_col = ['density_index','clutter_index','structure_index']
        index_limits = np.zeros((3,4))
        df_collect = []
        df_name_collect = []
        for i in range(3):
            index_name = name_of_col[i]
            q33 = df[index_name].quantile(1/3)
            q66 = df[index_name].quantile(2/3)
            df1 = df[df[index_name] <= q33]
            df2 = df[(df[index_name] > q33) & (df[index_name] <= q66)]
            df3 = df[df[index_name] > q66]
            # print(df1.describe())
            # print(df2.describe())
            # print(df3.describe())
            print(f"{index_name} 0-33%: Min = {df1[index_name].min()}, Max = {df1[index_name].max()}")
            print(f"{index_name} 33-66%: Min = {df2[index_name].min()}, Max = {df2[index_name].max()}")
            print(f"{index_name} 66-100%: Min = {df3[index_name].min()}, Max = {df3[index_name].max()}")
            index_limits[i,0] = df1[index_name].min()
            index_limits[i,1] = q33
            index_limits[i,2] = q66
            index_limits[i,3] = df3[index_name].max()        
        #now cut the df into 27 dfs based on the index_limits
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    df_collect.append(df[(df[name_of_col[0]] >= index_limits[0,i]) & (df[name_of_col[0]] <= index_limits[0,i+1]) &
                                          (df[name_of_col[1]] >= index_limits[1,j]) & (df[name_of_col[1]] <= index_limits[1,j+1]) & 
                                          (df[name_of_col[2]] >= index_limits[2,k]) & (df[name_of_col[2]] <= index_limits[2,k+1])])
                    df_name_collect.append(name_of_col[0]+' '+str(i) + ' '+name_of_col[1]+ ' '+str(j) + ' '+ name_of_col[2] + ' '+str(k))
        #now we have 27 dfs, we can do the evaluation
        sig_figs = 2

# Format the number with the specified number of significant figures
        for i in range(27):
            print("DF: ", i)
            
            for col_name in name_of_col:
                formatted_min = "{:.{p}g}".format(df_collect[i][col_name].min(), p=sig_figs)
                formatted_max = "{:.{p}g}".format(df_collect[i][col_name].max(), p=sig_figs)
                
            print(df_name_collect[i])
            print("number of data pts", df_collect[i][col_name].count(),"success = ", df_collect[i]["success"].sum())
            # print(df_collect[i].describe())




            # csv_reader = csv.reader(file)
            # for row in csv_reader:
            #     print(row)
            #     break
        # step 4: save the result file
        with open(full_path + '.txt', mode='w') as file:
            file.write("hello world") 
 

