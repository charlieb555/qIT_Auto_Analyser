from re import X
import pandas as pd
import matplotlib as plt
import matplotlib.pyplot as plt
import csv
from statistics import mean
import numpy as np
import os
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import sys
import inquirer
import statistics

def sort_data(df_results, raw_files, file_location):
    print('Data sorting...')
    df_z_prime = pd.DataFrame(columns=['Time', 'Z Prime'])

    controls_not_inc = dict()
    list_removals = list()
    for file in raw_files:
        if (file == 'removed_wells.csv'):
            df_removed = pd.read_csv(file_location + '/' + file)

            list_removals = df_removed['To Remove'].to_list()

            break

    for file in raw_files:
        try:


            file_parts = file.split()
                
            try:
                if (file_parts[1] == 'mins.csv' ) or (file_parts[1] == 'min.csv' ):
                    try:
                        df = pd.read_csv(file_location + '/' + file)

                    except:
                        
                        print('failed ' + file)
                        continue
                else:
                    continue
            except:
                continue


            
            name_split = file.split()
            time = name_split[0]


            if (len(df) == 46):
                df_2 = df.iloc[11:27]

            if (len(df) == 47):
                df_2 = df.iloc[12:28]

            if (len(df) == 48):
                df_2 = df.iloc[13:29]

            new_header = ['Rows', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']

            df_2 = df_2.drop(columns=['Unnamed: 25', 'Unnamed: 26'])
        #Possibly throwing up issue...
            #df_2 = df_2.set_axis(new_header, axis=1, inplace=False)
            df_2 = df_2.set_axis(new_header, axis=1)
            

            df_2 = df_2.set_index('Rows')
            df_2 = df_2.astype(int)
            #df_neg_cont = df_2.filter(['1', '2'])



            df_neg_cont = df_2.filter(['1','2'])


    #Removes any wells that are assinged in a 'removed_wells.csv' file:

            for well in list_removals:
            
                well_row = well[0]
                well_column = well[1:]


                df_2.loc[well_row, well_column] = np.nan



    #Finds average for negative and positive controls (rows 1,2 and 22,23).
            df_neg_cont = df_neg_cont.astype(int)

        #Removes any values that are outside of 3 standard deviations of the median.

            median_neg_cont = np.nanmedian(df_neg_cont)
            mean_neg_cont = np.nanmean(df_neg_cont)


            std_neg_cont = (np.nanstd(df_neg_cont))

            upper_lim_neg = median_neg_cont + (3 * std_neg_cont)
            lower_lim_neg = median_neg_cont - (3 * std_neg_cont)

            col_1_all = df_neg_cont['1'].tolist()
            col_2_all = df_neg_cont['2'].tolist()

            col_1_included = df_neg_cont['1'][(df_neg_cont['1'] > lower_lim_neg) & (df_neg_cont['1'] < upper_lim_neg)]
            col_1_included = col_1_included.to_list()

            col_2_included = df_neg_cont['2'][(df_neg_cont['2'] > lower_lim_neg) & (df_neg_cont['2'] < upper_lim_neg)]
            col_2_included = col_2_included.to_list()

            neg_list_included = col_1_included + col_2_included



            avg_neg_cont_final = statistics.mean(neg_list_included)

            df_2 = df_2.sub(avg_neg_cont_final)

            #Again removes any values that are outside of 3 standard deviations of the median for the positive control

            df_pos_cont = df_2.filter(['24','23'])



            median_pos_cont = np.nanmedian(df_pos_cont)
            mean_pos_cont = np.nanmean(df_pos_cont)
            std_pos_cont = (np.nanstd(df_pos_cont))
            upper_lim_pos = median_pos_cont + (3 * std_pos_cont)
            lower_lim_pos = median_pos_cont - (3 * std_pos_cont)

            #ONLY FOR MARINA DATA

            #lower_lim_pos = 100000

            col_23_included = df_pos_cont['23'][(df_pos_cont['23'] > lower_lim_pos) & (df_pos_cont['23'] < upper_lim_pos)]
            col_23_included = col_23_included.to_list()

            col_24_included = df_pos_cont['24'][(df_pos_cont['24'] > lower_lim_pos) & (df_pos_cont['24'] < upper_lim_pos)]
            col_24_included = col_24_included.to_list()

            cols_included = col_1_included + col_2_included + col_23_included + col_24_included


            #Now works out how many have been excluded.

            col_23_all = df_pos_cont['23'].tolist()
            col_24_all = df_pos_cont['24'].tolist()
            all_controls = col_1_all + col_2_all + col_23_all + col_24_all


            num_cols_excluded = len(all_controls) - len(cols_included)

            controls_not_inc[time] = num_cols_excluded

            pos_list_included = col_23_included + col_24_included


            
            avg_pos_cont_final = statistics.mean(pos_list_included)


            df_2 = df_2.divide(avg_pos_cont_final)

            z_prime = 1 - (((3*std_pos_cont) + (3*std_neg_cont)) / (mean_pos_cont - mean_neg_cont))

            #NEW WAY OF APPENDING TO DF

            #Gets length of DF and sets index number to 0 if first entry.
            
            if (len(df_z_prime) == 0):
                index_number = 0
            else:
                index_number = len(df_z_prime)

            df_z_prime.loc[index_number] = [time, z_prime]



        #This normalises the fluorescence measurements relative to the positive and negative controls for each compound at each timepoint and puts them into a dataframe.

            n = 0

            row = 0
            

            while row <16:
                while n < 24:
                    column = df_2.iloc[row]
                    
                    fluorescence = column.iloc[n]
                    n = n + 1


                    index_key = column.name + str(n)


                    #NEW WAY OF APPENDING TO DF

                    #Gets length of DF and sets index number to 0 if first entry.

                    if (len(df_results) == 0):

                        index_number_results = 0

                    else:

                        index_number_results = len(df_results)

                    #Adds results using loc.

                    df_results.loc[index_number_results] = [index_key, fluorescence, time]



                row = row + 1
                n = 0

            #df_results = df_results.sort_values(by=['Index Key', 'Time'], ascending=True)
            df_results = df_results.set_index(df_results['Index Key'])
        except Exception as e:
                print('EXCEPTION RAISED')
                print('FAILURE IN DATASORTING, likely an issue with file formating.\n')
                print(str(e))
                print('EXLCUDING TIMEPOINT: ' + str(time) + ' mins\n')
                continue
    print('Data sorting complete \n')
    return(df_results, df_z_prime, controls_not_inc)
