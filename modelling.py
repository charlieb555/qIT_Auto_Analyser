from cmath import nan
from re import X
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

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
from qit_plotter import plot_curve
from data_sorting import sort_data

def model(glu_lib_q, cmp_info_file, construct_name, glu_data_path, index_list, df_results):
    
    df_r2_rate = pd.DataFrame(columns=['Index Key', 'Rate', 'Half Life', 'R Sqrd', 'Initial Fluor < 0.4', 'Model', 'a', 'c'])
    
    print('Performing modelling...')

    for index in index_list:

        try:

            if (index[1:] == '1') or (index[1:] == '2') or (index[1:] == '23') or (index[1:] == '24'):
                continue

            df_one_index = df_results[(df_results['Index Key'] == index)]

            def func(x, a, b, c):

                return a * np.exp(-b * x) + c

            def linear_func(x, b, c):

                return (-b * x) + c

            df_one_index['Time'] = df_one_index['Time'].astype('int')

            df_one_index = df_one_index.sort_values(by=['Time'])

            fluor = df_one_index['Fluorescence']
            initial_fluor = (df_one_index['Fluorescence'].to_list())[0]

            xdata = df_one_index['Time'].to_numpy()

            ydata = df_one_index['Fluorescence'].to_numpy()

            if initial_fluor < 0.4:

            #Helps the modelling if it is rapidly quenched

                xdata = np.append(xdata,[0])
                ydata = np.append(ydata,[1])

                crasher = 'Yes'

            else:
                crasher = 'No'


            try:


                #First it tries to fit a curve with no pre-defined bounds.

                popt, pcov = curve_fit(func, xdata, ydata, maxfev=50000, p0=[0.9, 0.002, 0.2])


                n = len(xdata)

                y_pred = np.empty(n)

                for i in range(n):
                    y_pred[i] = func(xdata[i], popt[0], popt[1], popt[2])

                R2_score = r2_score(ydata, y_pred)

                
                if (R2_score < 0.5):

                    raise Exception
                if (popt[1] < 0.00001):
                    raise Exception
                if (popt[1] > 0.5):
            
                    raise Exception



            except:

                try:
                
                #If this doesn't work, bounds are inserted.

                    popt, pcov = curve_fit(func, xdata, ydata, maxfev=50000, bounds=((0.4, 0, -np.inf), (1.2, 0.5, 0.5)), p0=[0.9,0.002,0.2])

                    n = len(xdata)

                    y_pred = np.empty(n)

                    for i in range(n):
                        y_pred[i] = func(xdata[i], popt[0], popt[1], popt[2])
                    R2_score = r2_score(ydata, y_pred)
                    

                except:
                    #If both these fail, the experimental data is likely problematic.

                    print('Not able to find curve for ' + index)
                            
                    if (len(df_r2_rate) == 0):
                        index_number_df_r2_rate = 0
                    else:
                        index_number_df_r2_rate = len(df_r2_rate)
                    df_r2_rate.loc[index_number_df_r2_rate] = [index, None, None, None, crasher, 'Non-linear', None, None]
                    
                    continue
            
            #Is this needed???
            if (-1 > R2_score):
                
                try:
                    popt, pcov = curve_fit(func, xdata, ydata, maxfev=5000, bounds=((-np.inf, 0, -np.inf), (np.inf, 2, np.inf)))
                    n = len(xdata)

                    y_pred = np.empty(n)

                    for i in range(n):
                        y_pred[i] = func(xdata[i], popt[0], popt[1], popt[2])
                        R2_score = r2_score(ydata, y_pred)

                except:
                    print('Not able to find curve for ' + index)


                    if (len(df_r2_rate) == 0):
                        index_number_df_r2_rate = 0
                    else:
                        index_number_df_r2_rate = len(df_r2_rate)
                    df_r2_rate.loc[index_number_df_r2_rate] = [index, None, None, None, crasher, 'Non-linear', None, None]                   
                    continue

            

            x_curve = np.linspace(0,3000,100)
            
            
            rate = popt[1]
                    
            

            construct_name = str(construct_name)

            try:
                model = 'Non-linear'
                plot_curve(glu_lib_q, cmp_info_file, model, construct_name, glu_data_path, xdata, ydata, index, x_curve, *popt)
            except Exception as e:
                print('*** Failed to plot graph ***')
                print(e)

            if (len(df_r2_rate) == 0):
                index_number_df_r2_rate = 0
            else:
                index_number_df_r2_rate = len(df_r2_rate)

            half_life = -(np.log((0.5-popt[2])/(popt[0]))/(popt[1]))

                


            df_r2_rate.loc[index_number_df_r2_rate] = [index, popt[1], half_life, R2_score, crasher, model, popt[0], popt[2]]

        except Exception as e:
            print('ERROR IN MODELING')
            print(str(e))
            continue
    print('Modelling complete')
    
    return(df_r2_rate)
