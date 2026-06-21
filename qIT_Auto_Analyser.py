#!/usr/bin/python3

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
from qit_plotter import plot_curve
from data_sorting import sort_data
from modelling import model
from report_builder import write_report
import shutil

import io



#This code is for analysing the data obtained from the Quantitative irreversible tethering (qIT) assay
#Any suggestions for how to improve the code should be emailed to c.brown22@imperial.ac.uk :)


#First imports the files from the clariostar plate reader (needs to be in the format: 'x_mins.csv') and puts them into a dataframe.
#Note this will only work if .csv files are converted to .csv in mircosoft excel.

#The file name with the data in must be in the format 'files_x' where x = construct name.

#COULD ADD ENTER GSH FILE FUNCTIONALITY



def qit_analysis():

    print('************* PLEASE READ *************\n')
    print('Hello, welcome to the qIT Auto-analyser.\n -This software is used to analyse the results of qIT screens in a 384 well format using the standard output files from a Clariostar plate reader using mars software. \n -See the example data for an example of the file formating, any changes to this standard format (e.g. an extra line or an incomplete plate) will fail.\n -Input files should be in the form of .csv files and should always be named following the pattern - x mins.csv. ANY CHANGES TO THIS FORMAT WILL NOT WORK.\n')
    print('-If you would like your compound structures and descriptor information to appear in the report, an appropriately labelled datawarror file (e.g. DHC_library.dwar) must be present in the current working directory.')
    print('To prepare data to analyse, please do the following:')
    print('STEP 1. Convert all input files to csv files by opening in excel > save as > .csv. Name the file following the format - x mins.csv.')
    print('STEP 2. Combine these files into a new folder named files_XXX (XXX = a relevant name, e.g. files_IDH1TAP1)')
    print('STEP 3. Place this folder in the inputs folder (if you do not alreay have one, a folder titles inputs should be created in the same folder as this module)')
    print('STEP 4. Ensure the correct GSH data is in the folder marked Glu_data. \n')
    print('-If these instructions have been followed but there are still issues, please contact charliepeterbrown@gmail.com\n')
    construct_name = input("Which dataset would you like to analyse? (Please input the name of you used in the input folder, e.g. if my file was named files_IDH1 I would input IDH1)")

    file_location = 'inputs/' + 'files_' + construct_name
    raw_files = os.listdir(file_location)
    print(raw_files)
    df_results = pd.DataFrame(columns = ['Index Key', 'Fluorescence', 'Time'])

    glu_files = os.listdir('Glu_data')
    glu_files = [w[:-10] for w in glu_files]

    old_output_files_directory = 'outputs/' + construct_name + '_outputs'
    
    try:
        old_output_files = os.listdir(old_output_files_directory)

        if os.path.isdir(old_output_files_directory):
            shutil.rmtree(old_output_files_directory)
    except:
        pass
    
    question_glu_choice = [
    inquirer.List('glu_lib_q',
                    message="Is this GSH data you are trying to analyse? (If so, it will not be compared to past GSH data)",
                    choices= ['Yes', 'No'],
                ),
    ]
    
    glu_lib_q = inquirer.prompt(question_glu_choice)
    glu_lib_q = glu_lib_q['glu_lib_q']

    if glu_lib_q == 'Yes':
        glu_lib = construct_name

    if glu_lib_q == 'No':

        question_glu_lib = [
        inquirer.List('glu_lib',
                        message="Which library did you use?",
                        choices= glu_files,
                    ),
        ]

        glu_lib = inquirer.prompt(question_glu_lib)

        glu_lib = glu_lib['glu_lib']

        glu_data_path = 'Glu_data/' + glu_lib + '_rates.csv'
        print('GLU DATA PATH')
        print(glu_data_path)

    else:

        #If GSH data is beign analysed it will save this data into the Glu data folder.

        glu_data_path = 'Glu_data/' + construct_name + '_rates.csv'

#'glu_lib_2' holds the file containing Compound ID and Smiles strings.

    glu_files_2 = os.listdir('Glu_data/Compound_Information')
    print(glu_files_2)
    

    question_glu_lib_2 = [
    inquirer.List('glu_lib_2',
                    message="Which library did you use?",
                    choices= glu_files_2,
                ),
    ]

    glu_lib_2 = inquirer.prompt(question_glu_lib_2)

    glu_lib_2 = glu_lib_2['glu_lib_2']


    glu_data_path = 'Glu_data/' + glu_lib + '_rates.csv'
        

    #Puts all the data into a df

    output_dir = 'outputs/' + construct_name + '_outputs'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_dir = output_dir + '/'

    output = sort_data(df_results, raw_files, file_location)

    df_results = output[0]

    df_z_prime = output[1]

    controls_not_included = output[2]

    df_results.to_csv(output_dir + construct_name + '_RESULTS.csv')

#Next it goes through each well and calculates the rate of reaction (using non-linear regression and assuming pseudo-first order kinetics) and the R squared (an imperfect stastic here, but useful to understand how close the curve fits to the data).

    index_list= df_results['Index Key'].tolist()

    index_list = list( dict.fromkeys(index_list))

#Tries to model data using a non-linear regression model, returns a dataframe with results in.
    cmp_info_file = 'Glu_data/Compound_Information/' + glu_lib_2

    df_r2_rate = model(glu_lib_q, cmp_info_file, construct_name, glu_data_path, index_list, df_results)

    df_r2_rate.to_csv(output_dir + construct_name + 'RESULTS_r2_rate.csv')

#Opens gluathione data, calculates REF and saves results to CSV files.

    if glu_lib_q == 'No':
        df_glu_pre = pd.read_csv(glu_data_path)

            
        compound_info_df = pd.read_csv(cmp_info_file)

        try:

            df_glu = compound_info_df.merge(df_glu_pre, on=['Index Key', 'Compounds', 'Smile'])

        except:

            df_glu = compound_info_df.merge(df_glu_pre, on='Index Key')



        try:

            df_merged = df_glu.merge(df_r2_rate, on=['Index Key', 'Compounds', 'Smile'])

        except:
            df_merged = df_glu.merge(df_r2_rate, on='Index Key')


        if glu_lib_q == 'No':

            df_merged["Rate"] = pd.to_numeric(df_merged["Rate"], errors='coerce')

            df_merged["Rate Glu"] = pd.to_numeric(df_merged["Rate Glu"], errors='coerce')

            df_merged['REF'] = df_merged['Rate'] / df_merged['Rate Glu']




    if glu_lib_q == 'Yes':
    #Merges with compound info dataframe, renames 'Rate' to 'Rate Glu' and saves to appropriate folder.
        df_merged_1 = df_r2_rate
        compound_info_df = pd.read_csv(cmp_info_file)

        try:
            df_merged = compound_info_df.merge(df_merged_1, on=['Index Key', 'Compounds', 'Smile'])
        except:
            df_merged = compound_info_df.merge(df_merged_1, on='Index Key')

        df_gsh_output = df_merged[['Index Key', 'Rate', 'a', 'c']].copy()
        df_gsh_output.rename(columns={'Rate': 'Rate Glu'}, inplace=True)
        df_gsh_output.to_csv(glu_data_path)

    df_good_data = df_merged[(df_merged['R Sqrd']) > 0.7]


    df_merged.to_csv(output_dir + 'Final_results' + construct_name + '.csv')
    df_good_data.to_csv(output_dir + 'Good_data' + construct_name + '.csv')

    #Next saves GSH data to correct folder.

    write_report(df_z_prime, construct_name, df_merged, controls_not_included, glu_lib)


        

    return

qit_analysis()
