from re import X
#from turtle import pd
import io
import pandas as pd
import matplotlib as plt
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_curve(glu_lib_q, cmp_info_file, model, construct_name, glu_data_path, xdata, ydata, index, x_curve, *popt):

    output_dir = 'outputs/' + construct_name + '_outputs/'
    new_dir = output_dir + construct_name + '_graphs'

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    def func(x, a, b, c):

        return a * np.exp(-b * x) + c
    
    def linear_func(x, b, c):

        return (-b * x) + c

    if glu_lib_q == 'No':

        #If GSH data isn't present it will just plot the protein data by itself.

        df_glu_pre = pd.read_csv(glu_data_path)

        compound_info_df = pd.read_csv(cmp_info_file)
        try:
            df_glu = compound_info_df.merge(df_glu_pre, on=['Index Key', 'Compounds', 'Smile'])
        except:
            df_glu = compound_info_df.merge(df_glu_pre, on='Index Key')


        try:

            df_glu = df_glu.set_index('Index Key')

        except:

            df_glu['Index Key'] = df_glu['Unnamed: 1'] + df_glu['Unnamed: 2'].astype(str)

            df_glu = df_glu.drop(['Unnamed: 1','Unnamed: 2' ], axis=1)
            df_glu = df_glu.set_index('Index Key')


        index_slice = df_glu.loc[index]


        glu_rate = index_slice['Rate Glu']

        try:
            zcode = index_slice['Compounds']
        except:
            zcode = ' '
        

        try:
            model_glu = index_slice['Model glu']

        except:

            model_glu = 'Non-linear'

        try:
            glu_a = index_slice['a']
            glu_c = index_slice['c']

            glu_parameters = [glu_a, glu_rate, glu_c]

            plt.figure()

            plt.scatter(x = xdata, y = ydata, c = 'black')

            x_max = xdata.max()

            y_min = ydata.min()

            y_max = ydata.max()

            x_curve = np.linspace(0,x_max,100)


            if model == 'Non-linear':
                plt.plot(x_curve, func(x_curve, *popt), c = 'black', label = 'Protein')


            if model_glu == 'Non-linear':    
                plt.plot(x_curve, func(x_curve, *glu_parameters), c = 'red', label = 'Glutathione')

        except:

            glu_parameters = [1, glu_rate, 0]

            
            plt.figure()

            plt.scatter(x = xdata, y = ydata, c = 'black')

            x_max = xdata.max()

            y_min = ydata.min()

            y_max = ydata.max()

            x_curve = np.linspace(0,x_max,100)


            if model == 'Non-linear':
                plt.plot(x_curve, func(x_curve, *popt), c = 'black', label = 'Protein')


            if model_glu == 'Non-linear':    
                plt.plot(x_curve, func(x_curve, *glu_parameters), c = 'red', label = 'Glutathione')
        

    if glu_lib_q == 'Yes':

        compound_info_df = pd.read_csv(cmp_info_file)

        #Triggered if there's no GSH data
        compound_info_df = compound_info_df.set_index('Index Key')
        index_slice = compound_info_df.loc[index]

        try:
            zcode = index_slice['Compounds']
        except:
            zcode = ' '

        plt.figure()

        plt.scatter(x = xdata, y = ydata, c = 'red')


        x_max = xdata.max()

        y_min = ydata.min()

        y_max = ydata.max()

        x_curve = np.linspace(0,x_max,100)


        if model == 'Non-linear':
            plt.plot(x_curve, func(x_curve, *popt), c = 'red', label = 'Glutathione')
        
        #zcode = ' '



    plt.title(index + ' - ' + zcode)
    plt.legend()
    plt.xlabel('Time (mins)')
    plt.ylabel('Relative Fluorescence')


    y_min = ydata.min()

    y_max = ydata.max()

    if y_min > 0:
        y_min = 0
    else:
        y_min =ydata.min() - 0.1

    if y_max < 1.2:

        y_max = 1.2

    else:
        y_max = y_max + 0.2


    plt.ylim((y_min - 0.1), y_max)

    save_name = new_dir + '/' + index + 'graph.jpg'

    plt.savefig(save_name, format="jpg", bbox_inches="tight", dpi=150)



    return()

