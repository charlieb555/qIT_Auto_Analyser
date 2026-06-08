import PyPDF2
import matplotlib.pyplot as plt
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import rdkit
from PIL import Image, ImageDraw, ImageFilter
from fpdf import FPDF

#STILL NEED OPTION FOR NO GSH DATA...

def write_report(df_z_prime, name, df_all_descriptors, controls_not_included, glu_lib):
    print('Writing repprt...')
    output_dir = 'outputs/' + name + '_outputs/'

    output = output_dir + name + '_report.pdf'

    graphs = os.listdir(output_dir + name + '_graphs')
        
    write_first_page(name, df_z_prime, controls_not_included, output_dir)
    PDFmerge(graphs, output, name, df_all_descriptors, output_dir)
    print('Report complete')

def PDFmerge(graphs, output, name, df_all_descriptors, output_dir):
    #Creating pdf file merger object

    pdfMerger = PyPDF2.PdfMerger()

    z_prime_save_name = output_dir + name + '_graphs/' + name +  '_report.pdf'

    pdfMerger.append(z_prime_save_name)


    graphs.sort()
 
    # appending pdfs one by one

    for graph in graphs:

        index = graph[:-9]

        try:
            
            fused_name = index + 'fused.pdf'
            indexed_data = df_all_descriptors.loc[(df_all_descriptors["Index Key"] ==  index)]
            smile = indexed_data['Smile'].values[0]


            #This fuses the molecular structure image and graph together and saves as one pdf file.
            savelocation = output_dir + name + '_graphs/' + fused_name

            m = Chem.MolFromSmiles(smile)
            img = Draw.MolToImage(m)
            im1 = Image.open(output_dir + name + '_graphs/'+ index + 'graph.jpg')
            
            ws = Image.open('whitespace.jpg')

            ws.paste(img,(380,0))

            ws.paste(im1,(75,300))

            im_1 = ws.convert('RGB')
            
            im_1.save(savelocation, 'PDF', save_all=True)
        except:
            #If the image of the compound can't be generated, this part ensures it is left blank.

            savelocation = output_dir + name + '_graphs/' + fused_name

            ws = Image.open('whitespace.jpg')

            im1 = Image.open(output_dir + name + '_graphs/'+ index + 'graph.jpg')

            ws.paste(im1,(75,300))

            im_1 = ws.convert('RGB')
            
            im_1.save(savelocation, 'PDF', save_all=True)
            


        try:
            location = output_dir + name + '_graphs/' + fused_name
            pdfMerger.append(location)

        except:
            print('Error for ' + graph)

        # writing combined pdf to output pdf file
    with open(output, 'wb') as f:
        pdfMerger.write(f)
    
    files = os.listdir(output_dir + name + '_graphs')

    for file in files:
        if (file[2:] == 'fused.pdf'):

            os.remove(output_dir + name + '_graphs/' + file)

        elif (file[3:] == 'fused.pdf'):

            os.remove(output_dir + name + '_graphs/' + file)
        else:
            continue

def plot_logp(name, df_all_descriptors, output_dir):

    plt.figure()
    plt.scatter(y = df_all_descriptors['REF'], x = df_all_descriptors['MolLogP'])
    plt.ylim(0, 100)
    plt.xlabel('LogP')
    plt.ylabel('REF')
    p_save_name = output_dir + name + '_graphs/' + name + '_logP.pdf'
    plt.savefig(p_save_name, format="pdf", bbox_inches="tight")

def write_first_page(name, df_z_prime, controls_not_included, output_dir):
    filename = output_dir + name + '_graphs/' + name + '_report.txt'
    df_z_prime['Time'] = df_z_prime['Time'].astype(int)
    df_z_prime = df_z_prime.sort_values(by=['Time'])

    mean_z_prime = df_z_prime['Z Prime'].mean()

    if (mean_z_prime > 0.7):
        grade = 'First'
    if (0.5 < mean_z_prime < 0.7):
        grade = '2:1'
    if (0.3 < mean_z_prime < 0.5):
        grade = '2:2'
    if (0.3 > mean_z_prime ):
        grade = 'Fail :('


    key_list = list(controls_not_included.keys())
    key_list = [int(i) for i in key_list]
    key_list.sort()

    with open(filename, 'w') as f:
        f.write(name)
        f.write('\n')
        f.write('\n')
        f.write('Z Prime Values')
        f.write('\n')
        f.write('\n')

        dfAsString = df_z_prime.to_string(header=True, index=False)
        f.write(dfAsString)

        f.write('\n')
        f.write('\n')
        #f.write('Grade: ' + grade)
        f.write('\n')
        f.write('\n')

        f.write('Control wells not included:')

        f.write('\n')
        f.write('\n')

        f.write('Time (mins) - No. of Failures')

        f.write('\n')
        f.write('\n')

        for key in key_list:
            f.write(str(key) + ' - ' + str(controls_not_included[str(key)]))
            f.write('\n')



    pdf = FPDF()      
    # Add a page 
    pdf.add_page()  
    # set style and size of font  
    # that you want in the pdf 
    pdf.set_font("Arial", size = 15)
    # open the text file in read mode 
    f = open(filename) 
    # insert the texts in pdf 
    for x in f: 
        pdf.cell(50,5, txt = x, ln = 1, align = 'C') 
    # save the pdf with name .pdf 
    pdf.output(output_dir + name + '_graphs/' + name + '_report.pdf')
