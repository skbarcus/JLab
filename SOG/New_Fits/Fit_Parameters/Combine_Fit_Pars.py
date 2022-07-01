import numpy as np
import matplotlib.pyplot as plt
import glob, os

#Script to combine the SOG fit parameter output files into one larger text file.

#Set correct directory to find SOG fit parameter files.
os.chdir("/home/skbarcus/JLab/SOG/New_Fits/Fit_Parameters")
os.system("pwd")

#Find all files from the same group of fits.
Fit_Files = []
for file in glob.glob("*3H_Bootstrap24*.txt"):
    Fit_Files.append(str(file))

#print(He3_files)

#Loop over the files and print the contents to a new file. Only keep the header line info for the first file read in.

#Create output file to write SOG fit parameters to.
file_name = 'All_Fit_Pars_3H_Bootstrap24_6-23-2022.txt'
with open(file_name, 'w') as out:

    first_file=1
    for file in Fit_Files:
        file = str(file)
        print(file)
        #print(type(file))

        with open('/home/skbarcus/JLab/SOG/New_Fits/Fit_Parameters/'+file) as f:
            lines = f.readlines()
            if first_file!=1:
                del lines[0]    #Delete the duplicate column labels.

        for line in lines:
            out.write(str(line))
        first_file=0 #No longer the first file.
print('***********************************************************************')
print('Created combined file '+file_name+'.')
print('***********************************************************************')


