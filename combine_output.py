# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from netCDF4 import Dataset
import glob
import re
import os
### 

# combining ensemble of SURFEX outputs into a single outputfile 
# with additional dimension for ensemble memer

parentfolder = r'/home/jpnousu/CROSIM/RESULTS/ENSEMBLE_test/SIIKANEVA/CrocO_vortex_sandbox/arch_test_ol'
folders = glob.glob(f'{parentfolder}/mb*')
folders.sort()

files = []
for f in folders:
    files.append(glob.glob(f'{f}/pro/*.nc')[0])

files.sort()
# ADD TEST TO MAKE SURE SORTED THAT MB001 is the first

# OPEN FIRST FILE TO DERIVE DIMENSIONS AND VARIABLES
data = Dataset(files[0], 'r')

variables = data.variables.keys()

#------------------------------------------------------------#

outname = re.split('(/)', files[0])[-1]
firstmember = re.split('(/)', folders[0])[-1]
lastmember = re.split('(/)', folders[-1])[-1]
outputfolder = f'{parentfolder}/{firstmember}_{lastmember}'
outfilefull = f'{outputfolder}/{outname}'

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

# creating the netcdf file: dataset
new_output = Dataset(outfilefull, 'w', format='NETCDF4_CLASSIC')

# creating the netcdf file: 
# dimensions from the example results
# + dimension for different members

for d in data.dimensions.keys():
    temp = new_output.createDimension(d, data.dimensions[d].size)

new_output.createDimension('Ensemble_members', len(files))


list_of_variables = list(data.variables.keys())
list_of_variable_units = []
list_of_variable_longnames = []

# create list of units
for v in data.variables.keys():
    try:
        list_of_variable_longnames.append(data[v].long_name)
    except AttributeError:
        list_of_variable_longnames.append('')        
    try:
        list_of_variable_units.append(data[v].units)
    except AttributeError:
        list_of_variable_units.append('')

# creating variables
for i in range(len(list_of_variables)):
    variable = list_of_variables[i]
    try:
        units = data[variable].units
    except AttributeError:
        units = False
    dtype = data[variable].dtype
    dimensions = data[variable].dimensions + ('Ensemble_members',)
    longname = list_of_variable_longnames[i]
    try:
        fill_value = data[variable]._FillValue
    except AttributeError: # this happaned to 'Projection_Type'... bad fix now
        fill_value = -2147483647
    
    output_var = new_output.createVariable(variable, dtype, dimensions,
                                        fill_value=fill_value)
    if units != False:
        output_var.units = units
    output_var.long_name = longname


##------------------

# looping over the files to read data into right dimensions

for j in range(len(files)):
    member = j
    file = files[j]
    print(member)
    tempfile = Dataset(file, 'r')
    for e in range(len(files)):
        
        for i in list_of_variables:
            if tempfile[i].dimensions != ():
                dimensionlen = len(list(tempfile[i].dimensions))
            else:
                dimensionlen = 0
            if dimensionlen == 4:
                new_output[i][:,:,:,:,e] = tempfile[i][:,:,:,:]
            elif dimensionlen == 3:
                new_output[i][:,:,:,e] = tempfile[i][:,:,:]            
            elif dimensionlen == 2:
                new_output[i][:,:,e] = tempfile[i][:,:]       
            elif dimensionlen == 1:
                new_output[i][:,e] = tempfile[i][:]       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    