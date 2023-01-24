#!/usr/bin/env python
# -*- coding: utf-8 -*-   #con esto puedo poner tilde.

"""
Code to stick all spectra extracted from ALCHEMI's cubes by stracted_spec code.
"""

from matplotlib.gridspec import GridSpec #to split page to select area where further plot
import numpy as np #to get numerical function available
import matplotlib.pyplot as plt #to plot
import glob #to read data in a directory
from os import walk #to print archives in a certain directory given
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from coordconverter import * #Pedro Humire R. #this is required for the HMS2deg function.
import os #to create folders.
from matplotlib.backends.backend_pdf import PdfPages #multi pdf pages.
import scipy.constants as c

redshift=0.00081 #preferred z from NED IPAC
x_axis_save    		= []
x_axis_total   		= []
x_axis_total2	 	= []
x_axis_total3	 	= []
y_axis_save    		= [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
y_axis_total   		= [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
y_axis_mean_save    = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
y_axis_mean_total   = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
T_mb        		= [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
T_mb_mean   		= [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
overlaps_start 		= []
overlaps_end   		= []
overlaps_mean_start = []
overlaps_mean_end   = []
########################################################################
current_dir=os.getcwd()   #This line should be there. Otherwise os.makedirs will complain.
dir_input=current_dir+'/EXTRACTED_SPECS_mp_contsub/output/'
dir_output =current_dir+'/EXTRACTED_SPECS_mp_contsub/data/'
alchemi_data='/media/pedro/Nuevo_vol/Doctorado/alchemi/Archive/v0/'
alchemi_data=current_dir+'/'
########################################################################

###########################################################
# Create target directory & all intermediate directories if don't exists
for new_dir in [dir_output]:
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
        print("Directory " , new_dir ,  " Created ")
    else:    
        print("Directory " , new_dir ,  " already exists")   
###########################################################  

region_1		=('11.88341','-25.29118')
region_2		=('11.88449','-25.28895')
region_3		=('11.88669','-25.28932')
region_4		=('11.88739','-25.28888')
region_5		=('11.88838','-25.28817') #Galactic Center.
region_6		=('11.88888','-25.28771')
region_7		=('11.89018','-25.28702')
region_8		=('11.89176','-25.28650')
region_9		=('11.89236','-25.28674')
region_10		=('11.89265','-25.28551')
regions = [region_1,region_2,region_3,region_4,region_5,region_6,region_7,region_8,region_9,region_10]

### Regions for Yaoting
region_1		=('11.89229583','-25.28696000')
region_2		=('11.89156667','-25.28635750')
region_3		=('11.89020417','-25.28701750')
region_4		=('11.88890000','-25.28764861')
region_5		=('11.88807917','-25.28827972')
region_6		=('11.88883750','-25.28845167')
region_7		=('11.88737917','-25.28879583')
region_8		=('11.88667917','-25.28922611')
region_9		=('11.88468333','-25.28876722')
region_10		=('11.88442917','-25.28982861')
region_11		=('11.88385833','-25.29074667')
region_12		=('11.88306250','-25.29146361')
regions = [region_1,region_2,region_3,region_4,region_5,region_6,region_7,region_8,region_9,region_10,region_11,region_12]

n_reg   = len(regions)

###########################################################  

folder_to_read=['B3a']                                                     #if I want to select just some of the bands
together='no'                                                              #if 'yes', the program will read everything inside alchemi_data
                                                         #if 'no', the program will read just the folder folder_to_read inside alchemi_data
######################## Reading the data
folders = []     
if (together == 'yes'):                                 
	for (dirpath, dirnames, filenames) in walk(alchemi_data): #"or you could use os.walk() which will yield two lists for each directory
			folders.extend(dirnames)			  #it visits - splitting into files and dirs for you. If you only want the
			break					  #top directory you can just break the first time it yields."
									  #print 'folders that will be read:', folders
else:
	folders=folder_to_read


print ('folders that will be read:', folders)
n_fol=0
for folder_names in folders:   #folder_names: carpetas de alchemi: asd, B6e,B4f, B4c una por una

	print ('################# reading folder:', folder_names,n_fol)
	n_fol=n_fol+1
	cubes=(glob.glob(alchemi_data+folder_names+'/*.fits'))			  #read every fits file in dir_cube
	cubes.sort()		
				  #sort the names



	for cubesnames in cubes:
		if (('.contsub.cv01_6.cube' in cubesnames) == True): 	          #this is to select just data with 1.6 arcsec resolution
		#if (('.085000.contsub.cv01_6.cube' in cubesnames) == True): 	          #this is to select just data with 1.6 arcsec resolution


			data = ascii.read(dir_input+'ALCHEMI_sum_'+cubesnames[len(alchemi_data)+len(folder_names)+1:-5]+'_.dat', delimiter=' ')  #len(alchemi_data)+len(folder_names)+1 is to avoid the root of the file and the last slash "/" character
			data_mean = ascii.read(dir_input+'ALCHEMI_mean_'+cubesnames[len(alchemi_data)+len(folder_names)+1:-5]+'_.dat', delimiter=' ')  #len(alchemi_data)+len(folder_names)+1 is to avoid the root of the file and the last slash "/" character

			x_axis_save.append(data[0][:]) #the same for data and data_mean
			overlaps_start.append(data[0][0])            #start and end of x axis of each cube.
			overlaps_end.append(data[len(data)-1][0])
			overlaps_mean_start.append(data_mean[0][0])            #start and end of x axis of each cube.
			overlaps_mean_end.append(data_mean[len(data_mean)-1][0]) 


			for dim in range(0,len(data[0])-1):#range until -1 because xaxis do one dim less
				y_axis_save[dim].append(data[dim+1][:]) #dim+1 because data[0] is the xaxis
				y_axis_mean_save[dim].append(data_mean[dim+1][:])



for a in x_axis_save:   #https://www.geeksforgeeks.org/nested-list-comprehensions-in-python/ 
	for b in a:     #left the whole list in 1D
		x_axis_total.append(b)
		x_axis_total2.append(b)

for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 was the x axis originally
	for a in y_axis_save[dim]:
		for b in a:
			y_axis_total[dim].append(b)
            
for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 was the x axis originally
	for a in y_axis_mean_save[dim]:
		for b in a:
			y_axis_mean_total[dim].append(b)

#############################################################################################
for dim in range(0,len(data[0])-1):     #sort y axis with respect to a further sorted x axis 
	y_axis_total[dim]      = [x for _,x in sorted(zip(x_axis_total,y_axis_total[dim]))]
	y_axis_mean_total[dim] = [x for _,x in sorted(zip(x_axis_total,y_axis_mean_total[dim]))]
x_axis_total=np.sort(x_axis_total)	#then sort x axis definitely

#############################################################################################

for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 is the x axis
	for i in range(0,len(x_axis_total)):
		T_mb[dim].append(1.222*1E3*y_axis_total[dim][i]*(1.0/(1.6**2*(x_axis_total[i]*1e-9)**2))*1E3)

for dim in range(0,len(data_mean[0])-1):    
	for i in range(0,len(x_axis_total)):
		T_mb_mean[dim].append(1.222*1E3*y_axis_mean_total[dim][i]*(1.0/(1.6**2*(x_axis_total[i]*1e-9)**2))*1E3)

#################################################################################################################################
########################################## Here is where I delete repeated fluxes!! #############################################
#################################################################################################################################


for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 was the x axis originally

	for j in range(1,len(x_axis_total)):

		if (np.abs(x_axis_total[j]-x_axis_total[j-1]) == 0):

			y_axis_total[dim][j]=(y_axis_total[dim][j]+y_axis_total[dim][j-1])/2.
			T_mb[dim][j]=(T_mb[dim][j]+T_mb[dim][j-1])/2.

			y_axis_mean_total[dim][j]=(y_axis_mean_total[dim][j]+y_axis_mean_total[dim][j-1])/2.
			T_mb_mean[dim][j]=(T_mb_mean[dim][j]+T_mb_mean[dim][j-1])/2.

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


save_first_xy_value_duplicated=[] #x'D
one_time=0.0

new_x_axis=np.zeros(len(x_axis_total)+1)

new_T_mb_mean=np.zeros([15,len(T_mb_mean[0])])      
new_T_mb=np.zeros((15,len(T_mb[0])))      
new_y_axis_mean_total=np.zeros((15,len(y_axis_mean_total[0])))      
new_y_axis_total=np.zeros([15,len(y_axis_total[0])])      

new_x_axis_list				= []
new_T_mb_mean_list			= [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
new_T_mb_list				= [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
new_y_axis_mean_total_list	= [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
new_y_axis_total_list		= [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 is the x axis
	count=0	
	for j in range(1,len(x_axis_total)):

		if (np.abs(x_axis_total[j]-x_axis_total[j-1]) == 0):
			#print('no',j)							#uncomment to show which pixels aren't taken into account (duplicated)
			if one_time==0.0:
				save_first_xy_value_duplicated.append(j)
			one_time=one_time+1
		else:
			if dim==0:
				new_x_axis[count]=x_axis_total[j]

			new_y_axis_total[dim][count]=y_axis_total[dim][j]
			new_T_mb[dim][count]=T_mb[dim][j]
			new_y_axis_mean_total[dim][count]=y_axis_mean_total[dim][j]
			new_T_mb_mean[dim][count]=T_mb_mean[dim][j]

			#print('sí',j)							#uncomment to show which pixels are taken into account (no duplicated)
		count=count+1
	#print('dim original: ',len(x_axis_total),len(T_mb_mean[0]))
	#print('first x (and y) value duplicated:',save_first_xy_value_duplicated)



for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 is the x axis
	for i in range(0,len(new_x_axis)):
		if new_x_axis[i] !=0.0:
			if dim==0:
				new_x_axis_list.append(new_x_axis[i])

			new_y_axis_total_list[dim].append(new_y_axis_total[dim][i])
			new_T_mb_list[dim].append(new_T_mb[dim][i])
			new_y_axis_mean_total_list[dim].append(new_y_axis_mean_total[dim][i])
			new_T_mb_mean_list[dim].append(new_T_mb_mean[dim][i])

#print('dim final: ',len(new_x_axis_list),len(new_y_axis_total_list[0]), len(new_T_mb_mean_list[0]))

ndx_l=len(new_x_axis_list)  #non-duplicated x axis length
nd_x_axis_total=np.zeros([ndx_l]) #non-duplicated x axis

nd_y_axis_total=np.zeros((15,ndx_l)) #non-duplicated y axis
nd_T_mb=np.zeros((15,ndx_l)) #non-duplicated y axis
nd_y_axis_mean_total=np.zeros((15,ndx_l)) #non-duplicated y axis
nd_T_mb_mean=np.zeros((15,ndx_l)) #non-duplicated y axis



for dim in range(0,len(data[0])-1):     							#for Jy/beam (sum)
	count=0
	for a,b in zip(new_x_axis_list,new_y_axis_total_list[dim]):
		if dim==0:
			nd_x_axis_total[count]=float(a)
		nd_y_axis_total[dim][count]=float(b)
		count=count+1

for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 is the x axis
	count=0
	for a,b in zip(new_x_axis_list,new_T_mb_list[dim]):
		nd_T_mb[dim][count]=float(b)
		count=count+1

for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 is the x axis
	count=0
	for a,b in zip(new_x_axis_list,new_y_axis_mean_total_list[dim]):
		nd_y_axis_mean_total[dim][count]=float(b)
		count=count+1

for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 is the x axis
	count=0
	for a,b in zip(new_x_axis_list,new_T_mb_mean_list[dim]):
		nd_T_mb_mean[dim][count]=float(b)
		count=count+1

#print('dim final (check): ',len(nd_x_axis_total),len(nd_x_axis_total))



for dim in range(0,len(data[0])-1):     #range finish at -1 because 0 was the x axis originally

	store_data     = Table([ [float(i) for i in nd_x_axis_total], [float(i) for i in nd_y_axis_total[dim]]])              #Save the data SUM
	store_data_Tmb = Table([ [float(i) for i in nd_x_axis_total], [float(i) for i in nd_T_mb[dim]]])				#dim+1 because region's name
	ascii.write(store_data, dir_output+'ALCHEMI_stick_Jy_beam_noclip_sum_region_'+str(dim+1)+'_.dat', overwrite=True) #started at region_1
	ascii.write(store_data_Tmb, dir_output+'ALCHEMI_stick_Tmb_noclip_sum_region_'+str(dim+1)+'_.dat', overwrite=True) #not region_0

	store_data_aver     = Table([ [float(i) for i in nd_x_axis_total], [float(i) for i in nd_y_axis_mean_total[dim]]])
	store_data_Tmb_aver = Table([ [float(i) for i in nd_x_axis_total], [float(i) for i in nd_T_mb_mean[dim]]])
	ascii.write(store_data_aver, dir_output+'ALCHEMI_stick_Jy_beam_noclip_mean_region_'+str(dim+1)+'_.dat', overwrite=True)
	ascii.write(store_data_Tmb_aver, dir_output+'ALCHEMI_stick_Tmb_noclip_mean_region_'+str(dim+1)+'_.dat', overwrite=True)
	#print(dim)
	#print('x limits: ',np.min(nd_x_axis_total),np.max(nd_x_axis_total))

