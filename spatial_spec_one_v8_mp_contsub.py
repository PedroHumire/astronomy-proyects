#!/usr/bin/env python
# -*- coding: utf-8 -*-   #con esto puedo poner tilde.
import numpy as np #to get numerical function available
from matplotlib.gridspec import GridSpec #to split page to select area where further plot
import matplotlib.pyplot as plt #to plot
from astropy.io import fits #to read fits' files
#from pvextractor import Path
#from pvextractor import extract_pv_slice #to make pv diagrams
from matplotlib.backends.backend_pdf import PdfPages #if a want sabe results in a PDF's file
from astroquery.splatalogue import Splatalogue 
from astropy import units as u #to fix freq. units in splatalogue
from matplotlib import pyplot as mp #to plot
from coordconverter import * #Pedro Humire R.
import glob #to read data in a directory
import matplotlib.patches as patches
import os #to create folders.
from os import walk #to print archives in a certain directory given
from astropy.utils import data
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
#from spectral_cube import SpectralCube
from scipy.interpolate import griddata
from astropy.wcs import WCS
import math
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import gc
import timeit
import time
import multiprocessing

starttime = time.time()

start_time = timeit.default_timer()

sigma_clip=0.0

if sigma_clip == 0.0:
    clipping='no'
    print('################ NO sigma clipping was used #################')
else:
    clipping='yes'
    print('################ A sigma clip of {0:.1f} was used #################'.format(sigma_clip))

current_dir=os.getcwd()   #This line should be there. Otherwise os.makedirs will complain.
results_dir=current_dir+'/EXTRACTED_SPECS_mp_contsub/output/'
################################################################ CREATE DIRECTORY

# Create target directory & all intermediate directories if don't exists
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
    print("Directory " , results_dir ,  " Created ")
else:    
    print("Directory " , results_dir ,  " already exists")   
################################################################     
    

gs = GridSpec(1, 2)				# Number of page's divisions
fig = plt.figure(figsize=(16.5, 10))


speed_of_light = 299792.458
#alchemi_data='/media/pedro/Nuevo_vol/Doctorado/alchemi/Archive/v0/'
alchemi_data=current_dir+'/'

#folder_to_read=['B4a','B4b','B4c']
#together='yes'                   #if 'yes', the program will read everything inside alchemi_data
                                                                            #if 'no', the program will read just the folder folder_to_read inside alchemi_data

folder_to_read=['B3a']
#folder_to_read=['B5h']
together='no'  
"""
galaxy_center  =  HMS2deg(ra='00 47 33.18', dec='-25 17 16.93')      #TH2, e.g. in Mangum+2019. This is where I am extracting the spectrum.
#regions according to Ando+2017:
region_1       =  HMS2deg(ra='00 47 33.30', dec='-25 17 15.6')        #region_1 or clump 1 in Ando+2017 (Tab. 1) is Region 6 in Mangum+2019
region_2       =  HMS2deg(ra='00 47 33.16', dec='-25 17 17.2')             #ra=33.2994, dec=15.5952.
region_3       =  HMS2deg(ra='00 47 33.20', dec='-25 17 16.7')             #region_1 is the zone where Mangum+2019 found most Ch3Oh masers!
region_4       =  HMS2deg(ra='00 47 33.11', dec='-25 17 17.7')
region_5       =  HMS2deg(ra='00 47 33.12', dec='-25 17 18.2')
region_6       =  HMS2deg(ra='00 47 32.99', dec='-25 17 19.7')
region_7       =  HMS2deg(ra='00 47 32.94', dec='-25 17 20.2')
region_8       =  HMS2deg(ra='00 47 32.84', dec='-25 17 21.2')

regions = [galaxy_center,region_1,region_2,region_3,region_4,region_5,region_6,region_7,region_8]


"""
##My regions

region_1		=('11.88341','-25.29118')
region_2		=('11.88449','-25.28895')
region_3		=('11.88669','-25.28932')
region_4		=('11.88739','-25.28888')
region_5		=('11.88838','-25.28817')
region_6		=('11.88888','-25.28771')
region_7		=('11.89018','-25.28702')
region_8		=('11.89176','-25.28650')
region_9		=('11.89236','-25.28674')
region_10		=('11.89265','-25.28551')

regions = [region_1,region_2,region_3,region_4,region_5,region_6,region_7,region_8,region_9,region_10]


"""
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

"""

n_reg   = len(regions)                                                      #number of regions.
#Region without emission (fixed position) to obtain a rms estimation:
rms_region     =  HMS2deg(ra='00 47 33.646', dec='-25 17 25.83')
rms_radius     =  HMS2deg(dec='-25 17 23.58')                               #15 pix are 2.25 arcseconds
#galaxy_center  =   αJ2000=00h47m33.s3,δJ2000=−25◦17′23.′′15 from Martin+2010


def multiprocessing_func(folder_names):

        single_spec_sum  = np.zeros([1000])

        spec_sum  = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]        #making n-D list

        spec_mean = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]

        sum_spec  = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]

        #for folder_names in folders:   #folder_names: alchemi's directories.

        print('################# reading folder:', folder_names)

        cubes=(glob.glob(alchemi_data+folder_names+'/*.fits'))			  #read every fits file in alchemi_data or folder_to_read
        cubes.sort()		
                                #sort the names


        
        for cubesnames in cubes:				 
                if (('.contsub.cv01_6.cube' in cubesnames) == True): #Will select only continuum-subtracted data cubes.

        ########################################################## Reading the cube: ##########################################
                        cube = fits.open(cubesnames)
                        
        ########################################################## Getting header's parameters ################################
                        RESTFRQ   = cube[0].header['RESTFRQ']
                        CRVAL1    = cube[0].header['CRVAL1']
                        CDELT1    = cube[0].header['CDELT1']
                        CRPIX1    = cube[0].header['CRPIX1']
                        CRVAL2    = cube[0].header['CRVAL2']
                        CDELT2    = cube[0].header['CDELT2']
                        CRPIX2    = cube[0].header['CRPIX2']
                        NAXIS1    = cube[0].header['NAXIS1'] #total number of pixels in X axis 
                        NAXIS2    = cube[0].header['NAXIS2'] #total number of pixels in Y axis 
                        NAXIS3    = cube[0].header['NAXIS3'] #total number of pixels in freq. axis 
                        CRVAL3    = cube[0].header['CRVAL3'] #first frequency value (Hz)
                        CRPIX3    = cube[0].header['CRPIX3'] #first pixel value (=1 for ALCHEMI) 
                        CDELT3    = cube[0].header['CDELT3'] #frequency of each channel bin. 
                        BMAJ      = cube[0].header['BMAJ']   #maximum beam 
                        w = WCS(cube[0].header)
        ##########################################################

                        x1,x2=int(NAXIS1*0.35),int(NAXIS1*0.70)
                        y1,y2=int(NAXIS2*0.35),int(NAXIS2*0.70)

                        data = cube[0].data
                        data = data[0,:,:,:] #le quitamos una dimensión.
                        z_size,y_size,x_size=data.shape[:]  	#en python, las dimensiones están al revés de 
                                                                #como uno lo piensa :z,y,x en vez de x,y,z
                        del cube

                        if clipping != 'no':
                                data_rms=np.zeros([z_size,y_size,x_size])

                        data[:,:y1,:]=0
                        data[:,:,:x1]=0
                        data[:,-(y_size-y2):,:]=0
                        data[:,:,-(x_size-x2):]=0
                        
                        #data_out=data*np.Nan                    #creating NaN array to do nansum and nanmean afterwards.
                        data_out=np.zeros([n_reg,z_size,y_size,x_size],dtype='float32') #creating NaN array to do nansum and nanmean afterwards.
                        data_out_save=data_out


                        ##################################################
                        if clipping != 'no':  #Looking for radius of the region to obtain the noise
                                
                                count=0			
                                for i in range(0,x_size):
                                        lon, lat,zd,fa = w.all_pix2world(i, 0, 0,0,0)
                                        if np.isclose(lon, float(rms_region[0]), rtol=2e-06, atol=1e-08, equal_nan=False)==True:
                                                rms_reg_x_center = count   #central longitude position (in pixel) of the region to obtain the noise 
                                        count=count+1
                                        
                                count=0			
                                for j in range(0,y_size):
                                        lon, lat,zd,fa = w.all_pix2world(0, j, 0,0,0)
                                        if np.isclose(lat, float(rms_region[1]), rtol=1e-06, atol=1e-08, equal_nan=False)==True:
                                                rms_reg_y_center = count   #central latitude position (in pixel) of the region to obtain the noise 
                                        if np.isclose(lat, float(rms_radius), rtol=1e-06, atol=1e-08, equal_nan=False)==True:
                                                rms_radius = count
                                        count=count+1
                                
                                h_rms,k_rms = rms_reg_x_center,rms_reg_y_center
                                r_rms	  = np.abs(rms_radius-k_rms)  #radius of an offset positon minus central latitude position (both in latitude) of the noise region, to get the radius (in pixels)
                                print('rms radius en pix', r_rms)
                        
                                for j in range(0,y_size):  
                                        for i in range(0,x_size):
                                                if ((i-h_rms)**2 + (j-k_rms)**2 <= r_rms**2):
                                                        for k in range (0,z_size):
                                                                data_rms[k,j,i]=data[k,j,i]                                                         
                                                else:
                                                        for k in range (0,z_size):
                                                                data_rms[k,j,i]=0.0


                                sigma_cut = sigma_clip*data_rms[data_rms != 0.0].std()
                                print('sigma en la region',sigma_cut)
                                del data_rms
                                
                                
                                for j in range(0,y_size):  
                                        for i in range(0,x_size):
                                                for k in range (0,z_size):
                                                        if (data[k,j,i] < sigma_cut):  # sigma clip  #ERROR TODOS LOS DATOS BAJO 0 SE TRANSFORMAN EN 0!!!!
                                                                data[k,j,i]=0.0   

                        ##################################################

                        r = float(abs(BMAJ/CDELT2))/2. #resolution
                        h_c,k_c = [],[]
                        xaxis=(np.arange(NAXIS3)+1 - CRPIX3)*CDELT3+CRVAL3    #ALTERNATIVE WAY: xaxis=np.linspace(CRVAL3,CDELT3*(z_size-CRPIX3)+CRVAL3,z_size)
                
                        for region in regions: 
                        
                                count=0			
                                for i in range(0,x_size):
                                        lon, lat,zd,fa = w.all_pix2world(i, 0, 0,0,0)	
                                        if np.isclose(lon, float(region[0]), rtol=2e-06, atol=1e-08, equal_nan=False)==True:
                                                region_x_center_pix = count
                                        count=count+1

                                count=0			
                                for j in range(0,y_size):
                                        lon, lat,zd,fa = w.all_pix2world(0, j, 0,0,0)
                                        if np.isclose(lat, float(region[1]), rtol=1e-06, atol=1e-08, equal_nan=False)==True:
                                                region_y_center_pix = count
                                        count=count+1
                                
                                h_c.append(region_x_center_pix),k_c.append(region_y_center_pix) 

                        h_c,k_c = np.array(h_c, dtype=np.float32),np.array(k_c, dtype=np.float32)  #center of the region in pixels.

                        ##################################################  Cut datacube in interested CIRCULAR region 

                        for region in range(0,n_reg):                        
                                area_in_pix=0
                                for i in range(0,x_size):
                                        for j in range(0,y_size):
                                                if ((i-h_c[region])**2 + (j-k_c[region])**2 <= r**2):
                                                        data_out[region,:,j,i]=data[:,j,i]
                                                        area_in_pix=area_in_pix+1
                                                                

                                for ka in range(0,z_size):
                                        single_spec_sum[ka]  = np.sum(data_out[region,ka,:,:])   #Creation of the spectra. 
                                                                

                                sum_spec[region].append(single_spec_sum[:int(z_size)])

                                spec_sum[region]=np.reshape(sum_spec[region],(int(z_size),-1))  #from [1,2,3,4] to [[1],
                                                                                                                #                   [2],     
                                                                                                                #                   [3],
                                                                                                                #                   [4]])
                                spec_sum[region]=np.array(spec_sum[region], dtype=np.float32)             
                                spec_mean[region]=spec_sum[region]/area_in_pix
                                #print((spec_sum[0]/spec_mean[0])[0],area_in_pix)
                                data_out=data_out_save
                        
                        
                        ################################################## Saving data.
                        xaxis=np.array(xaxis, dtype=np.float32)
                        
                        data_sum  = Table([xaxis,spec_sum[0],spec_sum[1],spec_sum[2],spec_sum[3],spec_sum[4],spec_sum[5],spec_sum[6],spec_sum[7],spec_sum[8],spec_sum[9]])

                        data_mean  = Table([xaxis,spec_mean[0],spec_mean[1],spec_mean[2],spec_mean[3],spec_mean[4],spec_mean[5],spec_mean[6],spec_mean[7],spec_mean[8],spec_mean[9]])

                        np.savetxt(results_dir+'ALCHEMI_sum_'+cubesnames[len(alchemi_data)+len(folder_names)+1:-5]+'_.dat', data_sum)
                        np.savetxt(results_dir+'ALCHEMI_mean_'+cubesnames[len(alchemi_data)+len(folder_names)+1:-5]+'_.dat', data_mean)

                        del data
                        del data_out
                        del data_out_save
                        del data_sum
                        del data_mean
                        single_spec_sum=np.zeros([1000])
                                                                                                                                        
                        spec_sum  = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]        #making n-D list

                        spec_mean = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]

                        sum_spec  = [[], [], [], [], [], [], [], [], [], [], [], [], [], []]



################################################################ RUNNING THE CODE #####################################################


######################## Reading the data
folders = []
if (together == 'yes'):                                 
	for (dirpath, dirnames, filenames) in walk(alchemi_data): #"or you could use os.walk() which will yield two lists for each directory
		    folders.extend(dirnames)			  #it visits - splitting into files and dirs for you. If you only want the
		    break					  #top directory you can just break the first time it yields."
		  						  #print 'folders that will be read:', folders
else:
	folders=folder_to_read


print('folders that will be read:', folders)


if __name__ == '__main__':
    starttime = time.time()
    processes = []
    for i in folders:
        p = multiprocessing.Process(target=multiprocessing_func, args=(i,))
        processes.append(p)
        p.start()
        
    for process in processes:
        process.join()

                        
#print(timeit.default_timer() - start_time)

print('That took {} seconds'.format(time.time() - starttime))
print('       or {} minutes'.format((time.time() - starttime)/60.))
 
 
