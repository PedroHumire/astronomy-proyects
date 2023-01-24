#!/usr/bin/python

import numpy as np
from astropy.io import fits
import astropy.constants as c
import astropy.units as u
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec #to split page to select area where further plot
from scipy import ndimage as ndi


def remove_hd_dim(hd, dim=3):

    """
    Keep information from datacube header
    only related to dimension to keep
    
    hd: astropy.header instance
    
    dim: int
    dimensions to keep, default=3
    
    """

    if dim == 2:
        nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS2'],hd['NAXIS1']])).header    
        dimlist = ['1','2']
    if dim == 3:
        nhd = fits.PrimaryHDU(np.zeros([hd['NAXIS3'],hd['NAXIS2'],hd['NAXIS1']])).header
        dimlist = ['1','2','3']        
        
    for i in dimlist:
        for t in ['CRVAL','CRPIX','CDELT','CTYPE','CROTA','CUNIT']:
            if hd.get(t+i) != None:
                nhd[t+i] = hd[t+i]

    for t in ['BUNIT','BMAJ','BMIN','BPA','RESTFRQ']:
        if hd.get(t) != None:
            nhd[t] = hd[t]

    return nhd #new header

def make_mommaps(hdu, zpix ,momlist = [0]):

    """
    Function to generate moment maps
    
    hdu: astropy.hdu instance
    
    momlist: list of int
    moment maps to generate, default=[0,1,2]
    
    """
    
    data = hdu.data[zpix[0]:zpix[1],:,:]  #here I can select the region where my line of interest is, in pixels, Z axis.
    hd = hdu.header

    naxis3 = hd['NAXIS3']
    crpix3 = hd['CRPIX3']
    cdelt3 = hd['CDELT3']
    crval3 = hd['CRVAL3']

    vaxis = (np.arange(naxis3)+1 - crpix3)*cdelt3+crval3
    vaxis = vaxis * 1e-3
    
    velcube = np.zeros(data.shape)
    for v in range(data.shape[0]):
        velcube[v,:,:] = vaxis[v]

    mhd = hd.copy()
    mhd = remove_hd_dim(hd, dim=2)

    moms = []

    if 0 in momlist:
        #print("Generate moment 0...")
        mom = np.nansum(data, axis = 0)*np.abs(vaxis[1] - vaxis[0])      
        mhd['BUNIT'] = 'K.km/s'
        mom = fits.PrimaryHDU(mom,mhd)
        moms.append(mom)

    
    if 1 in momlist:
        #print("Generate moment 1...")
        mom = np.nansum(data*velcube, axis=0)/np.nansum(data, axis=0)
        mhd['BUNIT'] = 'km/s'
        mom = fits.PrimaryHDU(mom,mhd)
        moms.append(mom)

    
    if 2 in momlist:
        #print("Generate moment 2...")        
        mom = np.sqrt((np.nansum(data*(velcube**2), axis=0)/np.nansum(data, axis=0)) -                    (np.nansum(data*velcube, axis=0)/np.nansum(data, axis=0))**2.)
        mhd['BUNIT'] = 'km/s'
        mom = fits.PrimaryHDU(mom,mhd)
        moms.append(mom)
    return moms

##################################################################################################################################################
 
