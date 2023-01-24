#!/usr/bin/env python
# -*- coding: utf-8 -*-   #con esto puedo poner tilde.

import numpy as np
from astropy.io import fits

"""
Code to perform a sigma clip to a one selected file.
"""


def snr_mask(hdu, fchan, snr):
    
    """
    Generate a mask based on the SNR
    per spectrum across the datacube
    
    hdu: astropy.hdu instance
    
    fchan: int
    Number of free-line-channels to
    at the beginning and end of the
    datacube to consider to generate
    a RMS map
    
    snr: float
    Sigma clipping SNR
    """
    
    data = hdu.data
    hd = hdu.header

    free1 = data[0:fchan,:,:]					#cut the 10 pix in the edges
    free2 = data[data.shape[0]-fchan:-1,:,:]
    freeT = np.concatenate((free1,free2), axis=0) 
    rmsmap = np.nanstd(freeT,axis=0)

    rmsmap[rmsmap == 0] = np.nan
    s2ncube = np.zeros(data.shape)

    for v in range(data.shape[0]):
        s2ncube[v,:,:] = data[v,:,:]/rmsmap[:,:]

    mask = np.zeros(data.shape, dtype='int32')
    mask[s2ncube >= snr] = 1
	
    return mask

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


#for ALCHEMI:

#dir_fits='/media/pedro/DATA_HDD/Doctorado/alchemi/Archive/v0/'


#dir_fits='/aux/vbackup2a/backups/phumire/alchemi/Archive/v0/B3b/'
#fname='ngc253.B3b.sc4_1.12mC12mE.088300.contsub.cv01_6.cube'
#fname='ngc253.B3b.sc4_1.12mC.088300.contsub.cv15.cube'
#fname='ngc253.B3b.sc4_1.12mC.088300.contsub.cv04.cube'


def sigma_clip_for_one_mp(dir_fits,fname,sc):														#sigma clip level
	# Open datacube and 
	hdu = fits.open(dir_fits+fname[7:10]+'/'+fname+'.fits')[0] #fname[7:10] is always the alma Band in ALCHEMI's nomenclature.
	#hdu = fits.open(dir_fits+fname+'.fits')[0] #fname[7:10] is always the alma Band in ALCHEMI's nomenclature.
	hdu.data = hdu.data.squeeze()  #drop fourth dimension (polarization axis)
	hdu.header = remove_hd_dim(hdu.header, dim=3)
	# Sigma clipping:
	mask1 = snr_mask(hdu, 4, sc)
	hdu.data[mask1 == 0] = np.nan
	hdu.writeto(fname+'_sc'+str(sc)+'_.fits')





