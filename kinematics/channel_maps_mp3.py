#!/usr/bin/python

from matplotlib.backends.backend_pdf import PdfPages #multi pdf pages
import matplotlib.pyplot as plt
from astropy.io import ascii
import aplpy
from matplotlib.gridspec import GridSpec #to split page to select area where further plot
from scipy.interpolate import interp1d
import numpy as np #to get numerical function available
import matplotlib.pyplot as plt #to plot
import scipy.constants as c #for Q_rot
from my_moments import *
from coordconverter import * #Pedro K. Humire #to get HMS2deg and deg2HMS functions.
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from astropy.io import fits
from astropy.wcs import WCS
from sigma_clip_for_one_mp import *
import os
import glob
"""
Code to perform channel maps of CO J=2-1 as in Sakamoto+2006
NGC 253 redshift: 0.000864 (SIMBAD)
Splatalogue (CDMS): 230.538, 230.33898712 GHz (redshifted)
B6b: 230.496 GHz (central freq.)

This code make figures centred at SB1, SB1 is to the southwest (Sakamoto et al. 2006)
Now it performs also figures centred at SB2 and, more useful, the Galaxy centre.
"""

def make_channel_maps(cube_name, cenfreq, channels, n_chan_maps, molnam,resol,transition):

        chan_width  = (channels[1]-channels[0])/n_chan_maps
        print('Channel with:', chan_width)

        #This is fixed:
        SB= 'GC' #Super bubble name. If SB is GC, then all the galaxy is shown.
        ALCHEMI_img_radius=0.004
        #ALCHEMI_img_radius=0.012
        fil = 4 #filas
        col = 4 #columnas

        hdu       = fits.open(cube_name)[0]         #read the cube
        ########################################################## Getting header's parameters ################################
        RESTFRQ   = hdu.header['RESTFRQ']
        CRVAL1    = hdu.header['CRVAL1']
        CDELT1    = hdu.header['CDELT1']
        CRPIX1    = hdu.header['CRPIX1']
        CRVAL2    = hdu.header['CRVAL2']
        CDELT2    = hdu.header['CDELT2']
        CRPIX2    = hdu.header['CRPIX2']
        NAXIS1    = hdu.header['NAXIS1'] #total number of pixels in X axis 
        NAXIS2    = hdu.header['NAXIS2'] #total number of pixels in Y axis 
        NAXIS3    = hdu.header['NAXIS3'] #total number of pixels in freq. axis 
        CRVAL3    = hdu.header['CRVAL3'] #first frequency value (Hz)
        CRPIX3    = hdu.header['CRPIX3'] #first pixel value (=1 for ALCHEMI) 
        CDELT3    = hdu.header['CDELT3'] #frequency of each channel bin. 
        BMAJ      = hdu.header['BMAJ']   #maximum beam 
        w = WCS(hdu.header)

        #CDELT3 = -1*np.abs(CDELT3)
        #print('##############################################',CDELT3)

        ####################### Get frequency (Hz) from header and choosen channels: ####################
        #xaxis=(np.arange(NAXIS3)+1 - CRPIX3)*CDELT3+CRVAL3 #original.
        xaxis_nu  = (np.arange(n_chan_maps+1)+1 - CRPIX3)*CDELT3*chan_width+(CRVAL3+CDELT3*(channels[0] + CRPIX3)) #considering,e.g., 100 and 120, it is 20+1 (n_chan_maps+1). I will need it in order to get the intermediate frequency and then the intermediate velocity!

        xaxis_vel = ((cenfreq-xaxis_nu*1E-9)/cenfreq)*c.c.value*1E-3 

        vel_mean = np.zeros([len(xaxis_vel)-1])

        for i in range(0,len(xaxis_vel)-1):
                vel_mean[i]=(xaxis_vel[i]+xaxis_vel[i+1])/2.  #vel_mean is the mean velocity in between each channel. The one that correspond to the right velocity in each moment0 map presented in the channel maps figures.
                
                
        #######Only used when SB= 'GC' is selected#######
        #SB1=HMS2deg(ra='00 47 32.65', dec='-25 17 24.2') ## according to me
        SB1=HMS2deg(ra='00 47 32.53', dec='-25 17 24.0') ##according to Sakamoto+2006
        SB2=HMS2deg(ra='00 47 34.64', dec='-25 17 09.5')
        
        ##according to new data B6 in HNCO 11-10 (myself)
        
        SBa =HMS2deg(ra='00 47 32.9', dec='-25 17 22.5') 
        SBb =HMS2deg(ra='00 47 32.55', dec='-25 17 24.0')
        SBc1=HMS2deg(ra='00 47 32.37', dec='-25 17 23.0')
        SBc2=HMS2deg(ra='00 47 32.45', dec='-25 17 24.5')        
        SBc3=HMS2deg(ra='00 47 32.55', dec='-25 17 25.5')
        
        #From Nanase (see below): bubbles that are not identified in our data                                                                                                                                    
        
        b1 = HMS2deg(ra='00 47 32.455',dec='-25 17 20.749')
        b2 = HMS2deg(ra='00 47 32.981',dec='-25 17 20.493')
        b3 = HMS2deg(ra='00 47 33.585',dec='-25 17 15.808')
        b4 = HMS2deg(ra='00 47 33.928',dec='-25 17 14.747')
        b5 = HMS2deg(ra='00 47 32.8',dec='-25 17 22.0')
        b6 = HMS2deg(ra='00 47 33.2',dec='-25 17 20.0')
        b7 = HMS2deg(ra='00 47 32.2',dec='-25 17 21.0')
        b8 = HMS2deg(ra='00 47 32.5',dec='-25 17 27.0')
        
        #New superbubles (CO 2-1 at 230.338GHz)
        
        #Mickey Mouse (-22 km/s)
        
        MMle = HMS2deg(ra='00 47 32.025',dec='-25 17 26.195') #MM nose
        MMre = HMS2deg(ra='00 47 32.04',dec='-25 17 25.3') #MM right ear      ok
        MMn  = HMS2deg(ra='00 47 31.95',dec='-25 17 25.8') #MM left ear       ok
        MMlc = HMS2deg(ra='00 47 31.97987',dec='-25 17 26.5845') #MM right cheek
        MMrc = HMS2deg(ra='00 47 32.086',dec='-25 17 26.04') #MM left cheek
        MMc  = HMS2deg(ra='00 47 32.05',dec='-25 17 26.9') #MM chin (barbilla) ok
        
        
        """
        Table 3 (Harada+21, Dec 10)
        Coordinates of Analyzed Positions
        Position
        R.
        A.(ICRS)
        Decl.
        (ICRS) Remarks
        00h
        47m −25°17′
        M2 32.29 s 19.10″ Near a superbubble and an SNR
        M3 31.94 29.10″
        M4 32.79 s 21.10″ Near the base of an outflow (SW
        streamer), super hot core
        M5 32.97 s 19.50″ IR Core, in a superbubble, hard X-ray
        source, super hot core
        M7 33.32 s 15.50″ Super hot core
        M8 33.65 s 13.10″
        M9 33.94 s 11.10″ Near a superbubble
        M10 34.15 s 12.30″ Near a superbubble
        A1 32.57 s 27.07″ Near a superbubble
        A2 32.85 s 23.42″ On the outflow SW streamer
        A3 32.96 s 18.53″ Shell of a superbubble
        A4 33.10 s 22.27″ Shell of a superbubble
        A5 33.09 s 19.01″ Shell of a superbubble
        A6 33.26 s 18.05″
        A7 33.52 s 14.51″ Near a superbubble
        A8 33.69 s 14.41″ Near a superbubble
        
        
        Fig2: Superbubble locations are shown with green solid
        circles for ones that were identified in our data, while those that we do not
        identify in our data but were reported in Krieger et al. (2019) are shown with
        dotted green circles with only approximate sizes and locations. 
        
        """
        
        
        #################################################

        if SB == 'SB1':
                GC     =  HMS2deg(ra='00 47 32.65', dec='-25 17 24.2')  #position of SB1 according to me+2020/21
                zoom	= 0.5  #amplification w.r.t. ALCHEMI_img_radius in the figures. More is less ##SB1
                min_p	= 0.02  #minimum percentage w.r.t. maximum value, to be shown, for logarithmic scale in figures  ##SB1
                #channels    = 99,115  #SB1
        if SB == 'SB2':
                GC     =  HMS2deg(ra='00 47 34.64', dec='-25 17 09.5')  #position of SB2 according to Sakamoto+2006, their table 1.
                zoom	= 0.7  #amplification w.r.t. ALCHEMI_img_radius in the figures. More is less ##SB2
                min_p	= 0.05  #minimum percentage w.r.t. maximum value, to be shown, for logarithmic scale in figures ##SB2
                #channels    = 125,141 #SB2
        if SB == 'GC':
                #GC		=HMS2deg(ra='00 47 33.14', dec='-25 17 17.52') (Leroy+2015 footnote (a), Table3, the center)
                #GC		=HMS2deg(ra='00 47 34.04', dec='-25 17 18.52') (Leroy+2015 region 5)
                GC		=HMS2deg(ra='00 47 33.14',dec='-25 17 17.52') #(center assumed in Leroy+2015*, from Muller-Sanchez et al. 2010, see Table 1 in URL below)
                                                                        #https://iopscience.iop.org/article/10.1088/0004-637X/801/1/25/pdf
                #GC		=('11.88838','-25.28817') #Galactic Center. TH2 (Leroy+2017 region 5)
                zoom	= 1.1  #amplification w.r.t. ALCHEMI_img_radius in the figures. More is less ##SB1
                min_p	= 0.04  #CO #minimum percentage w.r.t. maximum value, to be shown, for logarithmic scale in figures  ##SB1
                min_p	= 0.09  #HCN #minimum percentage w.r.t. maximum value, to be shown, for logarithmic scale in figures  ##SB1
                #min_p	= 0.01  #HCN #minimum percentage w.r.t. maximum value, to be shown, for logarithmic scale in figures  ##SB1






        ###########################################################################################

        fig = plt.figure(figsize=(14, 9)) #(x,y) 

        #fig.suptitle('CO J=2-1 Channel maps showing '+str(SB), size=20)                              #Plot total spectrum   
        fig.suptitle(molnam+' J='+transition+' Channel maps at '+resol+'\"', size=20)                              #Plot total spectrum   

        gs = GridSpec(fil, col)		# Number of page's divisions (y,x) coordinates
        gs.update(wspace=0.001, hspace=0.001) # set the spacing between axes. Separacion. Separation.

        fl_rb = cm.rainbow
        fl_rb.set_under(color='blue')


        #Get min max intensity!
        minvalue,maxvalue=[],[]
        for i in range(0,int(n_chan_maps)):
                mommaps1 = make_mommaps(hdu,zpix=([int(channels[1]-(i+1)*chan_width),int(channels[1]-i*chan_width)]))
                minvalue.append(np.max(mommaps1[0].data)*min_p)
                maxvalue.append(np.max(mommaps1[0].data))

        total_min = np.min(minvalue)
        total_max = np.max(maxvalue)


        count=0
        for i in range(0,int(n_chan_maps)):

                if CDELT3 > 0:
                        mommaps1 = make_mommaps(hdu,zpix=([int(channels[1]-(i+1)*chan_width),int(channels[1]-i*chan_width)]))
                else:
                        mommaps1 = make_mommaps(hdu,zpix=([int(channels[0]+(i)*chan_width),int(channels[0]+(i+1)*chan_width)]))


                if i <col:
                        ax1 = aplpy.FITSFigure(mommaps1[0], figure=fig, subplot=list(gs[0,i].get_position(fig).bounds))
                        ax1.recenter(x=float(GC[0]), y=float(GC[1]), radius=ALCHEMI_img_radius*zoom)
                        #ax1.show_colorscale(cmap=fl_rb, vmin=np.max(mommaps1[0].data)*min_p, vmax=np.max(mommaps1[0].data), stretch='power')  #change to logarithmic scale
                        #ax1.show_colorscale(cmap=fl_rb, vmin=total_max*min_p, vmax=total_max, stretch='power')  #change to logarithmic scale
                        ax1.show_colorscale(cmap='nipy_spectral') #linear scale
                        #ax1.show_colorscale(cmap='nipy_spectral') #linear scale
                        #ax1.show_contour(mommaps1[0],colors='black', levels=np.sqrt(2)**np.arange(10)*np.max(mommaps1[0].data)*min_p, linewidths=0.5, returnlevels='yes') 
                        #ax1.add_beam()
                        levels=ax1.show_contour(mommaps1[0],colors='white', levels=np.sqrt(2)**np.arange(12)*np.max(mommaps1[0].data)*min_p, linewidths=0.5, overlap=True, returnlevels=True, alpha=0)
                        ax1.show_contour(mommaps1[0],colors='black', levels=levels[5:len(levels)], linewidths=0.5, overlap=False, returnlevels=True)
                        ax1.set_nan_color('blue')
                        #ax1.beam.set_color('white')
                        #ax1.beam.set_hatch('/')
                        plt.gca().axes.get_xaxis().set_visible(False)
                        aplpy.Ticks.set_tick_direction(ax1, 'in')

                elif i <col*2:
                        ax1 = aplpy.FITSFigure(mommaps1[0], figure=fig, subplot=list(gs[1,i-col*2].get_position(fig).bounds))
                        ax1.recenter(x=float(GC[0]), y=float(GC[1]), radius=ALCHEMI_img_radius*zoom)
                        #ax1.show_colorscale(cmap=fl_rb, vmin=np.max(mommaps1[0].data)*min_p, vmax=np.max(mommaps1[0].data), stretch='power')  #change to logarithmic scale
                        #ax1.show_colorscale(cmap=fl_rb, vmin=total_max*min_p, vmax=total_max, stretch='power')  #change to logarithmic scale
                        ax1.show_colorscale(cmap='nipy_spectral') #linear scale
                        #ax1.show_contour(mommaps1[0],colors='black', levels=np.sqrt(2)**np.arange(10)*np.max(mommaps1[0].data)*min_p, linewidths=0.5)
                        levels=ax1.show_contour(mommaps1[0],colors='white', levels=np.sqrt(2)**np.arange(12)*np.max(mommaps1[0].data)*min_p, linewidths=0.5, overlap=True, returnlevels=True, alpha=0)
                        ax1.show_contour(mommaps1[0],colors='black', levels=levels[4:len(levels)], linewidths=0.5, overlap=False, returnlevels=True)
                        #ax1.show_colorscale(cmap='nipy_spectral') #linear scale
                        #ax1.add_beam()
                        ax1.set_nan_color('blue')
                        #ax1.beam.set_color('white')
                        #ax1.beam.set_hatch('/')
                        plt.gca().axes.get_xaxis().set_visible(False)
                        aplpy.Ticks.set_tick_direction(ax1, 'in')
                elif i <col*3:
                        ax1 = aplpy.FITSFigure(mommaps1[0], figure=fig, subplot=list(gs[2,i-col*3].get_position(fig).bounds))
                        ax1.recenter(x=float(GC[0]), y=float(GC[1]), radius=ALCHEMI_img_radius*zoom)
                        #ax1.show_colorscale(cmap=fl_rb, vmin=np.max(mommaps1[0].data)*min_p, vmax=np.max(mommaps1[0].data), stretch='power')  #change to logarithmic scale
                        #ax1.show_colorscale(cmap=fl_rb, vmin=total_max*min_p, vmax=total_max, stretch='power')  #change to logarithmic scale
                        ax1.show_colorscale(cmap='nipy_spectral') #linear scale
                        #ax1.show_contour(mommaps1[0],colors='black', levels=np.sqrt(2)**np.arange(10)*np.max(mommaps1[0].data)*min_p, linewidths=0.5) 
                        #ax1.add_beam()
                        levels=ax1.show_contour(mommaps1[0],colors='white', levels=np.sqrt(2)**np.arange(12)*np.max(mommaps1[0].data)*min_p, linewidths=0.5, overlap=True, returnlevels=True, alpha=0)
                        ax1.show_contour(mommaps1[0],colors='black', levels=levels[4:len(levels)], linewidths=0.5, overlap=False, returnlevels=True)
                        ax1.set_nan_color('blue')
                        #ax1.beam.set_color('white')
                        #ax1.beam.set_hatch('/')
                        plt.gca().axes.get_xaxis().set_visible(False)
                        aplpy.Ticks.set_tick_direction(ax1, 'in')
                #elif i < int(n_chan_maps)-1:
                else:
                        ax1 = aplpy.FITSFigure(mommaps1[0], figure=fig, subplot=list(gs[3,i-col*4].get_position(fig).bounds))
                        ax1.recenter(x=float(GC[0]), y=float(GC[1]), radius=ALCHEMI_img_radius*zoom)
                        #ax1.show_colorscale(cmap=fl_rb, vmin=np.max(mommaps1[0].data)*min_p, vmax=np.max(mommaps1[0].data), stretch='power')  #change to logarithmic scale
                        #ax1.show_colorscale(cmap=fl_rb, vmin=total_max*min_p, vmax=total_max, stretch='power')  #change to logarithmic scale
                        ax1.show_colorscale(cmap='nipy_spectral') #linear scale
                        #ax1.show_contour(mommaps1[0],colors='black', levels=np.sqrt(2)**np.arange(10)*np.max(mommaps1[0].data)*min_p, linewidths=0.5) 
                        levels=ax1.show_contour(mommaps1[0],colors='white', levels=np.sqrt(2)**np.arange(12)*np.max(mommaps1[0].data)*min_p, linewidths=0.5, overlap=True, returnlevels=True, alpha=0)
                        #print(levels[8:len(levels)],'########################################################')
                        ax1.show_contour(mommaps1[0],colors='black', levels=levels[4:len(levels)], linewidths=0.5, overlap=False, returnlevels=True)
                        if count==0:
                                ax1.add_beam()
                                ax1.beam.set_color('white')
                                ax1.beam.set_hatch('/')
                        ax1.set_nan_color('blue')
                        plt.gca().axes.get_xaxis().set_visible(False)
                        aplpy.Ticks.set_tick_direction(ax1, 'in')
                        count=count+1
                        
                if CDELT3 > 0.0:
                        ax1.add_label(0.86, 0.95, str(int(vel_mean[int(n_chan_maps)-1-i]-258.8))+'km/s ', relative=True, color='white', weight='bold', size='medium') #adds velocity labels.
                else:
                        ax1.add_label(0.86, 0.95, str(int(vel_mean[int(i)]-258.8))+'km/s ', relative=True, color='white', weight='bold', size='medium') #adds velocity labels.
                ax1.show_markers(float(GC[0]),float(GC[1]), marker='*', s=50, c='black', zorder=50) #zorder high number to show it's in front of figures.

                """
                if SB == 'GC':
                    ax1.show_markers(float(SB1[0]),float(SB1[1]), marker='o', s=20, c='red', zorder=50) #from Sakamoto+06
                    ax1.show_markers(float(SB2[0]),float(SB2[1]), marker='o', s=20, c='red', zorder=50) 
                        
                    ax1.show_markers(float(SBa[0]), float(SBa[1]), marker='*', s=20, c='yellow', zorder=50) #from myself (HNCO 11-10)
                    ax1.show_markers(float(SBb[0]), float(SBb[1]), marker='*', s=20, c='yellow', zorder=50)
                    ax1.show_markers(float(SBc1[0]), float(SBc1[1]), marker='*', s=20, c='yellow', zorder=50)
                    ax1.show_markers(float(SBc2[0]), float(SBc2[1]), marker='*', s=20, c='yellow', zorder=50)
                    ax1.show_markers(float(SBc3[0]), float(SBc3[1]), marker='*', s=20, c='yellow', zorder=50)                   
                            
                    ax1.show_markers(float(b1[0]), float(b1[1]), marker='*', s=20, c='white', zorder=50) #Nanase
                    ax1.show_markers(float(b2[0]), float(b2[1]), marker='*', s=20, c='white', zorder=50)
                    ax1.show_markers(float(b3[0]), float(b3[1]), marker='*', s=20, c='white', zorder=50)
                    ax1.show_markers(float(b4[0]), float(b4[1]), marker='*', s=20, c='white', zorder=50)
                    ax1.show_markers(float(b5[0]), float(b5[1]), marker='*', s=20, c='white', zorder=50)
                    ax1.show_markers(float(b6[0]), float(b6[1]), marker='*', s=20, c='white', zorder=50)
                    ax1.show_markers(float(b7[0]), float(b7[1]), marker='*', s=20, c='white', zorder=50)
                    ax1.show_markers(float(b8[0]), float(b8[1]), marker='*', s=20, c='white', zorder=50)   
                """
                    
                ax1.show_markers(float(SB1[0]),float(SB1[1]), marker='o', s=40, c='red', zorder=50) #from Sakamoto+06
                ax1.show_markers(float(SB2[0]),float(SB2[1]), marker='o', s=40, c='red', zorder=50) 
                    
                ax1.show_markers(float(SBa[0]), float(SBa[1]), marker='*', s=40, c='yellow', zorder=50) #from myself (HNCO 11-10)
                ax1.show_markers(float(SBb[0]), float(SBb[1]), marker='*', s=40, c='yellow', zorder=50)
                ax1.show_markers(float(SBc1[0]), float(SBc1[1]), marker='*', s=40, c='yellow', zorder=50)
                ax1.show_markers(float(SBc2[0]), float(SBc2[1]), marker='*', s=40, c='yellow', zorder=50)
                ax1.show_markers(float(SBc3[0]), float(SBc3[1]), marker='*', s=40, c='yellow', zorder=50)                   
                        
                ax1.show_markers(float(b1[0]), float(b1[1]), marker='*', s=40, c='white', zorder=50) #Nanase
                ax1.show_markers(float(b2[0]), float(b2[1]), marker='*', s=40, c='white', zorder=50)
                ax1.show_markers(float(b3[0]), float(b3[1]), marker='*', s=40, c='white', zorder=50)
                ax1.show_markers(float(b4[0]), float(b4[1]), marker='*', s=40, c='white', zorder=50)
                ax1.show_markers(float(b5[0]), float(b5[1]), marker='*', s=40, c='white', zorder=50)
                ax1.show_markers(float(b6[0]), float(b6[1]), marker='*', s=40, c='white', zorder=50)
                ax1.show_markers(float(b7[0]), float(b7[1]), marker='*', s=40, c='white', zorder=50)
                ax1.show_markers(float(b8[0]), float(b8[1]), marker='*', s=40, c='white',zorder=50)  
                
                ax1.show_markers(float(MMle[0]), float(MMle[1]), marker='o', s=60, edgecolor='red', linestyle='--', zorder=50) #Mickey mouse
                ax1.show_markers(float(MMre[0]), float(MMre[1]), marker='o', s=60, edgecolor='red', linestyle='--', zorder=50) #Mickey mouse
                ax1.show_markers(float(MMn[0]), float(MMn[1]), marker='o', s=60, edgecolor='red', linestyle='--', zorder=50) #Mickey mouse
                ax1.show_markers(float(MMlc[0]), float(MMlc[1]), marker='o', s=60, edgecolor='red', linestyle='--', zorder=50) #Mickey mouse
                ax1.show_markers(float(MMrc[0]), float(MMrc[1]), marker='o', s=60, edgecolor='red', linestyle='--', zorder=50) #Mickey mouse
                ax1.show_markers(float(MMc[0]), float(MMc[1]), marker='o', s=60, edgecolor='red', linestyle='--', zorder=50) #Mickey mouse

                        
                """
                MMle = HMS2deg(ra='00 47 32.04',dec='-25 17 25.3') #MM nose
                MMre = HMS2deg(ra='00 47 32.04',dec='-25 17 25.3') #MM right ear
                MMn  = HMS2deg(ra='00 47 31.95',dec='-25 17 25.8') #MM left ear
                MMlc = HMS2deg(ra='00 47 32.06',dec='-25 17 26.3') #MM right cheek
                MMrc = HMS2deg(ra='00 47 32.00',dec='-25 17 26.4') #MM left cheek
                MMc  = HMS2deg(ra='00 47 32.00',dec='-25 17 26.8') #MM chin (barbilla)
                """    

                if (i == 0 or i == col or i == col*2 or i == col*3 or i == col*3+1 or i == col*3+2 or i == col*3+3 or i == col*3+4):
                        plt.gca().coords[0].set_ticks_visible(True)
                else:
                        plt.gca().coords[0].set_axislabel(' ')
                        plt.gca().coords[1].set_axislabel(' ')
                        plt.gca().coords[1].set_ticklabel_visible(False)
                if (i == col*3+1 or i == col*3+2 or i == col*3+3 or i == col*3+4):
                        plt.gca().coords[1].set_axislabel(' ')
                        plt.gca().coords[1].set_ticklabel_visible(False)
                        
                        
                
                if i == int(n_chan_maps)-1:
                        #from mpl_toolkits.axes_grid1 import make_axes_locatable
                        #divider = make_axes_locatable(list(gs[3,i-col*4].get_position(fig).bounds))
                        #cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
                        #fig.colorbar(cm.ScalarMappable(), use_gridspec=True, pad=-0.2)
                        #fig.colorbar(cm.ScalarMappable(), use_gridspec=True, fraction=0.050,pad=0.0)
                        #fig.colorbar(cm.ScalarMappable(), ax=fig.add_axes([0.8, 0.1, 0.03, 0.8]))
                        #ax1.add_colorbar()
                        #ax1.colorbar.set_axis_label_text(r'[Jy beam$^{-1}$ km s$^{-1}$]')
                        #fig, ax = plt.subplots(figsize=(4,4))
                        #divider = make_axes_locatable(ax)
                        #cax = divider.new_vertical(size="5%", pad=0.7, pack_start=True)
                        #fig.add_axes(cax)
                        
                        #cax = fig.add_axes([gs[3,i-col*4].get_position(fig).x1+0.005,gs[3,i-col*4].get_position(fig).y0,0.02,gs[3,i-col*4].get_position(fig).y1-gs[3,i-col*4].get_position(fig).y0])
                        cax = fig.add_axes([gs[3,i-col*4].get_position(fig).x1+0.003,gs[3,i-col*4].get_position(fig).y0,0.01,gs[0,i-col*4].get_position(fig).y1-gs[3,i-col*4].get_position(fig).y0])  #with the last rest I rest the position of the last row so the upper part of the colorbar doesnt continue until the upper part of the figure but until the upper part of the channel map position.
                        #fig.colorbar(cm.ScalarMappable(cmap=fl_rb), use_gridspec=True, pad=-0.2, cax =cax)
                        import matplotlib as mpl
                        norm=mpl.colors.PowerNorm(vmin=total_min,vmax=total_max, gamma=np.sqrt(2))
                        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=fl_rb, norm=norm), boundaries=np.arange(total_min,total_max),cax=cax, format='%.1e')
                        cbar.set_label(r'[Jy beam$^{-1}$ km s$^{-1}$]', rotation=270, labelpad=30, size=20)
                        #plt.colorbar(plt.cm.ScalarMappable(cmap=fl_rb, norm=norm), boundaries=np.linspace(total_min,total_max,10, endpoint=True),cax=cax, format='%.2e')
                        #np.linspace(2.0, 3.0, num=5, endpoint=False)
                        
        return fig


        plt.savefig(str(SB)+'_channel_maps_'+molnam+'_'+res+'_.pdf')



##############################################################################################################################################

#dir_fits='/aux/vbackup2a/backups/phumire/alchemi/Archive/v0/'

#Perform a sigma clip (sc):
sc=5.0

#dir_fits=cwd = os.getcwd()+'/'
dir_fits='/home/pedro/Documentos/Nanase/high_res_cube/cubes/'  #note that the final B6a,b,c...f is not added 
                                                               #because it is automatically added inside the 
                                                               #function sigma_clip_for_one_mp (see below).
                                                               #This makes important to create folders such 
                                                               #as B6e, B6f, from where to put the cubes.
#first run it as:
#maps = {1: {'cube_name' : 'ngc253.B6b.cent.12m7m.230300.contsub.cv0_4.cube.fits',
#Then, depending on the sc level, add:
maps = {1: {'cube_name' : 'ngc253.B6b.cent.12m7m.230300.contsub.cv0_4.cube_sc'+str(sc)+'_.fits',
            'cenfreq': 230.538,
            #'channels': [40,88], This is half fine
            'channels': [96,144],              #this should be a multiple of 16 :/
            'n_chan_maps': 16.,
            'molnam': 'CO',
            'resol': '0.3',
            'transition' :'2-1'}
        }


"""        
maps = {1: {'cube_name' : 'ngc253.B3b.sc4_1.12mC.088300.contsub.cv15.cube_sc3_.fits',
            'cenfreq': 88.6316,
            'channels': [207,239],
            'n_chan_maps': 16.,
            'molnam': 'HCN',
            'resol': '15',
            'transition' :'1-0'},
        2: {'cube_name' : 'ngc253.B5a.7m.177800.contsub.cv12.cube_sc3_.fits',
            'cenfreq': 177.2599,
            'channels': [30,62],
            'n_chan_maps': 16.,
            'molnam': 'HCN',
            'resol': '12?',
            'transition' :'2-1'},
        3: {'cube_name' : 'ngc253.B6h.sc4_1.7m.265700.contsub.cv15.cube_sc3_.fits',
            'cenfreq': 265.8864,
            'channels': [100,132],
            'n_chan_maps': 16.,
            'molnam': 'HCN',
            'resol': '15',
            'transition' :'3-2'},
        4: {'cube_name' : 'ngc253.B7m.sc4_0.7m.354800.contsub.cv15.cube_sc3_.fits',
            'cenfreq': 354.5054,
            'channels': [137,169],
            'n_chan_maps': 16.,
            'molnam': 'HCN',
            'resol': '15',
            'transition' :'4-3'}
        }
"""
cwd = os.getcwd()+'/'



#cubes=(glob.glob(alchemi_data+folder_names+'/*.fits'))

for i in range(1,len(maps)+1):
        if (maps[i]['cube_name'] in glob.glob('*.fits')) == False:
                sigma_clip_for_one_mp(dir_fits,maps[i]['cube_name'][:-5],sc)


#CO J=2-1
#cube_name, cenfreq = 'ngc253.B6b.sc4_0.12m7m.230300.contsub.cv01_6.cube_sc3_.fits', 230.538, [96,144], 16.,'CO','1_6'
#cube_name, cenfreq = 'ngc253.B6b.sc4_0.7m.230300.contsub.cv15.cube_sc3_.fits', 230.538, [96,144], 16.,'CO','15'
#cube_name, cenfreq = 'ngc253.B6b.sc4_0.7m.230300.contsub.cv10.cube_sc3_.fits', 230.538, [96,144], 16.,'CO','10'

#HCN J=1-0
#cube_name,cenfreq,channels,n_chan_maps,molnam,resol,transition = 'ngc253.B3b.sc4_1.12mC.088300.contsub.cv04.cube_sc3_.fits', 88.6316, [214,230], 16., 'HCN','0.4','1-0'




with PdfPages('channel_maps_mp3_sc_'+str(sc)+'.pdf') as pdf:
        for i in range(1,len(maps)+1):
                fig = make_channel_maps(maps[i]['cube_name'], maps[i]['cenfreq'], maps[i]['channels'], maps[i]['n_chan_maps'], maps[i]['molnam'],maps[i]['resol'],maps[i]['transition'])
                pdf.savefig(fig)
                plt.show()
                plt.close()



