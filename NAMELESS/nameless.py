#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 10:10:41 2026

@author: coulombe
"""
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy import interpolate
from scipy.signal import correlate
from scipy.ndimage import median_filter
from scipy import stats
import pickle

import pdb
import time
import glob

from stpipe import datamodel

import NAMELESS.run_jwst_steps as run_jwst

from jwst import datamodels
import pastasoss



#%% NAMELESS class

class nameless(object):
    
    '''
    Data reduction code for JWST NIRISS/SOSS observations 
    '''
    
    # initialization of the nameless reduction object
    def __init__(self,obs_folder=None,red_folder=None,nseg=None,segs=None,handle=None,
                 s1_steps  = ['GS','DQ','SAT','SB','RP','LIN','DC','JD','RF','GAIN'],
                 s2_steps  = ['WCS','FF','ST'],
                 cst_steps = ['BP','BG','CR','1F','EXT']):
        
        
        # directory where the uncal files are stored
        if obs_folder is None:
            raise Exception('You need to provide an observation folder!')
            
        else:
            self.obs_dir = 'Observations/'+obs_folder 
            
        # directory where the reduction outputs will be stored
        if red_folder is None:
            self.red_dir = 'Reductions/'+obs_folder
        
        else:
            self.red_dir = 'Reductions/'+red_folder
            
        # directory where the ancillary files are stored
        self.ancil_dir = 'Auxiliary'
        
        #TODO: make a function that creates folders instead of doing all this in __init__
        # creating reduction folder if it doesn't already exist
        if not os.path.exists(self.red_dir):
            print('Creating folder for reduction outputs')
            os.makedirs(self.red_dir)
            
        # creating stage 1 and 2 output folders if they don't exist
        self.s1_red_dir = self.red_dir+'/stage1_outputs'
        if not os.path.exists(self.s1_red_dir):
            print('Creating folder for stage 1 reduction outputs')
            os.makedirs(self.s1_red_dir)
            
        self.s2_red_dir = self.red_dir+'/stage2_outputs'
        if not os.path.exists(self.s2_red_dir):
            print('Creating folder for stage 2 reduction outputs')
            os.makedirs(self.s2_red_dir)
            
        self.prep_dir = self.red_dir+'/prep_step_outputs'
        if not os.path.exists(self.prep_dir):
            print('Creating folder containing products needed for custom steps')
            os.makedirs(self.prep_dir)
            
        self.cst_dir = self.red_dir+'/custom_steps_outputs'
        if not os.path.exists(self.cst_dir):
            print('Creating folder containing custom steps outputs')
            os.makedirs(self.cst_dir)
            
        self.ext_dir = self.red_dir+'/extraction_outputs'
        if not os.path.exists(self.ext_dir):
            print('Creating folder containing extraction outputs')
            os.makedirs(self.ext_dir)
            
        # getting the handle of the JWST files, typically of format "jwXXXXYYY001_04101_00001"
        if handle is None:
            uncal_files = glob.glob(self.obs_dir+'/*uncal.fits')
            # assuming all files have the same handle, which should be the case
            self.obs_handle = uncal_files[0].split('/')[-1][:25] 
            
        else:
            self.obs_handle = handle
        
        if nseg is None:
            try:
                # number of data segments
                self.nseg = len(uncal_files)
                segs = []
                for i in range(self.nseg):
                    segs.append(uncal_files[i].split('/')[-1][29:32])
            
                self.segs = np.sort(segs)
            except:
                raise Exception('If nseg is set to None, handle also needs to be set to None')
            
        else:
            self.nseg = nseg
            
            # if segs is not provided
            if segs is None:
                segs = []
                for i in range(self.nseg):
                    segs.append((3-len(str(i+1)))*'0'+str(i+1))
                    
                self.segs = np.array(segs)
                
            else:
                self.segs = np.array(segs)
           
        print('segment list: ', self.segs)
            
        # list of stage 1 steps of the jwst STScI pipeline
        self.s1_steps  = s1_steps
        
        # list of stage 2 steps of the jwst STScI pipeline
        self.s2_steps  = s2_steps
        
        # list of custom steps for additional correction
        self.cst_steps = cst_steps
        
        # setting the status of the steps, False = step not completed, True = step completed
        self.s1_status, self.s2_status, self.cst_status = dict(), dict(), dict()
        
        for step in self.s1_steps:
            self.s1_status[step] = False
            
        for step in self.s2_steps:
            self.s2_status[step] = False
            
        for step in self.s1_steps:
            self.cst_status[step] = False
            
        
    def run_reduction(self,run_s1,run_s2,run_prep,run_cst,run_ext,s1_settings,
                      s2_settings,prep_settings,cst_settings,ext_settings):

        # running stage 1 steps of the jwst pipeline
        if run_s1:                
            run_jwst.run_stage1(self,s1_settings)
        
        # running stage 2 steps of the jwst pipeline
        if run_s2:            
            run_jwst.run_stage2(self,s2_settings)
        
        # producing a few outputs to prepare for the NAMELESS custom data reduction steps 
        if run_prep:
            # list of things needed: rough wlc, wavelength solution, time array, median stack, trace profile
            self.prep4cst(prep_settings)
            
        # running custom data correction steps for bad pixels, background, cosmic rays, 1/f, and extraction
        if run_cst:
            # bad pixel correction
            self.badpixcorr(cst_settings)
        
            # background subtraction
            self.bg_corr(cst_settings)
            
            # cosmic ray correction
            self.cosmic_rays(cst_settings)
            
            # 1/f correction
            self.oneoverf(cst_settings)
    
            self.saveclean()
            
        # running the spectral extraction step and saving the light curves
        if run_ext:
            self.extraction(ext_settings)
            
        
    # steps to prepare for the NAMELESS custom reduction steps
    #TODO includes:
    def prep4cst(self,settings):
        
        # let's collate the data from all individual segments into single arrays
        data, err, dq, t_int, pwcpos = self.collate(self.s2_red_dir,'sourcetypestep')
        
        self.data_s2, self.err_s2, self.dq_s2 = data, err, dq 
        self.t_int = t_int
        
        #TODO: save these arrays
        #np.save(self.prep_dir+'/dq_array.npy',self.dq_s2)
        
        self.pwcpos = pwcpos
        
        # creating a DQ mask, where all DQ flags of 0 are considered valid counts
        self.dq_mask = np.zeros(self.dq_s2.shape)
        self.dq_mask[np.where(self.dq_s2>0)] = 1
        
        # computing median stack of the stage 2 data
        self.medstack_s2 = np.nanmedian(self.data_s2,axis=0)
        
        #TODO: Also save median stack and rough white light curve
        
        # saving a figure of the median stack
        if settings['plotting']: #TODO: could even make a figure for plotting 2D images like this, if we repeat it often enough
            fig, ax = plt.subplots(1,1,figsize=(8,5))
            ax.imshow(self.medstack_s2,aspect='auto',vmin=0,vmax=np.nanmax(self.medstack_s2)/10)
            
            ax.set_xlabel('Spectral Direction')
            ax.set_ylabel('Dispersion Direction')
            fig.savefig(self.prep_dir+'/stage_2_median_stack.pdf',dpi=settings['dpi'])
            plt.close(fig)
        
        
        wlc_range = settings['wlc_pixs']
        # considering only a tiny portion of the detector to compute a rough "white" light curve
        cut_data = self.data_s2[:,wlc_range[1][0]:wlc_range[1][1],wlc_range[0][0]:wlc_range[0][1]]
        
        self.rough_wlc = np.nansum(cut_data,axis=(1,2))
        self.rough_wlc/= np.nanmedian(self.rough_wlc)
        
        if settings['plotting']: 
            fig, ax = plt.subplots(1,1,figsize=(8,6))
            ax.errorbar(self.t_int,self.rough_wlc,fmt='o',color='w',mec='k')
            
            ax.set_xlabel('Time [BJD_TDB]')
            ax.set_ylabel('Normalized Counts')
            
            ax2 = ax.twiny()
            ax2.set_xlim(0,len(self.t_int)-1)
            ax2.set_xlabel('Integration Number')
            
            fig.savefig(self.prep_dir+'/stage_2_rough_wlc.pdf',dpi=settings['dpi'])
            plt.close(fig)
            
            
        # next, we want to get the wavelength solution
        self.get_wave_sol(settings)
        
        # and then we want to compute the trace position 
        self.get_trace_pos(settings)


        
    # function that takes multiple jwst .fits files and collates the data, errors, dq flags, 
    # and timestamps into single arrays 
    def collate(self,direc,filetype):
        
        for i,seg in enumerate(self.segs):
            #TODO: assuming for now that this is used to collate files after s2, make it more flexible
            data_file = direc+'/'+self.obs_handle+'-seg'+seg+'_nis_1_'+filetype+'.fits'
            dm = datamodels.open(data_file)
            
            #TODO: keep metadata?
            
            if i==0:
                data, err, dq  = dm.data, dm.err, dm.dq
                t_int  = dm.int_times['int_mid_BJD_TDB'] # getting mid-integration times in BJD_TDB
                
                # NIRISS pupil wheel position
                pwcpos = dm.meta.instrument.pupil_position 
                
            else:
                # concatenating new segment data with existing array
                data, err, dq = (np.concatenate((data,dm.data),axis=0),
                                np.concatenate((err,dm.err),axis=0),
                                np.concatenate((dq,dm.dq),axis=0))
                t_int  = np.concatenate((t_int,dm.int_times['int_mid_BJD_TDB']))
                  
        return data, err, dq, t_int, pwcpos
        #TODO: pwcpos and other important infos could be stored earlier on, to think about

    # function that takes the cleaned data (after custom steps) and puts them back into
    # its constituent segments before saving as jwst datamodels on disk
    def saveclean(self):
        direc = self.s2_red_dir
        filetype = 'sourcetypestep'
        cleandirec = self.cst_dir
        cleanfiletype = 'clean'

        for i, seg in enumerate(self.segs):
            # Read the saved file after s2 to get the segment size and metadata right
            # TODO: assuming for now that this is used to collate files after s2, make it more flexible
            data_file = direc + '/' + self.obs_handle + '-seg' + seg + '_nis_1_' + filetype + '.fits'
            dm = datamodels.open(data_file)
            # the starting and ending integrations for that segment
            istart = dm.meta.exposure.integration_start
            iend = dm.meta.exposure.integration_end
            # replace the datamodel content by the values from cst
            dm.data = np.copy(self.data_cst[istart-1:iend])
            dm.err = np.copy(self.err_cst[istart-1:iend])
            dm.dq = np.copy(self.dq_cst[istart-1:iend])
            # save on disk in the cst directory with clean file name
            dm.save(cleandirec + '/' + self.obs_handle + '-seg' + seg + '_nis_1_' + cleanfiletype + '.fits')

        return

    
    # function that calls pastasoss to get wavelength solution and trace position
    def get_wave_sol(self,settings):
        
        # getting trace position and wavelenths from pastasoss for orders 1 and 2
        trace_o1 = pastasoss.get_soss_traces(pwcpos=self.pwcpos, order='1', interp=True)
        trace_o2 = pastasoss.get_soss_traces(pwcpos=self.pwcpos, order='2', interp=True)
        
        # wavelengths
        self.wave_o1, self.wave_o2 = trace_o1.wavelength[1:-1], trace_o2.wavelength
        
        # trace positions
        self.pasta_trace1 = np.array([trace_o1.x,trace_o1.y])[:,1:-1]
        self.pasta_trace2 = np.array([trace_o2.x,trace_o2.y])
    
        # applying lower cutoff on order 2 since flux becomes neglibile at some point
        # xmin_o2 corresponds to cutting off wavelengths ~> 1.05, which are redundant with order 1
        self.wave_o2      = self.wave_o2[np.where((trace_o2.x>=settings['xmin_o2'])&(trace_o2.x<=settings['xmax_o2']))]
        self.pasta_trace2 = self.pasta_trace2[:,np.where((trace_o2.x>=settings['xmin_o2'])&(trace_o2.x<=settings['xmax_o2']))[0]]
    
    
    def get_trace_pos(self,settings): #method='CCF',sampfac=100,nknots1=11,nknots2=7,plotting=True,dpi=300):
        
        method = settings['trace_method']
        
        # simply using the values returned by pastasoss for the trace positions
        if method=='pastasoss':
            self.trace1 = self.pasta_trace1
            self.trace2 = self.pasta_trace2
            
            y_o1, y_o1_smooth = self.trace1[1], self.trace1[1]
            y_o2, y_o2_smooth = self.trace2[1], self.trace2[1]
            
        # doing cross-correlation with a reference trace profile
        elif method=='CCF': #TODO: this should be a separate function
            SP = fits.open(self.ancil_dir+'/jwst_niriss_specprofile_0022.fits')
            ST = fits.open(self.ancil_dir+'/jwst_niriss_spectrace_0023.fits')
            
            sampfac = settings['sampfac']
            nknots1, nknots2 = settings['nknots_o1'], settings['nknots_o2']
            
            # model order 1 trace profile
            trace_prof1 = SP[1].data[20:-20,20:-20]
            
            # taking one column from the model trace as a reference
            ind_ref = 2038 # column # 2038, could be any of them
            mod_col = trace_prof1[:,ind_ref] 
            abs_pos = ST[1].data['Y'][np.argmin(np.abs(ST[1].data['X']-ind_ref-0.5))]
            #TODO: could also use each column of the model trace for each column of the detector
            
            y      = np.arange(len(mod_col)) + 0.5
            y_fine = np.linspace(0,len(mod_col),len(mod_col)*sampfac)
            
            int_mod_col  = interpolate.interp1d(y,mod_col,kind='cubic',fill_value='extrapolate')
            mod_col_fine = int_mod_col(y_fine)
            
            # y positions of order 1 trace
            y_o1 = np.zeros(self.wave_o1.size)
            y_cc = np.linspace(-len(mod_col_fine),len(mod_col_fine),2*len(mod_col_fine)-1)/sampfac
            
            # computing maximum of cross correlation function for each column of order 1
            for i,x in enumerate(self.pasta_trace1[0]):
                data_col      = np.copy(self.medstack_s2[:,int(x)])
                
                # let's set everything that's >30 pixels away from the pastasoss guess to 0
                data_col[0:int(self.pasta_trace1[1,i]-30)]  = np.nanmedian(self.medstack_s2[:,int(x)])
                data_col[int(self.pasta_trace1[1,i]+30):-1] = np.nanmedian(self.medstack_s2[:,int(x)])
                
                int_data = interpolate.interp1d(y[~np.isnan(data_col)],data_col[~np.isnan(data_col)],
                                                kind='cubic',fill_value='extrapolate')
                
                data_col_fine = int_data(y_fine)
                cc = correlate(data_col_fine,mod_col_fine)
                
                y_o1[i] = y_cc[np.argmax(cc)]
                
            y_o1 += abs_pos
            
            
            # smoothing the profile, adapted from transitspectroscopy package routine
            #TODO: define xvals from pasta_trace1, cleaner
            xend, xstart = self.pasta_trace1[0,0], self.pasta_trace1[0,-1]
            
            knots = np.arange(xend,xstart, (xstart - xend) / np.double(nknots1))[1:]
            tck   = interpolate.splrep(np.double(self.pasta_trace1[0]),y_o1, t=knots)
            
            y_o1_smooth = interpolate.splev(self.pasta_trace1[0], tck)
                
            
            # repeating the same thing for order 2
            y_o2 = np.zeros(self.wave_o2.size)
            
            for i,x in enumerate(self.pasta_trace2[0]):
                data_col = np.copy(self.medstack_s2[:,int(x)])
                
                data_col[0:int(self.pasta_trace2[1,i]-30)]  = np.nanmedian(self.medstack_s2[:,int(x)])
                data_col[int(self.pasta_trace2[1,i]+30):-1] = np.nanmedian(self.medstack_s2[:,int(x)])
                
                data_col[-5:] = np.nanmedian(data_col)
                
                int_data = interpolate.interp1d(y[~np.isnan(data_col)],data_col[~np.isnan(data_col)],
                                                kind='cubic',fill_value='extrapolate')
                
                data_col_fine = int_data(y_fine)
                cc = correlate(data_col_fine,mod_col_fine)
                
                y_o2[i] = y_cc[np.argmax(cc)]
                
                
            y_o2 += abs_pos
            
            xend, xstart = self.pasta_trace2[0,0], self.pasta_trace2[0,-1]
            
            knots = np.arange(xend,xstart, (xstart - xend) / np.double(nknots2))[1:]
            tck   = interpolate.splrep(np.double(self.pasta_trace2[0]),y_o2, t=knots)
            
            y_o2_smooth = interpolate.splev(self.pasta_trace2[0], tck)           
        
            
            #TODO: again, set arrays of x1 and x2 instead of always calling pasta_trace
            self.trace1 = np.array([self.pasta_trace1[0],y_o1_smooth])
            self.trace2 = np.array([self.pasta_trace2[0],y_o2_smooth])
            
            #TODO: also save these arrays as .txt files
            SP.close()
            ST.close()
            
        else:
            #TODO: also add convolution method
            print('Other methods have not yet been implemented!')
            pdb.set_trace()
            
            
        #TODO: put at the end of the function
        if settings['plotting']:
            fig, ax = plt.subplots(1,1,figsize=(8,6))

            ax.imshow(self.medstack_s2,aspect='auto',vmin=0.,vmax=np.nanmax(self.medstack_s2/10))
            
            ax.plot(self.pasta_trace1[0],y_o1,'-r',lw=2)
            ax.plot(self.pasta_trace1[0],y_o1_smooth,'--b',lw=2)
            
            ax.plot(self.pasta_trace2[0],y_o2,'-r',lw=2)
            ax.plot(self.pasta_trace2[0],y_o2_smooth,'--b',lw=2)
            fig.savefig(self.prep_dir+'/trace_positions.pdf',dpi=settings['dpi'])
            plt.close(fig)
    
    #TODO: all the data reduction functions should be callable without using the nameless function
    
    def badpixcorr(self,settings):
        
        data, dq_mask = self.data_s2, self.dq_mask
        
        med_data = np.nanmedian(data,axis=0)
        
        med_dq = np.zeros(self.dq_mask.shape[1:])
        med_dq[np.mean(dq_mask,axis=0)>settings['frac_int']] = 1
        
        data_corr = np.zeros(data.shape)
        
        #TO-DO: save outlier map (med_dq) for reference
        for i in range(len(data)):
            #TODO: take out step counting, or change for something that looks better
            print('\rBad Pixel Step: '+str(i+1)+'/'+str(len(data)),end='')
            
            # integration to be corrected
            corr_int = np.copy(data[i])
            
            inds_gp    = np.where((med_dq==0) & (~np.isnan(corr_int))) # indices of good pixels
            y_gp, x_gp = inds_gp
            z_gp       = corr_int[inds_gp] # values of good pixels
            
            inds_bp    = np.where((med_dq==1) | (np.isnan(corr_int))) # indices of bad pixels
            y_bp, x_bp = inds_bp
            
            #TODO: check if this is needed
            #corr_int[inds_bp] = -1e3 # set bad pixels to some arbitary non-NaN value
            
            corr_int[inds_bp] = interpolate.griddata((x_gp,y_gp),z_gp,(x_bp,y_bp),method=settings['interp_bp'],fill_value=0.)
            
            data_corr[i] = corr_int
            
        #TODO: add option to set flags of corrected pixels back to 0, if we trust them
        # would also require to look how to fix the errors as well
        
        if settings['plotting']:
            fig, axes = plt.subplots(1,2,figsize=(12,5))
            
            axes[0].imshow(med_data,aspect='auto',vmin=0,vmax=np.nanmax(med_data)/10)
            axes[1].imshow(np.nanmedian(data_corr,axis=0),aspect='auto',vmin=0,vmax=np.nanmax(med_data)/10)
            
            for ax in axes:
                ax.set_xlabel('Spectral Direction')
                ax.set_ylabel('Dispersion Direction')
                
            fig.savefig(self.cst_dir+'/bad_pixel_correction.pdf',dpi=settings['dpi'])
            plt.close(fig)
        
        self.data_cst, self.err_cst, self.dq_cst = data_corr, self.err_s2, self.dq_s2
    
    
    def bg_corr(self,settings):
        
        #TODO: add option for 96 substrip model
        bg_mod = np.load(self.ancil_dir+'/model_background256.npy')
        
        #TODO: add option to not consider the pixels with DQ flags
        med_data = np.nanmedian(self.data_cst,axis=0)
        
        grad_bg   = np.gradient(bg_mod,axis=1) # finding the position of the background jump for each row
        grad_bg[:,10:-10] = 0. # setting the edge derivatives to 0
        
        pos_jump = np.argmax(grad_bg,axis=1) 
        
        inds1, inds2 = settings['inds1_bg'], settings['inds2_bg']
        
        xmin1,xmax1,ymin1,ymax1 = inds1[0][0],inds1[0][1],inds1[1][0],inds1[1][1]
        xmin2,xmax2,ymin2,ymax2 = inds2[0][0],inds2[0][1],inds2[1][0],inds2[1][1]
        
        
        bg_scaling1 = np.quantile(med_data[ymin1:ymax1,xmin1:xmax1]/bg_mod[ymin1:ymax1,xmin1:xmax1],settings['bg_quant'])
        bg_scaling2 = np.quantile(med_data[ymin2:ymax2,xmin2:xmax2]/bg_mod[ymin2:ymax2,xmin2:xmax2],settings['bg_quant'])
        
        scaled_bg = np.zeros(bg_mod.shape)
        
        for i in range(len(bg_mod)):
            scaled_bg[i,:pos_jump[i]] = bg_mod[i,:pos_jump[i]] * bg_scaling1
            scaled_bg[i,pos_jump[i]:] = bg_mod[i,pos_jump[i]:] * bg_scaling2
            
        
        #TODO: add plotting option for a single row
        self.data_cst -= scaled_bg[np.newaxis]
        
        if settings['plotting']:
            fig, axes = plt.subplots(1,2,figsize=(12,5))
            
            axes[0].imshow(med_data,aspect='auto',vmin=0,vmax=np.nanmax(med_data)/10)
            axes[1].imshow(np.nanmedian(self.data_cst,axis=0),aspect='auto',vmin=0,vmax=np.nanmax(med_data)/10)
            
            for ax in axes:
                ax.set_xlabel('Spectral Direction')
                ax.set_ylabel('Dispersion Direction')
                
            fig.savefig(self.cst_dir+'/background_subtraction.pdf',dpi=settings['dpi'])
            plt.close(fig)
        
    
    
    def cosmic_rays(self,settings):
        
        # let's do an initial 1/f correction taking only the median of each frame
        med_data = np.nanmedian(self.data_cst,axis=0)
        
        wlc_range = settings['wlc_pixs']
        rough_wlc = np.nansum(self.data_cst[:,wlc_range[1][0]:wlc_range[1][1],wlc_range[0][0]:wlc_range[0][1]],axis=(1,2))
        rough_wlc = rough_wlc / np.nanmedian(rough_wlc)
        
        rough_wlc = median_filter(rough_wlc,settings['wlc_window'])
        
        astro_sig = med_data[np.newaxis] * rough_wlc[:,np.newaxis,np.newaxis]
        
        old_data = np.copy(self.data_cst)
        
        data = self.data_cst - astro_sig
        
        med_col = np.nanmedian(data,axis=1) # 1/f noise
        
        data -= med_col[:,np.newaxis]
        
        running_median = median_filter(data,(settings['CR_window'],1,1))
        
        resids = data - running_median
        
        resids[:,:5],resids[:,-5:]  = 0.,0. # not considering the reference pixels 
        
        stdev = stats.iqr(resids,rng=(16,84),axis=0)/2.
        
        #TODO: also save array that corresponds to the CR flags, and update DQs
        nsig = settings['nsig_CR']
        data[np.where(np.abs(resids)>nsig*stdev)] = running_median[np.where(np.abs(resids)>nsig*stdev)]

        # now let's add back the 1/f and astrophysical signal
        
        data += med_col[:,np.newaxis]
        data += astro_sig
        
        self.data_cst = data
        
        if settings['plotting']:
            inds = np.array([1,2,3]) * len(data) / 3. - 1
            inds = inds.astype(int)
            
            fig, axes = plt.subplots(3,2,figsize=(12,12))
            
            for i,ind in enumerate(inds):
            
                axes[i,0].imshow(old_data[ind],aspect='auto',vmin=0,vmax=np.nanmax(data[ind])/10)
                axes[i,1].imshow(data[ind],aspect='auto',vmin=0,vmax=np.nanmax(data[ind])/10)
            
                for ax in axes[i]:
                    ax.set_xlabel('Spectral Direction')
                    ax.set_ylabel('Dispersion Direction')
                
            fig.savefig(self.cst_dir+'/CR_correction.pdf',dpi=settings['dpi'])
            plt.close(fig)
            


    def oneoverf(self,settings):
        
        #TODO: clean up this function

        #TODO: add option to consider a different median stack as a function of time (in case of tilt event)
        
        median_ints = np.arange(settings['median_ints'][0],settings['median_ints'][1])

        med_stack = np.nanmedian(self.data_cst[median_ints],axis=0)[np.newaxis]

        #TODO: fix so that we also have err_cst and dq_cst values
        data, err, dq = np.copy(self.data_cst), np.copy(self.err_s2), np.copy(self.dq_s2)

        err[dq!=0]         = np.inf # setting all errors of pixels with non zero dq flags to infinity
        err[np.isnan(err)] = np.inf # setting error values that are nan or 0 to infinity
        err[err==0]        = np.inf

        err2 = np.copy(err)

        # setting error of pixels >ap from the order 1 center to infinity
        x1,y1 = self.trace1
        
        ap = settings['ap_1f']
        
        for i in range(x1.size):
            err[:,int(y1[i]+ap):,int(x1[i])] = np.inf
            if not (y1[i]-ap)<0:
                err[:,:int(y1[i]-ap),int(x1[i])] = np.inf
            
        
        # computing the values of scaling and 1/f noise of the median stack that minimizes the chi2
        ssx  = np.nansum(med_stack/(err**2.),axis=1)
        ssx2 = np.nansum((med_stack/err)**2.,axis=1)
        ssy  = np.nansum(data/(err**2.),axis=1)
        ssxy = np.nansum(med_stack*data/(err**2.),axis=1)
        ssw  = np.nansum(1/(err**2.),axis=1)
        
        ssxssw = ssx/ssw
        ssyssw = ssy/ssw
        
        slope     = (ssxssw*ssy - ssxy) / (ssxssw*ssx - ssx2) # scaling   
        intercept = ssyssw - slope*ssxssw # 1/f noise
        
        
        # and we now repeat this for order 2! 
        data_corr1 = data - intercept[:,np.newaxis] # data corrected for first round of 1/f 
        
        x2,y2 = self.trace2
        
        mask1 = np.zeros(data.shape[1:])
        mask2 = np.zeros(data.shape[1:])
        
        for i in range(x1.size):
            mask1[int(y1[i]+ap):,int(x1[i])] = 1
            if not (y1[i]-ap)<0:
                mask1[:int(y1[i]-ap),int(x1[i])] = 1
                
        ind2 = 0
        for i in range(x1.size):
            if not np.any(x1[i]==x2):
                mask2[:,int(x1[i])] = 1
                mask2[:,int(x1[i])] = 1

            else:
                mask2[int(y2[ind2]+ap):,int(x2[ind2])] = 1
                mask2[:int(y2[ind2]-ap),int(x2[ind2])] = 1
                
                ind2+=1 
                
                
        mask2[np.where((mask1+mask2) == 0)] = 1
        mask2[np.where(mask2==1)] = np.inf
        
        err2 *= (1 + mask2[np.newaxis])
        
        ssx_2  = np.nansum(med_stack/(err2**2.),axis=1)
        ssx2_2 = np.nansum((med_stack/err2)**2.,axis=1)
        ssy_2  = np.nansum(data_corr1/(err2**2.),axis=1)
        ssxy_2 = np.nansum(med_stack*data_corr1/(err2**2.),axis=1)
        ssw_2  = np.nansum(1/(err2**2.),axis=1)
        
        ssxssw_2 = ssx_2/ssw_2
        ssyssw_2 = ssy_2/ssw_2
        
        slope2     = (ssxssw_2*ssy_2 - ssxy_2) / (ssxssw_2*ssx_2 - ssx2_2)
        intercept2 = ssyssw_2 - slope2*ssxssw_2
        
        
        if settings['plotting']:
            fig,ax = plt.subplots(2,2,figsize=(12,8))
            
            ax[0,0].imshow(slope,aspect='auto',vmin=0.998,vmax=1.002)
            ax[0,1].imshow(slope2,aspect='auto',vmin=0.998,vmax=1.002)
            
            ax[1,0].imshow(intercept,aspect='auto',vmin=-0.5,vmax=0.5)
            ax[1,1].imshow(intercept2,aspect='auto',vmin=-0.5,vmax=0.5)
            
            for i in range(2):
                ax[0,i].set_xlabel('Column Number')
                ax[1,i].set_xlabel('Column Number')
                
                ax[0,i].set_ylabel('Integration Number')
                ax[1,i].set_ylabel('Integration Number')
            
            fig.savefig(self.cst_dir+'/scaling_and_1overf_noise.pdf',dpi=settings['dpi'])
            plt.close(fig)
        
        
        intercept2[np.where(np.isnan(intercept2))] = 0.
        mask2[np.where(mask2 == np.inf)] = 1.
        
        resid_frames = data - (med_stack * slope[:,np.newaxis])        

        # multiply resid_frames by a mask of the trace: with nans in it
        resid_corrs  = resid_frames - intercept[:,np.newaxis] - intercept2[:,np.newaxis] * np.abs(mask2[np.newaxis]-1.)
        
        frames_corr  = resid_corrs + (med_stack * slope[:,np.newaxis])
        
        data_corr   = np.copy(data)
        
        data_corr  -= (intercept[:,np.newaxis] + intercept2[:,np.newaxis] * np.abs(mask2[np.newaxis]-1.))
        
        
        self.data_cst = data_corr
        
        if settings['plotting']:
            fig,axes = plt.subplots(5,2,figsize=(12,16),sharex=True)
        
            int_no = [int(len(data)/3),int(2*len(data)/3)]
            
            im0 = axes[0,0].imshow(data[int_no[0]],aspect='auto',vmin=0,vmax=10)
            im0 = axes[0,1].imshow(data[int_no[1]],aspect='auto',vmin=0,vmax=10)
            clb = plt.colorbar(im0,ax=axes[0])
            clb.set_label('Rate/Int',fontsize=8)
            axes[0,0].set_title('1st Integration')
            axes[0,1].set_title('Last Integration')
            
            im1 = axes[1,0].imshow(med_stack[0],aspect='auto',vmin=0,vmax=10)
            im1 = axes[1,1].imshow(med_stack[0],aspect='auto',vmin=0,vmax=10)
            clb = plt.colorbar(im1,ax=axes[1])
            clb.set_label('Median Integration',fontsize=8)
            
            im2 = axes[2,0].imshow(resid_frames[int_no[0]],aspect='auto',vmin=-1,vmax=1)
            im2 = axes[2,1].imshow(resid_frames[int_no[1]],aspect='auto',vmin=-1,vmax=1)
            clb = plt.colorbar(im2,ax=axes[2])
            clb.set_label('Integration - Median',fontsize=8)
        
            im3 = axes[3,0].imshow(resid_corrs[int_no[0]],aspect='auto',vmin=-1,vmax=1)
            im3 = axes[3,1].imshow(resid_corrs[int_no[1]],aspect='auto',vmin=-1,vmax=1)
            clb = plt.colorbar(im3,ax=axes[3])
            clb.set_label('C-by-C Subtraction',fontsize=8)
            
            im4 = axes[4,0].imshow(frames_corr[int_no[0]],aspect='auto',vmin=0,vmax=10)
            im4 = axes[4,1].imshow(frames_corr[int_no[1]],aspect='auto',vmin=0,vmax=10)
            clb = plt.colorbar(im4,ax=axes[4])
            clb.set_label('Median + C-by-C Corrected',fontsize=8)
            
            fig.savefig(self.cst_dir+'/1overf_noise_correction.pdf',dpi=settings['dpi'])
            print('Produced and saved figure')
            plt.close()
        

    
    def extraction(self,settings):
        
        data, err, dq = np.copy(self.data_cst), np.copy(self.err_s2), np.copy(self.dq_s2)
        
        #TODO: include option to not consider pixels with dq flags in extraction
        
        # extraction of order 1 data
        x1, y1 = self.trace1[0], self.trace1[1]
        
        aperture_radius = np.zeros(x1.size) + settings['ap_radius']
        
        flux_arr1, err_arr1 = self.aperture_extraction_vr(data,
                                                          x1,
                                                          y1,
                                                          aperture_radius = aperture_radius,
                                                          error_data = err)
        
        
        
        # repeating extraction for order 2
        x2, y2 = self.trace2[0], self.trace2[1]        
        
        aperture_radius = np.zeros(x2.size) + settings['ap_radius']

        flux_arr2, err_arr2 = self.aperture_extraction_vr(data,
                                                          x2,
                                                          y2,
                                                          aperture_radius = aperture_radius,
                                                          error_data = err)
        
        # instead of using the errors from the extraction, let's compute the point-to-point scatter
        noise_arr1 = 0.5*(flux_arr1[0:-2]+flux_arr1[2:]) - flux_arr1[1:-1]
        noise_arr1 = (stats.iqr(noise_arr1,axis=0,rng=(16,84))/2) / np.sqrt(3.0/2.0)
        noise_arr1 = np.zeros(flux_arr1.shape) + noise_arr1[np.newaxis] 
        
        noise_arr2 = 0.5*(flux_arr2[0:-2]+flux_arr2[2:]) - flux_arr2[1:-1]
        noise_arr2 = (stats.iqr(noise_arr2,axis=0,rng=(16,84))/2) / np.sqrt(3.0/2.0)
        noise_arr2 = np.zeros(flux_arr2.shape) + noise_arr2[np.newaxis] 
        
        #TODO: Do light curve figures, and also save in xarray format
        
        primary_hdu = fits.PrimaryHDU()
        
        # Order 1
        hdu_f1 = fits.ImageHDU(flux_arr1, name='FLUX_O1')
        hdu_e1 = fits.ImageHDU(noise_arr1,  name='ERR_O1')
        
        hdu_fn1 = fits.ImageHDU(flux_arr1/np.nanmedian(flux_arr1,axis=0), name='NORM_FLUX_O1')
        hdu_en1 = fits.ImageHDU(noise_arr1/np.nanmedian(flux_arr1,axis=0),  name='NORM_ERR_O1')
        
        hdu_w1 = fits.ImageHDU(self.wave_o1, name='WAVE_O1')

        
        # Order 2 
        hdu_f2 = fits.ImageHDU(flux_arr2, name='FLUX_O2')
        hdu_e2 = fits.ImageHDU(noise_arr2,  name='ERR_O2')
        
        hdu_fn2 = fits.ImageHDU(flux_arr2/np.nanmedian(flux_arr2,axis=0), name='NORM_FLUX_O2')
        hdu_en2 = fits.ImageHDU(noise_arr2/np.nanmedian(flux_arr2,axis=0),  name='NORM_ERR_O2')
        
        hdu_w2 = fits.ImageHDU(self.wave_o2, name='WAVE_O2')
        
        
        # Time 
        hdu_t = fits.ImageHDU(self.t_int, name='TIME')
        
        hdul = fits.HDUList([primary_hdu, hdu_f1, hdu_e1, hdu_fn1, hdu_en1, hdu_w1,
                             hdu_f2, hdu_e2, hdu_fn2, hdu_en2, hdu_w2, hdu_t])
        
        hdul.writeto(self.ext_dir+'/2D_light_curves_extraction.fits',overwrite=True)
        
        
   
    
    # aperture extraction function slightly adapted from N. Espinoza's transitspectroscopy pipeline
    def aperture_extraction_vr(self,data, x, y, aperture_radius, error_data=None):
        """
        This function takes as inputs two arrays (x,y) that follow the trace,
        and returns the added flux over the defined aperture radius (and its error, if an error matrix
        is given as well), substracting in the way any background between the aperture radius and the
        background radius. The background is calculated by taking the median of the points between the
        aperture_radius and the background_radius.
        Parameters
        ----------
        data: ndarray
            Image from which the spectrum wants to be extracted
        x: ndarray
            Array with the x-axis of the trace (i.e., the columns, wavelength direction)
        y: ndarray
            Array with the y-axis of the trace (i.e., rows, spatial direction)
        aperture_radius: float
            Distance from the center of the trace at which you want to add fluxes.
        error_data: ndarray
            Image with the errors of each pixel value on the data ndarray above
        """
        
        # Create array that will save our fluxes:
        flux = np.zeros((len(data),len(x)))
        if error_data is not None:
            flux_error = np.zeros((len(data),len(x)))
            
        max_column = data.shape[1] - 1
    
        for i in range(len(x)):
            # Cut the column with which we'll be working with:
            column = data[:,:,int(x[i])]
            if error_data is not None:
                variance_column = error_data[:,:,int(x[i])]**2
                
            # Define limits given by the aperture_radius and background_radius variables:                
            left_side_ap  = np.max([y[i] - aperture_radius[i], 0])
            right_side_ap = np.min([max_column, y[i] + aperture_radius[i]])
    
    
            # Perform aperture extraction of the background-substracted column, being careful with pixelization
            # at the edges. First, deal with left (up) side:
            l_decimal, l_integer = np.modf(left_side_ap)
            l_integer = int(l_integer)
            if l_decimal < 0.5:
                l_fraction = (0.5 - l_decimal) * column[:,np.min([l_integer, max_column])]
                l_limit = l_integer + 1
                if error_data is not None:
                    l_fraction_variance = ((0.5 - l_decimal)**2) * variance_column[:,np.min([l_integer, max_column])]
            else:
                l_fraction = (1. - (l_decimal - 0.5)) * column[:,np.min([l_integer + 1, max_column])]
                l_limit = l_integer + 2
                if error_data is not None:
                    l_fraction_variance = ((1. - (l_decimal - 0.5))**2) * variance_column[:,np.min([l_integer + 1, max_column])]
    
            # Now right (down) side:
            r_decimal, r_integer = np.modf(right_side_ap)
            r_integer = int(r_integer)
            if r_decimal < 0.5:
                r_fraction = (1. - (0.5 - r_decimal)) * column[:,np.min([max_column, r_integer])]
                r_limit = r_integer
                if error_data is not None:
                    r_fraction_variance = ((1. - (0.5 - r_decimal))**2) * variance_column[:,np.min([max_column, r_integer])]
            else:
                r_fraction = (r_decimal - 0.5) * column[:,np.min([max_column, r_integer + 1])]
                r_limit = r_integer + 1
                if error_data is not None:
                    r_fraction_variance = ((r_decimal - 0.5)**2) * variance_column[:,np.min([max_column, r_integer + 1])]
    
            # Save total flux in current column:
            flux[:,i] = l_fraction + r_fraction + np.sum(column[:,l_limit:r_limit],axis=1)

            if error_data is not None:
                # Total error is the sum of the variances:
                flux_error[:,i] = np.sqrt(np.sum(variance_column[:,l_limit:r_limit],axis=1) + l_fraction_variance + \
                                        r_fraction_variance)
        if error_data is not None:
            return flux, flux_error
        else:
            return flux
        
    # add loading and saving functions for the nameless object
    
    # add create folders function 
    
    # have a function to do PCA