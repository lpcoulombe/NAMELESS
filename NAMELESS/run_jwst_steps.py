#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 14:04:24 2025

@author: coulombe
"""


from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_spec2
from jwst import datamodels

import pdb

#%%

# ====== Stage 1 steps ======

def run_stage1(self,settings):

    
    for seg in self.segs:
        uncal_file = self.obs_dir+'/'+self.obs_handle+'-seg'+seg+'_nis_uncal.fits'
        
        # can add option to save intermediate steps at some point, also to override reference fils
        
        # group scale step
        results = calwebb_detector1.group_scale_step.GroupScaleStep.call(uncal_file,
                                                                         output_dir=self.s1_red_dir,
                                                                         save_results=settings['save_gs'])
        self.s1_status['GS'] = True
        
        # data quality step
        results = calwebb_detector1.dq_init_step.DQInitStep.call(results,
                                                                 output_dir=self.s1_red_dir,
                                                                 save_results=settings['save_dq'])
        self.s1_status['DQ'] = True
        
        # saturation step
        results = calwebb_detector1.saturation_step.SaturationStep.call(results,
                                                                        output_dir=self.s1_red_dir,
                                                                        save_results=settings['save_sat'])
        self.s1_status['SAT'] = True
    
    
        # superbias subtraction step
        results = calwebb_detector1.superbias_step.SuperBiasStep.call(results,
                                                                      output_dir=self.s1_red_dir,
                                                                      save_results=settings['save_sb'])
        self.s1_status['SB'] = True
        
        # reference pixel correction step
        results = calwebb_detector1.refpix_step.RefPixStep.call(results,
                                                                output_dir=self.s1_red_dir,
                                                                save_results=settings['save_rp'])
        self.s1_status['RP'] = True
        
        # linearity step
        results = calwebb_detector1.linearity_step.LinearityStep.call(results,
                                                                      output_dir=self.s1_red_dir,
                                                                      save_results=settings['save_lin'])
        self.s1_status['LIN'] = True
        
        # dark current step
        if not settings['skip_dc']:
            results = calwebb_detector1.dark_current_step.DarkCurrentStep.call(results,
                                                                               output_dir=self.s1_red_dir,
                                                                               save_results=settings['save_dc'])
            self.s1_status['DC'] = True
       
        # jump detection step
        if not settings['skip_jd']:
            results = calwebb_detector1.jump_step.JumpStep.call(results,
                                                                output_dir=self.s1_red_dir,
                                                                save_results=settings['save_jd'])
            self.s1_status['JD'] = True
        
        # ramp fitting step
        results = calwebb_detector1.ramp_fit_step.RampFitStep.call(results,
                                                                   output_dir=self.s1_red_dir,
                                                                   save_results=settings['save_rf'])
        self.s1_status['RF'] = True
        
        # gain scaling step
        if not settings['skip_gain']:
            results = calwebb_detector1.gain_scale.RampFitStep.call(results,
                                                                    output_dir=self.s1_red_dir,
                                                                    save_results=settings['save_gain'])
            self.s1_status['GAIN'] = True


# ====== Stage 2 steps ======        
    
def run_stage2(self,settings):
    
    for seg in self.segs:
        #TODO: add option to start from gain scale file instead of ramp fit
        rampfit_file = self.s1_red_dir+'/'+self.obs_handle+'-seg'+seg+'_nis_1_rampfitstep.fits'
        
        # world coordinate system step
        result = calwebb_spec2.assign_wcs_step.AssignWcsStep.call(rampfit_file,
                                                                  output_dir=self.s2_red_dir,
                                                                  save_results=settings['save_wcs'])
        self.s2_status['WCS'] = True
        
        # flatfield step
        result = calwebb_spec2.flat_field_step.FlatFieldStep.call(result, 
                                                                  output_dir=self.s2_red_dir,
                                                                  save_results=settings['save_ff'])
        self.s2_status['FF'] = True
        
        # source type step
        result = calwebb_spec2.srctype_step.SourceTypeStep.call(result,
                                                                output_dir=self.s2_red_dir,
                                                                save_results=settings['save_st'])
        self.s2_status['ST'] = True
