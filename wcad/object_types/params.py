#!/usr/bin/env python2
#
# Copyright 2015 by Ss Cyril and Methodius University in Skopje, Macedonia
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Branislav Gerazov, October 2015
#   Aleksandar Gjoreski, October 2015
#   Pierre-Edouard Honnet, October 2015
#

import numpy as np

class Params:
    def __init__(self):

        ####################
        # General parameters
        ####################

        self.sampling_rate = None  # defined in wave_input
        self.frame_shift = 5  # in milliseconds, so frameRate is 200Hz
        self.frame_rate = 1000 / self.frame_shift
        self.sample_shift = None  # defined in wave_input

        self.frame_size = 50
        self.frame_size_samples = None  # defined in wave_input

        self.f0min = 80
        self.f0max = 250

        # Results directory : 'keep_same' (default), 'delete_contents', 'make_separate'
        self.overwrite_results_dir = 'keep_same'
        # use system's libraries or the ones in the local bin/lib/ folder
        self.use_system_libs = True

        ###################
        # Weighing function
        ###################

        # Type of Weight: 'pov' (default), 'energy', 'pov_energy_mult','pov_energy_add', 'binary_both'
        # *all of them can be made binary with option self.weight_binary
        self.weight_type = 'pov_energy_mult'
        self.weight_alpha = 0.5  # for add: weight = alpha*pov + (1-alpha)*energy
        # normalizations that do nothing for the results : )
        self.weight_norm = True  # normalization of the weight func
        self.weight_norm_mult = False  # normalization of the pov and energy before multiplying in pov_energy_mult

        # binary weight
        self.weight_binary = False   # makes weight 0 or 1 depending on threshold:
        self.weight_th = 0.4  # threshold for determining binary value of the weight func (0.001 for real good atoms)
        self.weight_med_filt = 5  # implement a median filt to the binary weight

        # Silence everything below thresh - deprecated
        # sil_th = max(energy) * sil_th_coeff
        self.sil_th_coeff = 0

        ########################
        # Phrase atom extraction
        ########################

        self.phrase_type = 'fix_pos'  # only this supported

        # allow multiple phrase components - experimental
        self.multiphrase = False
        # frequency of the lowpass filter
        self.mulphras_high_freq = 0.5  # Hz
        # when looking for pauses in speech
        self.mulphras_silence_tresh = 0.1  # of the Amp, +/-
        self.mulphras_min_len = 0.075  # s
        # maximum distance when looking for pause to align with phrase curve lp minimum
        self.mulphras_silence_max_dist = 0.400  # seconds
        # count as long pause (phrase boundary) if is longer than
        self.mulphras_longpause_min_len = 0.200  # seconds
        # minimum lenght of phrase component (discarding shorter)
        self.mulphras_phrasecomp_min_len = 1  # seconds
        # padding to every (except last) phrase before sending to phrase extractor
        self.mulphras_post_pad_sec = 1  # seconds
        # (temporary) if debugging plot/show more info
        self.mulphras_debugging = False

        # reference function for fixing position: pov, energy, weight,
        # first (first to go over thresh), both (must be above thresh)
        # can also be annotations
        self.fix_pos_ref_func = 'annotations'

        # else:
        self.fix_pos_ref_func_th = 0.01  # threshold for seeking peak (ref_func is normalized)
        self.fix_pos_ref_func_t = 0  # s after threshold for seeking peak to start phrase - if 0 then threshold location

        # find the last point in the ref func that is > th and don't take it into account for the phrase
        self.fix_pos_ref_func_end = True
        self.fix_pos_ref_func_end_th = 0.001  # threshold for seeking end of utterance (smaller than at beginning)
        self.fix_pos_ref_func_end_t = 0.15   # s to discard before threshold from end
        self.fix_pos_ref_func_end_t_ext = 0  # s to extend for atom decomposition
        self.fix_pos_ref_func_end_s_ext = self.fix_pos_ref_func_end_t_ext * self.frame_rate   # samples

        self.theta_up = 0.5  # for mutant gama atoms
        # Use half atoms for the phrase component
        self.half_phrase = False

        # Thetas for phrase component
        ph_thetas = np.arange(0.1, 1, 0.1)
        ph_thetas = np.concatenate((ph_thetas, np.arange(1, 2, 0.2)))
        # ph_thetas = np.concatenate((ph_thetas, np.arange(2, 3.5, 0.5)))
        ph_thetas = np.concatenate((ph_thetas, np.arange(2, 10.5, 0.5)))
        self.phrase_thetas_range = ph_thetas

        # LP filtering before phrase extraction
        self.use_lowpass_phrase = False
        self.lp_fg = 2  # Hz

        self.boundaries_from_annotation = None

        #######################
        # Local atom extraction
        #######################

        self.atom_ex_pad_pre = 0.1  # resize before local atom extraction
        self.atom_ex_pad_post = 0.2

        # type of extraction: 'corr', 'wcorr' (default), 'womp' - weighted Orthogonal MP
        self.local_type = 'womp'

        # Atoms' sign: 'positive_only', 'both' (default)
        self.atoms_sign = 'both'

        # Multiply wCorr with Corr. Avoids deadlock when Corr is near zero
        self.mul_wcorr_corr = True

        # Dictionary parameters
        self.local_atoms_thetas = np.arange(0.01, 0.055, 0.005)
        # k can also be a list like [4, 5, 6]
        self.k = [6]

        # Trimming - for theta 0.0005-0.05 there are 98 atoms, but 16 after trimming with 0.98 (22 w/ 0.99)
        self.trim_dictionary = False
        # Maximum cross-correlation between dictionary elements
        self.dict_max_corr = 0.99

        # Stopping conditions
        # Use wcorr as stopping condition (default: None, not using)
        self.wcorr_limit = None
        # Number of atoms to extract (default: None, decision based on length)
        self.num_atoms = None
        # Minimum residue amplitude (default: 0.3)
        self.min_atom_amp = 0.3
        self.stop_at_first_min_atom = True
        # Minimum residue amplitude
        self.min_res_amp = 0
        # Avoid possible deadlock when atom's position is equal with the previous atom position
        self.break_deadlock = False

        ######################
        # Visualization/Saving
        ######################

        # Save plots
        self.save_plot = False
        # Save plots in figures folder
        self.save_plot_figures = True

        # Show plots
        self.show_plot = True
        # plot pov colored pitch or weight colored
        self.plot_pov_color = True
        # plot Praat textgrids
        self.plot_textgrids = True
        # plot pitch_max_in and out
        self.plot_pitch_max_in_out = True

        # Model
        # Save model
        self.save_model = False
