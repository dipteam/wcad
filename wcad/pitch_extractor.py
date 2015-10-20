#!/usr/bin/env python2
#
# Copyright 2015 by Ss Cyril and Methodius University in Skopje, Macedonia
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Branislav Gerazov, October 2015
#   Aleksandar Gjoreski, October 2015
#

import os
import subprocess
import tempfile
from scipy import signal
from .object_types.pitch_curve import *
import numpy as np


class PitchExtractor:
    def __init__(self, wave, params, paths):
        self.wave = wave
        self.params = params
        self.paths = paths
        # Make temp files
        temp_dir_prefix = paths.base_dir + os.sep + 'tmp-'
        self.temp_dir = tempfile.mkdtemp(prefix=temp_dir_prefix)
        self.kaldi_scp_file = self.temp_dir + os.sep + self.wave.basename + '.scp'
        self.kaldi_ark_file = self.temp_dir + os.sep + self.wave.basename + '.ark'

    def __del__(self):
        os.remove(self.kaldi_scp_file)
        os.remove(self.kaldi_ark_file)
        os.rmdir(self.temp_dir)

    @staticmethod
    def compute_weight(pov, energy, params):
        if params.weight_type != 'binary_both':
            if params.weight_type == 'pov':
                weight = pov
            elif params.weight_type == 'energy':
                weight = energy
            elif params.weight_type == 'pov_energy_mult':
                if params.weight_norm_mult is True:
                    pov = (pov - np.min(pov)) / (np.max(pov) - np.min(pov))
                    energy = (energy - np.min(energy)) / (np.max(energy) - np.min(energy))
                weight = pov * energy
            elif params.weight_type == 'pov_energy_add':
                alpha = params.weight_alpha
                weight = alpha*pov + (1-alpha)*energy

            if params.weight_norm is True:
                weight = (weight - np.min(weight)) / (np.max(weight) - np.min(weight))

            if params.weight_binary is True:  # threshold it to binary
                weight = np.where(weight > params.weight_th, 1, 0)
                if params.weight_med_filt:
                    weight = signal.medfilt(weight, params.weight_med_filt)

        else:  # binary decision based on both
            pov = (pov - np.min(pov)) / (np.max(pov) - np.min(pov))
            energy = (energy - np.min(energy)) / (np.max(energy) - np.min(energy))
            weight = np.where((pov > params.weight_th) & (energy > params.weight_th), 1, 0)
            if params.weight_med_filt:
                weight = signal.medfilt(weight, params.weight_med_filt)

        return weight

    def compute(self):
        kaldi_scp_file = self.kaldi_scp_file
        kaldi_ark_file = self.kaldi_ark_file

        # Write scp file for Kaldi
        text_file = open(kaldi_scp_file, "w")
        text_file.write(self.wave.filename + ' ' + self.wave.filepath + '\n')
        text_file.close()

        # Set Kaldi's parameters:
        kaldi_params = ' --sample-frequency=' + str(self.wave.fs)
        kaldi_params += ' --frame-length=' + str(self.params.frame_size)
        kaldi_params += ' --frame-shift=' + str(self.params.frame_shift)
        kaldi_params += ' --min-f0=' + str(self.params.f0min)
        kaldi_params += ' --max-f0=' + str(self.params.f0max)
        kaldi_params += ' scp:' + self.kaldi_scp_file
        kaldi_params += ' ark:' + self.kaldi_ark_file

        kaldi_libs = '' if self.params.use_system_libs else 'export LD_LIBRARY_PATH=' + self.paths.libs_dir + '; '

        subprocess.call(kaldi_libs + self.paths.kaldi + kaldi_params, shell=True)
        # agjoreski to-do: pipe the output directly to phthon, w/o using output file

        # Kaldi's result is stored as tuple (nccf, f0), with file's name in the header
        file_obj = open(kaldi_ark_file, 'rb')
        file_obj.seek(len(self.wave.filename) + 16)                             # skipping the header
        kaldi_pitch_t = np.dtype([('nccf', np.float32), ('f0', np.float32)])    # new type: (nccf, f0)
        kaldi_data = np.fromfile(file_obj, dtype=kaldi_pitch_t)
        file_obj.close()

        f0 = kaldi_data['f0']
        nccf = kaldi_data['nccf']

        # Convert NCCF to POV. According to Kaldi paper:
        # [ l = log(p(voiced)/p(unvoiced)) ]
        a = np.abs(nccf)
        l = -5.2 + 5.4*np.exp(7.5*(a - 1)) + 4.8*a - 2*np.exp(-10*a) + 4.2*np.exp(20*(a - 1)) 
        pov = 1./(1 + np.exp(-l))

        # Then we need to extend both sides because Kaldi starts from half window length
        energy = self.wave.energy
        len_energy = len(energy)
        len_f0 = len(f0)
        dif = len_energy - len_f0
        dif_h = dif // 2
        f0 = np.pad(f0, (dif_h, dif_h), 'edge')
        pov = np.pad(pov, (dif_h, dif_h), 'edge')

        if dif % 2:  # odd, then 1 is missing
            f0 = np.pad(f0, (1, 0), 'edge')
            pov = np.pad(pov, (1, 0), 'edge')

        # Calculate weight
        weight = self.compute_weight(pov, energy, self.params)

        f0_log = np.log(f0)

        # Temporary
        fs = self.params.frame_rate
        lp_fg = self.params.lp_fg  # default 0.5 Hz from Mixdorff
        lp_wg = lp_fg / (0.5 * fs)
        b, a = signal.butter(4, lp_wg, btype='lowpass')
        f0_log_filt = signal.filtfilt(b, a, f0_log, padtype='odd', padlen=len(f0_log)-1)

        # plt.plot(f0_log)
        # plt.plot(f0_log_filt)
        # plt.show()

        # phrase = Atom(f0_log_filtered, self.params.frame_rate)
        # return phrase
        # end if lowpass_phrase

        pitch = PitchCurve(f0, f0_log, f0_log_filt, pov, weight)

        return pitch




