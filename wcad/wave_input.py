#!/usr/bin/env python2
#
# Copyright 2015 by Ss Cyril and Methodius University in Skopje, Macedonia
# Copyright 2015 by Idiap Research Institute in Martigny, Switzerland 
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Branislav Gerazov, October 2015
#   Aleksandar Gjoreski, October 2015
#

import numpy as np
from scipy.io import wavfile

from wcad.energy_computer import EnergyComputer
from wcad.object_types.wave import *


class WaveInput:
    def __init__(self, wav_file, params):
        self.wav_file = wav_file
        self.params = params

    def read(self):
        fs, wave = wavfile.read(self.wav_file)

        self.params.sampling_rate = fs
        self.params.sample_shift = self.params.sampling_rate * self.params.frame_shift / 1000
        self.params.frame_size_samples = self.params.sampling_rate * self.params.frame_size / 1000

        if wave.dtype == 'int16':
            wave = np.float32(wave)
            wave /= 32768

        energy = EnergyComputer().compute(wave, self.params)

        filepath = self.wav_file
        filename = os.path.basename(filepath)        # ex 'sample.wav'
        basename = os.path.splitext(filename)[0]     # ex 'sample'

        return Wave(fs, wave, filepath, energy, filename=filename, basename=basename)
