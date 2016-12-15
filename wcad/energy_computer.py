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
#   Pierre-Edouard Honnet, October 2015
#   Phil Garner, October 2015
#

import numpy as np
import numpy.linalg as linalg


class EnergyComputer:
    def __init__(self):
        return None

    #  frame done without a for loop
    @staticmethod
    def frame(a, size=512, period=256, pad=True):
        if pad:
            # This ensures that frames are aligned in the centre
            a = np.pad(a, (size/2, size/2), 'edge')

        n_frames = (a.size - (size-period)) // period
        indf = period * np.arange(n_frames)
        indf = np.matrix(indf)
        indf = indf.transpose()
        indf = np.tile(indf, size)
        ind = indf + np.arange(size)

        frames = a[ind]  # sooo smooth : )
        #                  this actually gives an array!
        return frames

    # Frame energy
    @staticmethod
    def energy(a):
        e = linalg.norm(a, axis=1)**2  # they're in rows
        # in numpy versions <1.8 axis is not implemented
        # e = np.apply_along_axis(np.linalg.norm, 1, a) ** 2
        return e

    # Window
    # It's trivial, but helps the program look good
    @staticmethod
    def window(a, w):
        return a*w

    def compute(self, wav_input, params):
        f = self.frame(wav_input, params.frame_size_samples, params.sample_shift)
        f = self.window(f, np.hanning(params.frame_size_samples))
        e = self.energy(f)

        return e
