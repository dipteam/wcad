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
import numpy.linalg as linalg
import math


class Atom:
    """
    A superclass for all types of atoms. The method get_curve will be common
    to all classes of atoms and will have to be implemented separately.
    So it will be possible to make operations with the atom without knowing its class
    """
    def __init__(self, curve, fs=None, amp=1, position=0,
                 pitch_max_in=None, pitch_max_out=None, phrase_max_out=None):
        self.curve = curve
        self.fs = fs
        self.amp = amp
        self.position = position
        self.pitch_max_in = pitch_max_in
        self.pitch_max_out = pitch_max_out
        self.phrase_max_out = phrase_max_out

    def get_curve(self):
        if self.curve is not None:
            return self.curve
        else:
            raise Exception("Function not generated")

    def generate_curve(self):
        if self.curve is None:
            self.curve = self.get_curve()

    def regenerate_curve(self):
        self.curve = None
        self.curve = self.get_curve()

    def get_padded_curve(self, wanted_len):
        self.generate_curve()

        position = self.position
        curve = self.curve
        curve_len = len(self.curve)

        # Crop if position is negative
        if position < 0:
            crop_len = -position
            curve = curve[crop_len:]
            position = 0
            curve_len = len(curve)

        pre_pad = position
        if position + curve_len <= wanted_len:
            post_pad = wanted_len - (position + curve_len)
        elif pre_pad > wanted_len:
            print 'WARNING: position %d > wanted_len %d while padding the atom.' % (position, wanted_len)
            return np.zeros(wanted_len)
        else:
            # Crop the end
            post_pad = 0
            crop_len = position + curve_len - wanted_len
            curve = curve[:-crop_len]
        padded_curve = np.pad(curve, (pre_pad, post_pad), 'constant')

        # # Apply the gain:
        padded_curve *= self.amp

        return padded_curve

    def get_peak_position(self, wanted_len):
        """ May be rewritten not to use the get_padded_curve() method."""
        curve = self.get_padded_curve(wanted_len)
        return np.argmax(curve)


class GammaAtom(Atom):
    def __init__(self, k, theta, fs, amp=1, position=0, length=None):
        self.curve = None
        self.k = k
        self.theta = theta
        self.fs = fs
        self.amp = amp
        self.position = position
        self.length = length
        self.curve = self.get_curve()

    def get_curve(self):
        # If is already computed just return it
        if self.curve is not None:
            return self.curve

        length = 20     # maximum length in sec - see later if this need to be calculated differently
        k = self.k
        theta = self.theta

        # This is our time vector (just the length of the gamma atom):
        t = np.linspace(0, length, length*self.fs, endpoint=False)

        # np.vectorize is not really vectorized, it's just nicer way to loop
        gamma_function = np.vectorize(lambda tt:
                                      1/(math.gamma(k)*theta**k)*tt**(k-1)*math.exp(-tt/theta))
        gamma_atom = gamma_function(t)

        # # Now shorten the atoms
        thresh = 1e-5  # don't go above 1e-5 because atoms will be shorter than f0 (for Macedonian sentence at least :)
        gamma_atom_th_ind = np.where(gamma_atom > thresh)  # indexes of elements above thresh
        gamma_atom = gamma_atom[gamma_atom_th_ind]
        gamma_atom -= np.min(gamma_atom)

        gamma_atom /= linalg.norm(gamma_atom)  # norm-2 of 1

        self.length = len(gamma_atom)

        return gamma_atom


class mutantGammaAtom(Atom):
    def __init__(self, k, theta, theta_up, fs, amp=1, position=0, length=None):
        self.curve = None
        self.k = k
        self.theta = theta
        self.fs = fs
        self.amp = amp
        self.position = position
        self.length = length
        self.theta_up = theta_up
        self.curve = self.get_curve()

    def get_curve(self):
        # If is already computed just return it
        if self.curve is not None:
            return self.curve

        length = 20     # maximum length in sec - see later if this need to be calculated differently
        k = self.k
        theta = self.theta
        theta_up = self.theta_up

        # This is our time vector (just the length of the gamma atom):
        t = np.linspace(0, length, length*self.fs, endpoint=False)

        # np.vectorize is not really vectorized, it's just nicer way to loop
        gamma_function_up = np.vectorize(lambda tt: 1/(math.gamma(k)*theta_up**k)*tt**(k-1)*math.exp(-tt/theta_up))
        gamma_function_down = np.vectorize(lambda tt: 1/(math.gamma(k)*theta**k)*tt**(k-1)*math.exp(-tt/theta))
        gamma_atom_up = gamma_function_up(t)
        gamma_atom_down = gamma_function_down(t)

        # stick them together : )
        gamma_atom_up = gamma_atom_up[:np.argmax(gamma_atom_up)] / np.max(gamma_atom_up)
        gamma_atom_down = gamma_atom_down[np.argmax(gamma_atom_down):] / np.max(gamma_atom_down)
        gamma_atom = np.concatenate((gamma_atom_up, gamma_atom_down))

        gamma_atom /= linalg.norm(gamma_atom)  # this preserves array and eliminates for

        return gamma_atom
