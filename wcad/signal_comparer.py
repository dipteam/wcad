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
import copy


class SignalComparer():
    def __init__(self, sig1, sig2, *args):
        if len(args) == 1:
            w1 = None
            w2 = None
            w12 = args[0]
        elif len(args) == 2:
            w1 = args[0]
            w2 = args[1]
            w12 = None
        else:
            w1 = None
            w2 = None
            w12 = None
        if w1 is not None:
            assert(w1.shape[0] == w2.shape[0])
            assert(sig1.shape[0] == w1.shape[0])
        if w12 is not None:
            assert(sig1.shape[0] == w12.shape[0])
        self.sig1 = sig1
        self.sig2 = sig2
        self.w1 = w1
        self.w2 = w2
        self.w12 = w12

    def wrmse(self):
        """
        This calculates weighted RMSE using sig1 (or sig2) as reference
        given w1 and w2 which are the weight functions
        """
        sig1 = self.sig1
        sig2 = self.sig2
        w1 = self.w1
        w2 = self.w2
        w12 = self.w12
        if w12 is None:
            if np.any(w1) and np.any(w2):  # check for division w/ 0
                top = np.sum(w1 * w2 * (sig1 - sig2)**2)
                down = np.sum(w1 * w2)
                wrmse = np.sqrt(top / down)
            else:
                wrmse = 100  # random big number
        else:
            top = np.sum(w12 * (sig1 - sig2)**2)
            down = np.sum(w12)
            wrmse = np.sqrt(top / down)
        return wrmse

    def wcorr(self):
        """
        This calculates weighted correlation using sig1 (or sig2) as reference
        given w1 and w2 which are the weight functions
        """
        sig1 = self.sig1  # f0
        w1 = self.w1  # pov
        w2 = self.w2  # energy
        w12 = self.w12
        if w12 is None:
            w12 = w1 * w2

        sig2 = self.sig2  # atoms

        if sig2.ndim == 1:  # an array is passed
            top = np.sum(w12 * sig1 * sig2)
            down1 = np.sum(w12 * np.square(sig1))
            down2 = np.sum(w12 * np.square(sig2))
            if np.any(down1) and np.any(down2):  # check for division w/ 0
                wcorr = top / np.sqrt(down1*down2)
            else:
                wcorr = 0

        else:  # sig2 is a matrix
            top = np.sum(np.multiply(sig2, w12 * sig1), 1)  # gives a single column matrix
            down1 = np.sum(w12 * np.square(sig1))  # single number
            down2 = np.sum(np.multiply(w12, np.square(sig2)), 1)  # also
            with np.errstate(invalid='ignore'):
                wcorr = top / np.sqrt(down1*down2)  # this should work as is

        return wcorr

    def rmse(self):
        sig1 = self.sig1
        sig2 = self.sig2
        n = len(sig1)
        assert (len(sig1) == len(sig1))

        mse = np.sum((sig1-sig2)**2) / n
        rmse = np.sqrt(mse)

        return rmse

    def corr(self):
        sig1 = self.sig1
        sig2 = self.sig2
        assert (len(sig1) == len(sig1))

        corr = np.sum(sig1 * sig2)

        return corr

    def simple_wcorr(self):
        sig1 = self.sig1
        sig2 = self.sig2
        w = self.w12
        assert (len(sig1) == len(sig1))
        assert (len(w) == len(sig1))

        top = np.sum(w * sig1 * sig2)
        down = np.sqrt(np.sum(w * sig1**2) * np.sum(w * sig2**2))
        wcorr = float(top)/down

        return wcorr

    def simple_norm_wcorr(self):
        sig1 = copy.copy(self.sig1)
        sig2 = copy.copy(self.sig2)
        w = self.w12
        assert (len(sig1) == len(sig1))
        assert (len(w) == len(sig1))

        sig1 -= np.mean(sig1)
        sig2 -= np.mean(sig2)

        top = np.sum(w * sig1 * sig2)
        down = np.sqrt(np.sum(w * sig1**2) * np.sum(w * sig2**2))
        wcorr = float(top)/down

        return wcorr
