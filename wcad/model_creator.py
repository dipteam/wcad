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

import numpy as np

from .object_types.model import *


class ModelCreator():
    def __init__(self, phrase, atoms, pitch):
        self.phrase = phrase
        self.atoms = atoms
        self.pitch = pitch

    def compute(self):
        pitch_log = self.pitch.f0_log
        phrase = self.phrase
        atoms = self.atoms
        recon_len = len(self.pitch.f0_log)
        reconstruction = np.zeros(recon_len)

        if phrase is not None:
            reconstruction += phrase.curve         # padded_curve in loop when implemented

        for atom in atoms:
            padded_curve = atom.get_padded_curve(recon_len)
            reconstruction += padded_curve

        # to-do:
        # compute error, etc..

        model = Model(phrase, atoms, pitch_log, reconstruction)
        return model
