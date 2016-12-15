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

from .object_types.atom import *


class DictionaryGenerator:
    def __init__(self, params, paths):
        self.params = params
        self.corr_tresh = params.dict_max_corr
        self.k = params.k
        self.thetas = params.local_atoms_thetas
        self.fs = params.frame_rate
        self.trim_dictionary = params.trim_dictionary

    def compute(self):
        # Generating all atoms
        dictionary = []
        for k in self.k:
            for theta in self.thetas:
                new_atom = GammaAtom(k, theta, self.fs)
                dictionary.append(new_atom)

        if self.trim_dictionary:
            # The redundant should be removed
            for atom1 in dictionary:
                for atom2 in dictionary:
                    correlation = np.correlate(atom1.curve, atom2.curve, 'full')
                    if np.max(correlation) > self.corr_tresh and atom1 is not atom2:
                        dictionary.remove(atom2)

        return dictionary
