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

class Model:
    def __init__(self, phrase, atoms, pitch, reconstruction):
        self.phrase = phrase
        self.atoms = atoms
        self.original_pitch = pitch
        self.reconstruction = reconstruction

        # to-do: store model error, etc
