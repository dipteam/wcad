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

class Wave:
    def __init__(self, fs, wave, filepath, energy, filename=None, basename=None):
        self.fs = fs
        self.wave = wave
        self.energy = energy
        self.filepath = filepath

        self.filename = filename if filename is not None else os.path.basename(self.filepath)
        self.basename = basename if basename is not None else os.path.splitext(self.filename)[0]
