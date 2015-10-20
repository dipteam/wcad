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

class PitchCurve:
    def __init__(self, f0, f0_log, f0_log_filt, pov, weight):
        self.f0 = f0
        self.pov = pov
        self.f0_log = f0_log
        self.f0_log_filt = f0_log_filt
        self.weight = weight