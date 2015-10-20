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

from ConfigParser import ConfigParser
import os
from datetime import datetime
from time import time
import shutil


class Paths:
    def __init__(self, args, params):
        self.base_dir = os.path.abspath('.')

        # From Config file
        parser = ConfigParser()
        parser.read('paths.cfg')
        # self.libs_dir = parser.get('BasePaths', 'Libs')
        self.kaldi = parser.get('BasePaths', 'KaldiExtractor')

        # From command line
        self.wav = os.path.abspath(args[0])
        self.res_dir = os.path.abspath(args[1])
        wav_path, wav_tail = os.path.split(self.wav)
        wav_name, wav_ext = os.path.splitext(wav_tail)
        if params.overwrite_results_dir == 'make_separate':
            self.res_dir += '-' + wav_name
            timestamp = time()
            datetime_postfix = datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d-%H-%M-%S')
            self.res_dir += '-'+datetime_postfix
            os.makedirs(self.res_dir)
        elif os.path.isdir(self.res_dir):
            if params.overwrite_results_dir == 'delete_contents':
                shutil.rmtree(self.res_dir)
                os.makedirs(self.res_dir)
        else:   # keep_same
            if not os.path.exists(self.res_dir):
                os.makedirs(self.res_dir)

        # make separate folder for figures
        # self.figure = os.path.join(self.res_dir, wav_name)
        if params.save_plot_figures:
            self.figure = os.path.join(self.base_dir, 'figures')
            if not os.path.exists(self.figure):
                os.makedirs(self.figure)
            self.figure = os.path.join(self.figure, wav_name)
        else:
            self.figure = os.path.join(self.res_dir, wav_name)
        
        self.modelpath = os.path.join(self.res_dir, wav_name)

        self.textgrid_folder = os.path.join(wav_path, 'annotations/')
