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

from wcad import *
from optparse import OptionParser
import time
import shutil

(opts, args) = OptionParser().parse_args()
# args = ['audio/A01.wav', 'results']  # testing example
# args = ['audio/A01E.wav', 'results']  # testing example
# args = ['audio/A02.wav', 'results']  # testing example
# args = ['audio/A02E.wav', 'results']  # testing example
# args = ['audio/A03.wav', 'results']  # testing example
# args = ['audio/A03E.wav', 'results']  # testing example
# args = ['audio/A04.wav', 'results']  # testing example
# args = ['audio/A04E.wav', 'results']  # testing example
# args = ['audio/A05.wav', 'results']  # testing example
# args = ['audio/A05E.wav', 'results']  # testing example

params = Params()
paths = Paths(args, params)

start_t = time.time()
wave = WaveInput(paths.wav, params).read()
pitch = PitchExtractor(wave, params, paths).compute()

phrase = MultiphraseExtractor(pitch, wave, params, paths).compute()

dictionary = DictionaryGenerator(params, paths).compute()
atoms = AtomExtrator(wave, pitch, phrase, dictionary, params, paths).compute()

model = ModelCreator(phrase, atoms, pitch).compute()
print 'Model created in %s seconds' % (time.time() - start_t)

ModelSaver(model, params, paths).save()

ModelPlotterc(wave, model, pitch, params, paths).plot()

# clean up
if params.overwrite_results_dir == 'none':
    shutil.rmtree(paths.res_dir)

