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
#

import io
import json


class ModelSaver():
    def __init__(self, model, params, paths):
        self.model = model
        self.params = params
        self.paths = paths

    def __repr__(self):
        return '<ModelSaver(%s)>' % self.model

    @staticmethod
    def set_default(obj):
        if isinstance(obj, set):
            return list(obj)
        raise TypeError

    def save(self):
        if not self.params.save_model:
            return

        pitch_log = self.model.original_pitch
        phrase = self.model.phrase
        atoms = self.model.atoms
        reconstruction = self.model.reconstruction
        model_path = self.paths.modelpath + '.adim'
        # no pov for now / energy
        
        with io.open(model_path, 'w', encoding='utf-8') as f:
            pitch_log = pitch_log.tolist()
            reconstruction = reconstruction.tolist()

            curve1 = {}
            curve2 = {}
            curve1['original_f0'] = pitch_log
            curve2['reconstructed_f0'] = reconstruction
            f.write(unicode(json.dumps(curve1, ensure_ascii=False)))
            f.write(unicode('\n'))
            f.write(unicode(json.dumps(curve2, ensure_ascii=False)))
            f.write(unicode('\n'))

            # For phrase component(s)
            if phrase.__class__.__name__ == 'LinearAtom':
                phrase_entry = {}
                phrase_entry['atom_type'] = 'LinearAtom'
                phrase_entry['a'] = phrase.a
                phrase_entry['b'] = phrase.b
                phrase_entry['fs'] = phrase.fs
                phrase_entry['amp'] = phrase.amp
                phrase_entry['pos'] = phrase.position
                phrase_entry['length'] = phrase.length
                f.write(unicode(json.dumps(phrase_entry, ensure_ascii=False, indent=2, sort_keys=False)))
                f.write(unicode('\n'))
            elif phrase.__class__.__name__ == 'mutantGammaAtom':
                phrase_entry = {}
                phrase_entry['atom_type'] = 'mutantGammaAtom'
                phrase_entry['k'] = phrase.k
                phrase_entry['theta'] = phrase.theta
                phrase_entry['fs'] = phrase.fs
                phrase_entry['amp'] = phrase.amp
                phrase_entry['pos'] = phrase.position
                phrase_entry['length'] = phrase.length
                f.write(unicode(json.dumps(phrase_entry, ensure_ascii=False, indent=2, sort_keys=False)))
                f.write(unicode('\n'))
            else:
                print "Unknown type of phrase atom:", phrase.__class__.__name__
                f.write(unicode("{\n  "))
                f.write(unicode("Unknown type of phrase atom:"+phrase.__class__.__name__))
                f.write(unicode("\n}"))
                f.write(unicode('\n'))

            # For all the atoms
            n_atoms = len(atoms)
            for a in range(n_atoms):
                a_type = atoms[a].__class__.__name__
                if a_type == 'GammaAtom':
                    basic_entry = {}
                    basic_entry['atom_type'] = 'GammaAtom'
                    basic_entry['atom'] = a
                    basic_entry['k'] = atoms[a].k
                    basic_entry['theta'] = atoms[a].theta
                    basic_entry['fs'] = atoms[a].fs
                    basic_entry['amp'] = atoms[a].amp
                    basic_entry['pos'] = atoms[a].position
                    basic_entry['length'] = atoms[a].length
                    f.write(unicode(json.dumps(basic_entry, ensure_ascii=False, indent=2)))
                    f.write(unicode('\n'))
                else:
                    print "Unknown type of atom:", atoms[a].__class__.__name__
