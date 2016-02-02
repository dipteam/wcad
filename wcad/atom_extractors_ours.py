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

import time
from .signal_comparer import *


class AtomExtrator:
    def __init__(self, wave, pitch, phrase, dictionary, params, paths):
        self.wave = wave
        self.pitch = pitch
        self.phrase = phrase
        self.dictionary = dictionary
        self.params = params
        self.paths = paths

    def compute(self):
        the_chosen_one = self.params.local_type
        extractors = {
            'corr': OurAtomExtractor,
            'wcorr': OurWeightedAtomExtractor,
            'womp': OMPWeightedAtomExtractor
        }
        return extractors[the_chosen_one](self.wave, self.pitch, self.phrase,
                                          self.dictionary, self.params, self.paths).compute()


class OurAtomExtractor():
    """
    MPTK-like extractor. Uses not normalized cross-correlation.
    """

    def __init__(self, wave, pitch, phrase, dictionary, params, paths):
        self.wave = wave
        self.pitch = pitch
        self.phrase = phrase
        self.dictionary = dictionary
        self.params = params
        self.paths = paths  #

    def extrapolate(self):
        pre_pad = self.params.atom_ex_pad_pre * self.params.frame_rate
        post_pad = self.params.atom_ex_pad_post * self.params.frame_rate

        self.phrase.curve = np.pad(self.phrase.curve, (pre_pad, post_pad), 'minimum')
        self.pitch.f0 = np.pad(self.pitch.f0, (pre_pad, post_pad), 'minimum')
        self.pitch.f0_log = np.pad(self.pitch.f0_log, (pre_pad, post_pad), 'minimum')
        self.pitch.f0_log_filt = np.pad(self.pitch.f0_log_filt, (pre_pad, post_pad), 'minimum')
        self.pitch.pov = np.pad(self.pitch.pov, (pre_pad, post_pad), 'constant')  # zeros
        self.pitch.weight = np.pad(self.pitch.weight, (pre_pad, post_pad), 'constant')
        self.wave.energy = np.pad(self.wave.energy, (pre_pad, post_pad), 'constant')

    def deextrapolate(self, atoms):
        pre_pad = self.params.atom_ex_pad_pre * self.params.frame_rate
        post_pad = self.params.atom_ex_pad_post * self.params.frame_rate

        self.phrase.curve = self.phrase.curve[pre_pad:-post_pad]
        # self.phrase.length -= pre_pad + post_pad
        self.pitch.f0 = self.pitch.f0[pre_pad:-post_pad]
        self.pitch.f0_log = self.pitch.f0_log[pre_pad:-post_pad]
        self.pitch.f0_log_filt = self.pitch.f0_log_filt[pre_pad:-post_pad]
        self.pitch.pov = self.pitch.pov[pre_pad:-post_pad]
        self.pitch.weight = self.pitch.weight[pre_pad:-post_pad]
        self.wave.energy = self.wave.energy[pre_pad:-post_pad]

        # to-do: compute this while assigning position
        for atom in atoms:
            atom.position = atom.position - pre_pad

    def compute(self):
        self.extrapolate()
        # We want to fit the atoms to f0_diff
        f0_diff = self.pitch.f0_log - self.phrase.curve
        f0_diff_len = len(f0_diff)

        # Stopping parameters
        n_atoms = f0_diff_len//20 if self.params.num_atoms is None else self.params.num_atoms

        atoms = []
        print 'Extracting atoms in max ', n_atoms, ' iterations '
        start_time = time.time()

        for iteration in range(n_atoms):
            max_corr_atom = None
            max_corr = None

            for atom in self.dictionary:
                correlation = np.correlate(atom.curve, f0_diff, 'full')
                corr_ind = np.argmax(np.abs(correlation))  # index of maximum correlation
                atom_corr = correlation[corr_ind]

                if max_corr is None or np.abs(atom_corr) > np.abs(max_corr):
                    max_corr = atom_corr
                    max_corr_atom = atom
                    max_corr_atom.amp = np.correlate(atom.curve, f0_diff, 'full')[corr_ind]
                    position_max_corr = f0_diff_len+len(atom.curve)-1 - corr_ind - len(atom.curve)//2
                    max_corr_atom.position = position_max_corr - len(atom.curve)//2
                # end for atom in dictionary

            # Best atom found. Save it
            max_corr_atom.regenerate_curve()
            new_atom = copy.deepcopy(max_corr_atom)    # w/o this python uses the same reference later

            new_atom_padded = new_atom.get_padded_curve(f0_diff_len)
            f0_diff -= new_atom_padded

            if self.params.min_atom_amp and abs(max_corr_atom.amp) < self.params.min_atom_amp:
                print 'Atom amplitude is %f bellow %f. Not appending this one.' \
                      % (max_corr_atom.amp, self.params.min_atom_amp)
                continue

            atoms.append(new_atom)

            # # Substraction plotting
            # plt.subplot(211)
            # plt.plot(new_atom_padded, 'm', linewidth=2.0)
            # plt.plot(f0_diff+new_atom_padded, 'k')
            #
            # plt.subplot(212)
            # plt.plot(f0_diff+new_atom_padded, 'c--')
            # plt.plot(f0_diff, 'k')
            # plt.show()

            print 'iteration: ', iteration, 'atom.theta: ', max_corr_atom.theta, 'pos: ', \
                max_corr_atom.position, 'amp: ', max_corr_atom.amp

            # Check if stopping condition is reached
            if np.max(np.fabs(f0_diff)) < self.params.min_res_amp:
                print 'Residue amplitude bellow %f. Stopping.' % self.params.min_res_amp
                break

            if self.params.wcorr_limit is not None:
                # compute reconstruction
                recon_len = len(self.pitch.f0_log)
                reconstruction = np.zeros(recon_len)
                reconstruction += self.phrase.curve
                for atom in atoms:
                    padded_curve = atom.get_padded_curve(recon_len)
                    reconstruction += padded_curve
                # compute wcorr
                wcorr = SignalComparer(reconstruction, self.pitch.f0_log, self.pitch.weight).simple_norm_wcorr()
                if wcorr >= self.params.wcorr_limit:
                    print 'Reached wcorr greater than %f. Stopping.' % self.params.wcorr_limit
                    break
            # end for i in range(n_atoms)

        print "Extracted %d atoms in %s seconds" % (len(atoms), time.time() - start_time)

        self.deextrapolate(atoms)
        return atoms


class OurWeightedAtomExtractor(OurAtomExtractor):
    """
    Uses weighted cross-correlation.
    """

    def compute(self):
        # self.extrapolate()
        # pre_pad = self.params.atom_ex_pad_pre * self.params.frame_rate

        # We want to fit the atoms to weighted f0_diff
        f0 = self.pitch.f0_log
        f0_diff = self.pitch.f0_log - self.phrase.curve
        f0_diff_len = len(f0_diff)
        weight = self.pitch.weight
        if self.params.phrase_type == 'fix_pos':
            if hasattr(self.phrase, 'atoms'):
                start = self.phrase.atoms[0].pitch_max_in
                end = self.phrase.atoms[-1].pitch_max_out
            else:
                start = self.phrase.pitch_max_in
                end = self.phrase.pitch_max_out

            f0_diff[:start] = 0
            weight[:start] = 0
            f0_diff[end + self.params.fix_pos_ref_func_end_s_ext:] = 0
            weight[end + self.params.fix_pos_ref_func_end_s_ext:] = 0

        # Stopping parameters
        n_atoms = f0_diff_len//20 if self.params.num_atoms is None else self.params.num_atoms

        atoms = []
        print 'Extracting atoms in max ', n_atoms, 'iterations.'
        start_time = time.time()

        for iteration in range(n_atoms):
            max_wcorr = None
            max_wcorr_atom = None
            max_wcorr_index = None
            f0_diff_w = weight * f0_diff

            for atom in self.dictionary:
                top_corr = np.correlate(f0_diff_w, atom.curve, 'full')
                down1 = np.sqrt(np.sum(weight * (f0_diff**2)))
                down2 = np.sqrt(np.correlate(weight, atom.curve**2, 'full'))
                weighted_corr = np.where((down2 * down1) != 0, top_corr / (down1*down2), 0)
                corr = np.correlate(f0_diff, atom.curve, 'full')

                if self.params.mul_wcorr_corr:
                    weighted_corr = weighted_corr * corr

                wcorr_ind = np.argmax(np.abs(weighted_corr))  # index of maximum correlation
                atom_wcorr = weighted_corr[wcorr_ind]

                if max_wcorr is None or np.abs(atom_wcorr) > np.abs(max_wcorr):
                    max_wcorr = atom_wcorr
                    max_wcorr_index = wcorr_ind
                    max_wcorr_atom = atom
                    max_wcorr_atom.position = wcorr_ind - atom.curve.size + 1
                # end for atom in dictionary

            # Best atom found
            # Compute its amplitude
            correlation = np.correlate(f0_diff, max_wcorr_atom.curve, 'full')
            max_wcorr_atom.amp = correlation[max_wcorr_index]
            max_wcorr_atom.regenerate_curve()
            new_atom = copy.deepcopy(max_wcorr_atom)    # w/o this python uses the same reference later

            # Substract it from the signal
            new_atom_padded = new_atom.get_padded_curve(f0_diff_len)
            f0_diff -= new_atom_padded

            if self.params.min_atom_amp and abs(max_wcorr_atom.amp) < self.params.min_atom_amp:
                print 'Atom amplitude is %f bellow %f. Not appending this one.' \
                          % (max_wcorr_atom.amp, self.params.min_atom_amp)
                if self.params.stop_at_first_min_atom is True:
                    break
                else:
                    continue

            # Check possible deadlock
            if self.params.break_deadlock and len(atoms) > 0 and atoms[-1].position == new_atom.position:
                recon_len = len(self.pitch.f0_log)
                reconstruction = np.zeros(recon_len)
                reconstruction += self.phrase.curve
                for atom in atoms:
                    padded_curve = atom.get_padded_curve(recon_len)
                    reconstruction += padded_curve
                # compute wcorr
                wcorr = SignalComparer(reconstruction, self.pitch.f0_log, self.pitch.weight).simple_norm_wcorr()
                print 'wCorr before stop: %f.' % wcorr
                print 'Possible deadlock. Stopping.'
                break

            atoms.append(new_atom)

            print 'iteration: ', iteration, 'atom.theta: ', max_wcorr_atom.theta, 'pos: ', \
                max_wcorr_atom.position, 'amp: ', max_wcorr_atom.amp, 'wcorr_ind:', wcorr_ind

            # Check if stopping condition is reached
            if self.params.min_res_amp and np.max(np.fabs(f0_diff_w)) < self.params.min_res_amp:
                print 'Residue amplitude bellow %f. Stopping.' % self.params.min_res_amp
                break

            if self.params.wcorr_limit is not None:
                # compute reconstruction
                recon_len = len(self.pitch.f0_log)
                reconstruction = np.zeros(recon_len)
                reconstruction += self.phrase.curve
                for atom in atoms:
                    padded_curve = atom.get_padded_curve(recon_len)
                    reconstruction += padded_curve
                # compute wcorr
                wcorr = SignalComparer(reconstruction, self.pitch.f0_log, self.pitch.weight).simple_norm_wcorr()
                if wcorr >= self.params.wcorr_limit:
                    print 'Reached wcorr of %f, greater than %f. Stopping.' \
                          % (wcorr, self.params.wcorr_limit)
                    break

            # end for i in range(n_atoms)

        print "Extracted %d atoms in %s seconds" % (len(atoms), time.time() - start_time)

        # self.deextrapolate(atoms)
        return atoms


class OMPWeightedAtomExtractor(OurAtomExtractor):
    """
    Based on OurWeightedAtomExtractor. Includes OMP
    """

    def compute(self):
        # self.extrapolate()
        pre_pad = 0  # self.params.atom_ex_pad_pre * self.params.frame_rate
        # post_pad = self.params.atom_ex_pad_post * self.params.frame_rate
        # We want to fit the atoms to weighted f0_diff
        f0 = self.pitch.f0_log
        f0_diff = self.pitch.f0_log - self.phrase.curve
        f0_diff_len = len(f0_diff)
        weight = self.pitch.weight
        recon_len = len(self.pitch.f0_log)
        reconstruction = np.zeros(recon_len)
        reconstruction += self.phrase.curve

        wcorr = SignalComparer(f0, self.phrase.curve, weight).simple_norm_wcorr()
        # wcorr = SignalComparer(f0, self.phrase.curve, weight).rmse()
        # save WCORRs for performance assesment
        wcorrs = []
        wcorrs.append(wcorr)

        if self.params.phrase_type == 'fix_pos':
            if hasattr(self.phrase, 'atoms'):
                start = self.phrase.atoms[0].pitch_max_in
                end = self.phrase.atoms[-1].pitch_max_out
            else:
                start = self.phrase.pitch_max_in
                end = self.phrase.pitch_max_out

            f0_diff[:start+pre_pad] = 0
            weight[:start+pre_pad] = 0
            f0_diff[end + self.params.fix_pos_ref_func_end_s_ext+pre_pad:] = 0
            weight[end + self.params.fix_pos_ref_func_end_s_ext+pre_pad:] = 0

        # Stopping parameters
        n_atoms = f0_diff_len//20 if self.params.num_atoms is None else self.params.num_atoms

        atoms = []
        print 'Extracting atoms in max ', n_atoms, 'iterations.'
        start_time = time.time()

        f0_nophrase = np.copy(f0_diff)

        for iteration in range(n_atoms):
            max_wcorr = None
            max_wcorr_atom = None
            max_wcorr_index = None
            f0_diff_w = weight * f0_diff

            for atom in self.dictionary:
                top_corr = np.correlate(f0_diff_w, atom.curve, 'full')
                down1 = np.sqrt(np.sum(weight * (f0_diff**2)))
                down2 = np.sqrt(np.correlate(weight, atom.curve**2, 'full'))
                weighted_corr = np.where((down2 * down1) != 0, top_corr / (down1*down2), 0)
                corr = np.correlate(f0_diff, atom.curve, 'full')

                if self.params.mul_wcorr_corr:
                    weighted_corr = weighted_corr * corr

                wcorr_ind = np.argmax(np.abs(weighted_corr))  # index of maximum correlation
                atom_wcorr = weighted_corr[wcorr_ind]

                if max_wcorr is None or np.abs(atom_wcorr) > np.abs(max_wcorr):
                    max_wcorr = atom_wcorr
                    max_wcorr_index = wcorr_ind
                    max_wcorr_atom = atom
                    # position_max_corr = full_corr_len - wcorr_ind - len(atom.curve)//2
                    max_wcorr_atom.position = wcorr_ind - atom.curve.size + 1
                # end for atom in dictionary

            # Best atom found

            # Check possible deadlock
            if self.params.break_deadlock and len(atoms) > 0 and atoms[-1].position == max_wcorr_atom.position:
                print 'Possible deadlock. Stopping.'
                print 'max_wcorr was', max_wcorr, '| max_wcorr_index was', max_wcorr_index, \
                    '| max_wcorr_atom.position was', max_wcorr_atom.position
                break

            # Compute old amplitude
            correlation = np.correlate(f0_diff, max_wcorr_atom.curve, 'full')
            max_wcorr_atom.amp = correlation[max_wcorr_index]

            max_wcorr_atom.regenerate_curve()
            new_atom = copy.deepcopy(max_wcorr_atom)    # w/o this python uses the same reference later
            atoms.append(new_atom)

            # Compute all new amplitudes
            recons_len = len(self.pitch.f0_log)
            atoms_padded = [atom.get_padded_curve(recons_len, include_amplitude=False) for atom in atoms]
            atoms_padded = np.array(atoms_padded)

            new_amplitudes, _, _, _ = np.linalg.lstsq(atoms_padded.T, f0_nophrase)

            for atom_idx, atom in enumerate(atoms):
                atom.amp = new_amplitudes[atom_idx]
                atom.regenerate_curve()

            # Substract all atoms from the signal
            atoms_padded = [atom.get_padded_curve(recons_len, include_amplitude=False) for atom in atoms]
            atoms_padded = np.array(atoms_padded)
            f0_diff = f0_nophrase - np.dot(atoms_padded.T, new_amplitudes)

            # # redo reconstruction
            # reconstruction = self.phrase.curve + np.dot(atoms_padded.T, new_amplitudes)

            # # compute wcorr
            # wcorr = SignalComparer(reconstruction, f0, weight).simple_norm_wcorr()
            # # wcorr = SignalComparer(reconstruction, f0, weight).rmse()
            # wcorrs.append(wcorr)

            if self.params.min_atom_amp and abs(max_wcorr_atom.amp) < self.params.min_atom_amp:
                print 'Atom amplitude is %f bellow %f. Not appending this one.' \
                          % (max_wcorr_atom.amp, self.params.min_atom_amp)
                if self.params.stop_at_first_min_atom is True:
                    break
                else:
                    continue

            print 'iteration: ', iteration, 'atom.theta: ', max_wcorr_atom.theta, 'pos: ', \
                max_wcorr_atom.position, 'amp: ', max_wcorr_atom.amp, 'wcorr_ind:', wcorr_ind

            # Check if stopping condition is reached
            if self.params.min_res_amp and np.max(np.fabs(f0_diff_w)) < self.params.min_res_amp:
                print 'Residue amplitude bellow %f. Stopping.' % self.params.min_res_amp
                break

            if self.params.wcorr_limit is not None:
                # compute reconstruction
                recon_len = len(self.pitch.f0_log)
                reconstruction = np.zeros(recon_len)
                reconstruction += self.phrase.curve
                for atom in atoms:
                    padded_curve = atom.get_padded_curve(recon_len)
                    reconstruction += padded_curve
                # compute wcorr
                wcorr = SignalComparer(reconstruction, self.pitch.f0_log, self.pitch.weight).simple_norm_wcorr()
                if wcorr >= self.params.wcorr_limit:
                    print 'Reached wcorr of %f, greater than %f. Stopping.' \
                          % (wcorr, self.params.wcorr_limit)
                    break

            # end for i in range(n_atoms)

        print "Extracted %d atoms in %s seconds" % (len(atoms), time.time() - start_time)
        # self.deextrapolate(atoms)
        return atoms
