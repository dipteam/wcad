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

from __future__ import division
import time

import matplotlib.pyplot as plt
import scipy as sp
from scipy import signal

from .object_types.atom import *
from .object_types.phrase import *
from .object_types.pitch_curve import *
from .object_types.wave import *
from .signal_comparer import *
from nltk_textgrid.textgrid import *


class MultiphraseExtractor:
    def __init__(self, pitch, wave, params, paths):
        self.pitch = pitch
        self.wave = wave
        self.params = params
        self.paths = paths

    def find_silent_intervals(self, min_len, silence_tresh, medfilt=5):
        """
        Returns the interval boundaries [a, b] as a list of tuples.
        (a and b inclusive)
        """
        energy = np.copy(self.wave.energy)
        energy = np.sqrt(energy)
        energy = (energy-np.min(energy)) / (np.max(energy)-np.min(energy))
        # energy = sp.signal.medfilt(energy, medfilt)

        energy = energy

        num_samples_needed = min_len * self.params.frame_rate

        intervals = []
        in_interval = False
        samples_passed = 0
        interval_start = 0
        for i in range(len(energy)):
            if samples_passed >= num_samples_needed and in_interval \
                    and (energy[i] > silence_tresh or i == len(energy)-1):
                in_interval = False
                if not i == len(energy)-1:
                    intervals.append((interval_start, i-1))
                else:
                    intervals.append((interval_start, i))
                samples_passed = 0
                interval_start = i + 1
                continue
            if energy[i] > silence_tresh:
                in_interval = False
                samples_passed = 0
                interval_start = i + 1
                continue
            if energy[i] < silence_tresh:
                in_interval = True
                samples_passed += 1
                continue
        return intervals

    @staticmethod
    def point_in_intervals(point, intervals):
        """
        Returns Trues if the point exists in some interval in the list
        """
        for (a, b) in intervals:
            if a <= point <= b:
                return a, b
        return False

    def speech_in_interval(self, a, b):
        """
        Returns True if there is speech in the interval [a, b). Not used
        """
        voicing_treshold = 0.1
        voicing_min_len = 0.100  # seconds
        voicing_min_samples = voicing_min_len * self.params.frame_rate
        pov = self.pitch.pov[a:b]
        voicing = pov[pov > voicing_treshold]
        if len(voicing) > voicing_min_samples:
            return True
        else:
            return False

    def remove_short_phrasecomps(self, pauses):
        phrasecomp_min_len = self.params.mulphras_phrasecomp_min_len
        phrasecomp_min_len_smpls = phrasecomp_min_len * self.params.frame_rate
        pitch = self.pitch
        # edges will be removed later
        pause_edge1 = Pause((0, 2), is_long_pause=True, weight=999999)
        pause_edge2 = Pause((len(pitch.f0_log)-2, len(pitch.f0_log)), is_long_pause=True, weight=999999)
        pauses.insert(0, pause_edge1)
        pauses.append(pause_edge2)
        # to-do: this should be reconsidered
        # to-do: search for the overall best candidate to remove
        pauses_ok = False
        while not pauses_ok:
            pauses_ok = True
            for idx in range(len(pauses)):
                if idx-1 >= 0 and idx < len(pauses):
                    (a0, b0) = pauses[idx-1].interval
                    (a1, b1) = pauses[idx].interval
                    boundary0 = (b0+a0)/2
                    boundary1 = (b1+a1)/2
                    interval_left_len = boundary1 - boundary0
                    if interval_left_len < phrasecomp_min_len_smpls:
                        # remove one of the pauses
                        pauses_ok = False
                        if pauses[idx-1].weight < pauses[idx].weight:
                            pauses.pop(idx-1)
                        else:
                            pauses.pop(idx)
                            break
                if idx+1 < len(pauses):
                    (a1, b1) = pauses[idx].interval
                    (a2, b2) = pauses[idx+1].interval
                    boundary1 = (b1+a1)/2
                    boundary2 = (b2+a2)/2
                    interval_right_len = boundary2 - boundary1
                    if interval_right_len < phrasecomp_min_len_smpls:
                        # remove one of the pauses
                        pauses_ok = False
                        if pauses[idx].weight < pauses[idx+1].weight:
                            pauses.pop(idx)
                        else:
                            pauses.pop(idx+1)
                            break
        if pause_edge1 in pauses:
            pauses.remove(pause_edge1)
        if pause_edge2 in pauses:
            pauses.remove(pause_edge2)
        return pauses

    def compute(self):
        pitch = self.pitch
        wave = self.wave
        params = self.params
        paths = self.paths

        if params.phrase_type == 'none':
            zeros = np.zeros(len(pitch.f0_log))
            zero_phrase_atom = Atom(zeros, self.params.frame_rate)
            return Phrase(curve=zeros, atoms=zero_phrase_atom, pauses=[])

        f0_log = pitch.f0_log
        fs = params.frame_rate
        high_freq = params.mulphras_high_freq

        # Compute low-passed phrase
        nyq = 0.5 * fs
        high = high_freq / nyq
        b, a = sp.signal.butter(6, high, btype='lowpass')
        f0_log_lp = sp.signal.filtfilt(b, a, f0_log, padtype='odd', padlen=len(f0_log)-1)

        # Find pauses (silent intervals)
        silence_tresh = params.mulphras_silence_tresh
        min_len = params.mulphras_min_len
        silent_intervals = self.find_silent_intervals(min_len, silence_tresh, medfilt=5)
        # filter silent intervals at begin and end
        for silence in silent_intervals:
            (a, b) = silence
            if a == 0 or b == len(pitch.f0_log_filt)-1:
                silent_intervals.remove(silence)

        # Find minima and remove some / align others
        lp_minima_all = sp.signal.argrelmin(f0_log_lp)

        # a little patch...
        if not params.multiphrase:
            lp_minima_all = [[]]

        lp_min_pauses = []
        for minimum in lp_minima_all[0]:
            # if is already in silence it's ok
            interval_found = self.point_in_intervals(minimum, silent_intervals)
            if interval_found:
                new_pause = Pause(interval_found, is_lp_min=True, lp_min_location=minimum)
                lp_min_pauses.append(new_pause)
                continue
            # if is not, find closest silence and align the minimum to the silence
            # to-do: look for the second closest silence and see if is much longer
            silence_max_dist = params.mulphras_silence_max_dist
            silence_max_dist_samples = np.int(np.ceil(silence_max_dist * params.frame_rate))
            interval_found = None
            for i in range(silence_max_dist_samples):
                interval_found = self.point_in_intervals(minimum-i, silent_intervals)
                if interval_found:
                    break
                interval_found = self.point_in_intervals(minimum+i, silent_intervals)
                if interval_found:
                    break
            if interval_found:
                new_pause = Pause(interval_found, is_lp_min=True, lp_min_location=minimum)
                lp_min_pauses.append(new_pause)

        # If there is a long silence add it to the intervals
        # to-do: somehow compute longpause_min_len individually for the speaker
        longpause_min_len = params.mulphras_longpause_min_len
        long_pauses = []
        for silence in silent_intervals:
            num_samples_needed = longpause_min_len * params.frame_rate
            (a, b) = silence
            if (b-a) > num_samples_needed:
                new_pause = Pause(silence, is_long_pause=True)
                long_pauses.append(new_pause)

        # Now we have lp_min_pauses and long_pauses. They should be merged
        pauses = Pause.merge_sorted_pauses(lp_min_pauses, long_pauses)

        # Compute pause weight
        # to-do: this may go as getter in the Pause class
        for pause in pauses:
            (a, b) = pause.interval
            length = b-a
            is_long = 1 if pause.is_long_pause else 0
            is_lp_min = 1 if pause.is_lp_min else 0
            pause_weight = length * (is_long + is_lp_min)
            pause.weight = pause_weight

        # See if some phrase component is shorter than the threshold
        pauses = self.remove_short_phrasecomps(pauses)

        boundaries = [0]
        for pause in pauses:
            (a, b) = pause.interval
            boundaries.append((a+b)/2)
        boundaries.append(len(pitch.f0_log))

        phrase_intervals = []
        for i in range(len(boundaries)-1):
            phrase_intervals.append((boundaries[i], boundaries[i+1]))

        # Extract separate phrase component for every interval
        extracted_atoms = []

        # a little patch...
        if not params.multiphrase:
            phrase_intervals = [(0, len(pitch.f0_log))]

        phrase_intervals_padded = []
        post_pad_sec = params.mulphras_post_pad_sec
        post_pad = post_pad_sec * params.frame_rate
        for idx, (a, b) in enumerate(phrase_intervals):
            pitch_padded = copy.deepcopy(pitch)
            wave_padded = copy.deepcopy(wave)
            # allow overlapping for all except the last
            if idx != len(phrase_intervals)-1:
                pitch_padded.f0[b:] = 0
                pitch_padded.f0_log[b:] = 0
                pitch_padded.f0_log_filt[b:] = 0
                pitch_padded.pov[b:] = 0
                pitch_padded.weight[b:] = 0
                wave_padded.wave[b:] = 0
                wave_padded.energy[b:] = 0
                b += post_pad
            phrase_intervals_padded.append((a, b))
            # split everything
            pitch_part = PitchCurve(pitch_padded.f0[a:b], pitch_padded.f0_log[a:b],
                                    pitch_padded.f0_log_filt[a:b], pitch_padded.pov[a:b],
                                    pitch_padded.weight[a:b])
            wave_part = Wave(wave_padded.fs, wave_padded.wave[a:b], wave_padded.filepath,
                             wave_padded.energy[a:b])

            atom = FixPosPhraseExtractor(pitch_part, wave_part, params, paths).compute()
            atom.position += a
            extracted_atoms.append(atom)

        # If there is just one phrase atom end here
        if len(extracted_atoms) == 1:
            curve = extracted_atoms[0].curve
            return Phrase(curve=curve, atoms=extracted_atoms, pauses=[])

        # Find intersections, align to them
        intersections = []
        phrase_len = len(pitch.f0_log)
        for i in range(len(extracted_atoms)-1):
            (a1, b1) = phrase_intervals_padded[i]
            (a2, b2) = phrase_intervals_padded[i+1]
            atom1 = extracted_atoms[i].curve
            atom1_padded = np.pad(atom1, (a1, phrase_len-b1), 'constant')
            atom2 = extracted_atoms[i+1].curve
            atom2_padded = np.pad(atom2, (a2, phrase_len-b2), 'constant')
            sub = atom1_padded - atom2_padded
            diffs_sign = np.where(np.diff(np.sign(sub)))[0]
            diff_found = False
            for diff_sign in diffs_sign:
                if atom1_padded[diff_sign] != 0 and atom2_padded[diff_sign] != 0:
                    diff_found = True
                    intersection = diff_sign
                    break
            if not diff_found:
                raise Exception('Cannot find intersection between two phrases')

            if params.mulphras_debugging:
                plt.plot(atom1_padded)
                plt.plot(atom2_padded)
                plt.plot(intersection, atom1_padded[intersection], 'ro')
                plt.show()
            intersections.append(intersection)

        curve = np.zeros(phrase_len)
        for i in range(len(extracted_atoms)):
            atom = np.copy(extracted_atoms[i].curve)
            (a, b) = phrase_intervals_padded[i]
            atom_padded = np.pad(atom, (a, phrase_len-b), 'constant')
            if i == 0:
                atom_start = -1
            if 0 < i < len(extracted_atoms)-1:
                atom_start = intersections[i-1]
            if i == len(extracted_atoms)-1:
                atom_start = intersections[-1]

            atom_end = intersections[i] if i != len(extracted_atoms)-1 else phrase_len-1

            # +1 so the sample from the left is included.
            # the phrase should be smoother at the end, this avoids spikes
            atom_padded[:atom_start+1] = 0
            atom_padded[atom_end+1:] = 0

            curve += atom_padded

        if params.mulphras_debugging:
            plt.plot(curve, 'r')
            plt.plot(intersections, curve[intersections], 'ro')
            plt.show()

        # Return phrase intervals list
        phrase = Phrase(curve=curve, pauses=pauses, atoms=extracted_atoms)
        return phrase

class FixPosPhraseExtractor:
    def __init__(self, pitch, wave, params, paths):
        self.pitch = pitch
        self.wave = wave
        self.params = params
        self.paths = paths

    @staticmethod
    def make_dictionary(atom_type, k, thetas, theta_up, fs, half):
        """
        Given k, and a list of thetas and frameRate,
        generates a list of gamma atoms.
        """
        atoms = []
        for theta in thetas:
            if atom_type == 'gamma':
                new_atom = GammaAtom(k, theta, fs)
            elif atom_type in ('mutant_gamma', 'fix_pos'):
                new_atom = mutantGammaAtom(k, theta, theta_up, fs)

            if half:
                # use only 2nd half of phrase atom:
                # agjoreski: this will be better done with HalfGammaAtom class
                max_in = np.argmax(new_atom.curve)  # index of maximum element - might use it to try half atoms
                new_atom.curve = new_atom.curve[max_in:]
                # this is not needed anymore (because of the limited shifting):
                # new_atom.curve = new_atom.curve - min(new_atom.curve)  # bring it down to 0
            atoms.append(new_atom)
        return atoms

    def compute(self):
        f0_log = self.pitch.f0_log
        weight = self.pitch.weight
        # This gives us the basis phrase atoms:
        dictionary = self.make_dictionary(self.params.phrase_type, 2, self.params.phrase_thetas_range,
                                                     self.params.theta_up, self.params.frame_rate, self.params.half_phrase)

        filt_energy = sp.signal.medfilt(self.wave.energy, 5)
        max_wcorr_atom = None
        max_wcorr = None

        # find maximum pitch within given time after pov_norm goes above 0.5 - this is where we will place the atom
        if self.params.fix_pos_ref_func == 'annotations':
            annotations_file = self.paths.textgrid_folder + self.wave.basename + '.TextGrid'
            tiers = TextGrid(open(annotations_file).read())
            for tier in tiers:
                if tier.nameid == 'emphasis':
                    for i, row in enumerate(tier.simple_transcript):
                        if i == 0:
                            (start, end, label) = row
                            pitch_max_in = np.float(end)  # in time
                            pitch_max_in = int(pitch_max_in * self.params.frame_rate)  # in samples
                            print 'pitch_max_in =', pitch_max_in
                        if i == tier.size-1:
                            (start, end, label) = row
                            pitch_max_out = np.float(start)  # in time
                            pitch_max_out = int(pitch_max_out * self.params.frame_rate)  # in samples
                            print 'pitch_max_out =', pitch_max_out

        elif self.params.fix_pos_ref_func not in {'first', 'both'}:
            if self.params.fix_pos_ref_func == 'pov':
                ref_func = self.pitch.pov
            elif self.params.fix_pos_ref_func == 'energy':
                ref_func = self.wave.energy
            elif self.params.fix_pos_ref_func == 'weight':
                ref_func = self.pitch.weight
            elif self.params.fix_pos_ref_func == 'both':
                ref_func = self.pitch.weight

            ref_func_norm = (ref_func - np.min(ref_func)) / (np.max(ref_func) - np.min(ref_func))
            ref_func_start = np.where(ref_func_norm > self.params.fix_pos_ref_func_th)[0]
            ref_func_start = int(ref_func_start[0])
            if self.params.fix_pos_ref_func_t:  # if not zero search the time window after threshold
                pitch_max_in = np.argmax(f0_log[ref_func_start:
                                                ref_func_start+int(self.params.fix_pos_ref_func_t * self.params.frame_rate)])
                pitch_max_in += ref_func_start
            else:
                pitch_max_in = ref_func_start

            pitch_max_out = np.where(ref_func_norm > self.params.fix_pos_ref_func_end_th)[0]
            pitch_max_out = int(pitch_max_out[-1])
        else:
            # if first or both
            # first analyze energy
            ref_func = self.wave.energy

            ref_func_norm = (ref_func - np.min(ref_func)) / (np.max(ref_func) - np.min(ref_func))
            ref_func_start = np.where(ref_func_norm > self.params.fix_pos_ref_func_th)[0]
            ref_func_start = int(ref_func_start[0])
            if self.params.fix_pos_ref_func_t:  # if not zero search the time window after threshold
                pitch_max_in1 = np.argmax(f0_log[ref_func_start:
                                          ref_func_start+int(self.params.fix_pos_ref_func_t * self.params.frame_rate)])
                pitch_max_in1 += ref_func_start
            else:
                pitch_max_in1 = ref_func_start
            print 'Energy start in = ', pitch_max_in1

            pitch_max_out1 = np.where(ref_func_norm > self.params.fix_pos_ref_func_end_th)[0]
            pitch_max_out1 = int(pitch_max_out1[-1])

            # second pov
            ref_func = self.pitch.pov

            ref_func_norm = (ref_func - np.min(ref_func)) / (np.max(ref_func) - np.min(ref_func))
            ref_func_start = np.where(ref_func_norm > self.params.fix_pos_ref_func_th)[0]
            ref_func_start = int(ref_func_start[0])
            if self.params.fix_pos_ref_func_t:  # if not zero search the time window after threshold
                pitch_max_in2 = np.argmax(f0_log[ref_func_start:
                                          ref_func_start+int(self.params.fix_pos_ref_func_t * self.params.frame_rate)])
                pitch_max_in2 += ref_func_start
            else:
                pitch_max_in2 = ref_func_start
            print 'Pov start in = ', pitch_max_in2

            pitch_max_out2 = np.where(ref_func_norm > self.params.fix_pos_ref_func_end_th)[0]
            pitch_max_out2 = int(pitch_max_out2[-1])

            if self.params.fix_pos_ref_func == 'first':
                pitch_max_in = min(pitch_max_in1, pitch_max_in2)
                pitch_max_out = max(pitch_max_out1, pitch_max_out2)
            elif self.params.fix_pos_ref_func == 'both':
                pitch_max_in = max(pitch_max_in1, pitch_max_in2)
                pitch_max_out = min(pitch_max_out1, pitch_max_out2)

        phrase_max_out = pitch_max_out

        if self.params.fix_pos_ref_func_end_t:
            phrase_max_out = pitch_max_out - int(self.params.fix_pos_ref_func_end_t * self.params.frame_rate)

        start_time = time.time()

        # shorten f0_log and weight if we are discarding phrase_max_out
        if self.params.fix_pos_ref_func_end:
                f0_log = f0_log[: phrase_max_out]
                weight = weight[: phrase_max_out]

        print "Iterating gammas ",
        for atom in dictionary:
            print "#",
            sys.stdout.flush()   # don't buffer output, print hash immediately

            atom_max_in = np.argmax(atom.curve)  # index of maximum element in atom
            #                                      this we center with first voiced position in f0

            # make the atom the same length with f0
            atom_pad = atom.curve
            if atom_max_in > pitch_max_in:
                atom_pad = atom_pad[atom_max_in - pitch_max_in:]
            elif atom_max_in < pitch_max_in:
                atom_pad = np.pad(atom_pad, (pitch_max_in - atom_max_in,), 'constant')

            if len(atom_pad) > len(f0_log):
                atom_pad = atom_pad[:len(f0_log)]
            elif len(atom_pad) < len(f0_log):
                atom_pad = np.pad(atom_pad, (0, len(f0_log) - len(atom_pad)), 'constant')

            if self.params.fix_pos_ref_func_end:
                atom_pad = atom_pad[: phrase_max_out]

            atom_wcorr = np.sum(f0_log * atom_pad * weight)
            down1 = np.sqrt(np.sum(weight * f0_log**2))
            down2 = np.sqrt(np.sum(atom_pad**2 * weight))
            atom_wcorr = atom_wcorr / (down1*down2)

            if max_wcorr is None or atom_wcorr > max_wcorr:
                max_wcorr = atom_wcorr
                max_wcorr_atom = atom
                max_wcorr_atom.position = pitch_max_in - atom_max_in

            ###################################

        print "done in %s seconds" % (time.time() - start_time)

        print "max(wCorr) = %f" % max_wcorr

        f0_log = self.pitch.f0_log
        print 'phrase atom theta is ', max_wcorr_atom.theta
        print 'phrase atom position is ', max_wcorr_atom.position

        phrasecomp = max_wcorr_atom.get_padded_curve(len(f0_log))

        if self.params.phrase_type == 'fix_pos':
            if self.params.fix_pos_ref_func_end:
                phrasecomp /= linalg.norm(phrasecomp[pitch_max_in: phrase_max_out])
                phrasecomp *= np.dot(phrasecomp[pitch_max_in: phrase_max_out], f0_log[pitch_max_in: phrase_max_out])
            else:
                phrasecomp /= linalg.norm(phrasecomp[pitch_max_in:])
                phrasecomp *= np.dot(phrasecomp[pitch_max_in:], f0_log[pitch_max_in:])
        else:
            phrasecomp /= linalg.norm(phrasecomp)
            phrasecomp *= np.dot(phrasecomp, f0_log)

        print pitch_max_in, pitch_max_out
        phrase = Atom(phrasecomp, self.params.frame_rate,
                      pitch_max_in=pitch_max_in, pitch_max_out=pitch_max_out, phrase_max_out=phrase_max_out)

        # take out < pitch_max_in
        if self.params.weight_type == 'binary_both' or self.params.weight_binary:
            self.pitch.weight[: pitch_max_in] = 0

        return phrase
