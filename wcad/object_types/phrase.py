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


class Phrase:
    def __init__(self, curve=None, minima_lp=None, pauses=None, atoms=None, intersections=None):
        """
        :param curve: sum of all phrase components
        :param pauses: list of Pauses
        :param atoms: list of the separate Atoms, overlapping
        :param intersections: list of point where the overlapping of the atoms takes place
        :return: None
        """
        self.curve = curve
        self.minima_lp = minima_lp
        self.pauses = pauses
        self.atoms = atoms
        self.intersections = intersections

    def get_lp_intervals(self):
        phrase_lp_intervals = []
        for i in range(len(self.minima_lp)-1):
            phrase_lp_intervals.append((self.minima_lp[i], self.minima_lp[i+1]))
        return phrase_lp_intervals


class Pause:
    def __init__(self, interval=None, is_lp_min=False, is_long_pause=False,
                 lp_min_location=None, weight=None):
        """
        :param interval: tuple [start, end) of the silence
        :param is_lp_min: is it a minumum of a low-passed pitch curve
        :param is_long_pause: is it a long pause
        :param lp_min_location: the exact location of the low-pass minimum
        """
        self.interval = interval
        self.is_lp_min = is_lp_min
        self.is_long_pause = is_long_pause
        self.lp_min_location = lp_min_location
        self.weight = weight

    @staticmethod
    def pauses_overlapping(pause1, pause2):
        (a1, b1) = pause1.interval
        (a2, b2) = pause2.interval
        if (a1 <= b2) and (b1 >= a2):
            return True
        return False

    @staticmethod
    def insert_pause(new_pause, list_pauses):
        (a1, b1) = new_pause.interval
        for idx in range(len(list_pauses)):
            (a2, b2) = list_pauses[idx].interval
            if a1 <= a2:
                # if the pauses are overlapping merge them and check if the folowing should be merged
                while idx < len(list_pauses) \
                        and Pause.pauses_overlapping(new_pause, list_pauses[idx]):
                    new_pause = Pause.merge_pauses(new_pause, list_pauses[idx])
                    list_pauses.remove(list_pauses[idx])
                    idx += 1
                # add the newone + all the remaining from list
                list_pauses.insert(idx, new_pause)

    @staticmethod
    def merge_pauses(pause1, pause2):
        (a1, b1) = pause1.interval
        (a2, b2) = pause2.interval
        start = np.min([a1, a2])
        end = np.max([b1, b2])
        is_lp_min = pause1.is_lp_min or pause2.is_lp_min
        is_long_pause = pause1.is_long_pause or pause2.is_long_pause
        lp_min_location = None
        if pause1.lp_min_location is not None:
            lp_min_location = pause1.lp_min_location
        if pause2.lp_min_location is not None:
            lp_min_location = pause2.lp_min_location

        merged_pause = Pause(interval=(start, end), is_lp_min=is_lp_min,
                             is_long_pause=is_long_pause, lp_min_location=lp_min_location)
        return merged_pause

    @staticmethod
    def merge_sorted_pauses(list1, list2):
        if len(list1) == 0:
            return list2
        if len(list2) == 0:
            return list1
        (a1, b1) = list1[0].interval
        (a2, b2) = list2[0].interval
        # Always start with the first interval
        if a2 < a1:
            temp_list = list1
            list1 = list2
            list2 = temp_list
        # Insert every interval from list2 to list1
        for pause in list2:
            Pause.insert_pause(pause, list1)

        return list1