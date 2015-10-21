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

import colorsys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np

from nltk_textgrid.textgrid import *


class ModelPlotterc():
    def __init__(self, wave, model, pitch, params, paths):
        self.wave = wave
        self.model = model
        self.pitch = pitch
        self.params = params
        self.energy = wave.energy
        self.paths = paths

    @staticmethod
    def get_color(ncolor):
        """
        Returns ncolor colors
        """
        for hue in range(ncolor):
            hue = 1. * hue / ncolor
            col = [int(x) for x in colorsys.hsv_to_rgb(hue, 1.0, 230)]
            yield '#{0:02x}{1:02x}{2:02x}'.format(*col)

    def plot(self):

        f0 = self.model.original_pitch
        f0_log_filt = self.pitch.f0_log_filt
        phrase = self.model.phrase
        pov = self.pitch.pov
        weight = self.pitch.weight
        energy = self.energy
        params = self.params
        fs = self.params.frame_rate
        ts = 1. / fs
        t = np.arange(0, len(f0))
        t = t * ts

        scaled_pov = 1 - (pov - np.min(pov))/(np.max(pov) - np.min(pov))  # 1 - is for coloring
        scaled_energy = (energy - np.min(energy))/(np.max(energy) - np.min(energy))
        scaled_weight = 1 - (weight - np.min(weight))/(np.max(weight)-np.min(weight))  # 1 - is for coloring

        plt.figure(figsize=(16.0, 9.0))
        gs = gridspec.GridSpec(3, 1, height_ratios=[4, 4, 1])
        plt.subplots_adjust(hspace=0.35)

        if hasattr(phrase, 'atoms'):
            pitch_max_in = phrase.atoms[0].pitch_max_in
            pitch_max_out = phrase.atoms[-1].pitch_max_out
        else:
            pitch_max_in = phrase.pitch_max_in
            pitch_max_out = phrase.pitch_max_out

        x_start = t[pitch_max_in]
        x_end = t[pitch_max_out]
        print pitch_max_in, pitch_max_out
        print 'x_start = ', x_start
        print 'x_end = ', x_end

        # if there are annotations load them
        if params.plot_textgrids:
            annotations_file = self.paths.textgrid_folder + self.wave.basename + '.TextGrid'
            tiers = TextGrid(open(annotations_file).read())

        # Plot wave + pitch + phrase + reconstruction
        plt.subplot(gs[0])
        plt_ymax = np.max([np.max(f0[pitch_max_in:pitch_max_out]),
                           np.max(phrase.curve[pitch_max_in:pitch_max_out])])
        plt_ymin = np.min([np.min(f0[pitch_max_in:pitch_max_out]),
                           np.min(phrase.curve[pitch_max_in:pitch_max_out])]) - 0.2

        # plot annotations
        if params.plot_textgrids:
            ax = plt.gca()
            for tier in tiers:
                if tier.nameid == 'emphasis':
                    for row in tier.simple_transcript:
                        (start, end, label) = row
                        start = np.float(start)
                        end = np.float(end)
                        # print 'row-start:', start, 'row-end', end, 'row-label', label
                        plt.vlines([start, end], plt_ymin, plt_ymax, color='0.8', lw=2, alpha=0.4)
                        if '+' in label:
                            plt.fill_betweenx([plt_ymin, plt_ymax], x1=start, x2=end,
                                              facecolor='r', alpha=0.4)
                # plot words or phones
                if tier.nameid == 'words':
                    for row in tier.simple_transcript:
                        (start, end, label) = row
                        start = np.float(start)
                        end = np.float(end)
                        # print 'row-start:', start, 'row-end', end, 'row-label', label
                        plt.vlines([start, end], plt_ymin, plt_ymax, color='0.1', lw=2, alpha=0.8)
                        # ax_range = (plt_xmax - plt_xmin)
                        x = 0.5 * (start + end)  # - plt_xmin)/ax_range
                        y = plt_ymin + 0.1 * (plt_ymax - plt_ymin)
                        ax.text(x, y, label, size=12, ha='center', va='center')
                # if tier.nameid == 'phones':
                #     for row in tier.simple_transcript:
                #         (start, end, label) = row
                #         start = np.float(start)
                #         end = np.float(end)
                #         # print 'row-start:', start, 'row-end', end, 'row-label', label
                #         # plt.vlines([start, end], plt_ymin, plt_ymax, color='0.5', lw=1, alpha=0.2)
                #         # ax_range = (plt_xmax - plt_xmin)
                #         x = 0.5 * (start + end)  # - plt_xmin)/ax_range
                #         y = plt_ymin + 0.1 * (plt_ymax - plt_ymin)
                #         ax.text(x, y, label, size=12, ha='center', va='center')

                if tier.nameid == 'lexstress':
                    for row in tier.simple_transcript:
                        (start, end, label) = row
                        start = np.float(start)
                        end = np.float(end)
                        # print 'row-start:', start, 'row-end', end, 'row-label', label
                        # plt.vlines([start, end], plt_ymin, plt_ymax, color='0.7', lw=2, alpha=0.4)
                        if '+' in label:
                            plt.fill_betweenx([plt_ymin, plt_ymax], x1=start, x2=end,
                                              facecolor='lemonchiffon', alpha=0.6)

        # Plot wave in background
        wave_shifted = np.copy(self.wave.wave)
        wave_shifted *= (plt_ymax - plt_ymin)
        wave_shifted += plt_ymin + (plt_ymax-plt_ymin)/2
        t_wave = np.arange(0, len(wave_shifted))
        t_wave = t_wave * (1./ self.wave.fs)
        plt.plot(t_wave, wave_shifted, 'b', alpha=0.2)

        # Plot phrase component(s)
        plt_phrase, = plt.plot(t, phrase.curve, 'r', linewidth=3, alpha=0.8)

        # Plot reconstruction
        reconstruction = self.model.reconstruction
        plt_recons, = plt.plot(t, reconstruction, 'm', linewidth=4, alpha=0.8)

        colmap_step = 1
        if params.plot_pov_color:
            color_scale = scaled_pov
        else:
            color_scale = scaled_weight

        for i in range(0, len(f0) - colmap_step, colmap_step):
            plt.plot(t[i:i+colmap_step+1], f0[i:i+colmap_step+1],
                     color=(0.0, 1-color_scale[i], color_scale[i]),
                     linewidth=3, alpha=0.8)

        # Plot filtered pitch POV-coloured
        if params.use_lowpass_phrase:
            for i in range(0, len(f0_log_filt) - colmap_step, colmap_step):
                plt.plot(t[i:i+colmap_step+1], f0_log_filt[i:i+colmap_step+1],
                         color=(0.0, 1-color_scale[i], color_scale[i]),
                         linewidth=3)

        plt.xlim([x_start, x_end])
        plt.ylim([plt_ymin, plt_ymax])
        plt.title(self.wave.filename, size=11)
        plt.ylabel('log(f0)')

        # Pitch max in and out
        if params.phrase_type == 'fix_pos' and params.plot_pitch_max_in_out:
            plt.axvspan(0, x_start, alpha=0.15, color='black')
            plt.axvspan(x_end, len(t), alpha=0.15, color='black')

        # Plot phrase boundaries
        if hasattr(phrase, 'pauses'):
            for pause in phrase.pauses:
                (a, b) = pause.interval
                plt.axvspan(t[a], t[b], alpha=0.1, color='black')
                if pause.is_lp_min:
                    plt.axvline(x=t[pause.lp_min_location], color='#9eb2ec', alpha=0.7, linewidth=1.5, linestyle='dashed')

            plt_f0 = mpatches.Patch(color='b')
            if params.use_lowpass_phrase:
                plt.legend([plt_f0, plt_f0, plt_phrase, plt_recons],
                           ['f0 & weight', 'f0_lp', 'phrase comp', 'reconstruction'],
                           loc='best', ncol=3, fontsize=9, frameon=False)  # not true if weight is used for plotting
            else:
                plt.legend([plt_f0, plt_phrase, plt_recons],
                           ['f0 & weight', 'phrase comp', 'reconstruction'],
                           loc='best', ncol=3, fontsize=9, frameon=False)

        # Plot atoms
        plt.subplot(gs[1])

        n_atom = len(self.model.atoms)
        color = self.get_color(n_atom)
        for atom in self.model.atoms:
            acolor = next(color)
            plt.plot(t, atom.get_padded_curve(len(f0)), color=acolor, linewidth=2.5)

        # Pitch max in and out
        if params.phrase_type == 'fix_pos' and params.plot_pitch_max_in_out:
            plt.axvspan(0, x_start, alpha=0.15, color='black')
            plt.axvspan(x_end, len(t), alpha=0.15, color='black')

        # Plot phrase boundaries
        if hasattr(phrase, 'pauses'):
            for pause in phrase.pauses:
                (a, b) = pause.interval
                plt.axvspan(t[a], t[b], alpha=0.1, color='black')

        plt.xlim([x_start, x_end])
        plt.ylabel('atom amplitude')

        plt.xlim([x_start, x_end])
        plt.gca().xaxis.grid(True)
        plt.grid()

        # Plot weight
        plt.subplot(gs[2])

        plt_pov, = plt.plot(t, 1 - scaled_pov, '#bf8700')

        plt_energy, = plt.plot(t, scaled_energy, 'm')

        for i in range(0, len(weight), colmap_step):
            plt.plot(t[i:i+1+colmap_step], 1-scaled_weight[i:i+1+colmap_step],
                     color=(0, 1-scaled_weight[i], scaled_weight[i]),
                     linewidth=2)

        # Plot pitch max in and out
        if params.phrase_type == 'fix_pos' and params.plot_pitch_max_in_out:
            plt.axvspan(0, x_start, alpha=0.15, color='black')
            plt.axvspan(x_end, len(t), alpha=0.15, color='black')

        # Plot phrase boundaries
        if hasattr(phrase, 'pauses'):
            for pause in phrase.pauses:
                (a, b) = pause.interval
                plt.axvspan(t[a], t[b], alpha=0.15, color='black')

        plt.xlim([x_start, x_end])
        plt.xlabel('time [s]')
        plt.ylabel('weight')
        plt.yticks([0, 0.5, 1])
        plt.legend([plt_pov, plt_energy], ['pov', 'energy'], loc=1, fontsize=9, frameon=False)

        # Save and/or show
        if params.save_plot:
            savefile = self.paths.figure + '.png'
            plt.savefig(savefile, bbox_inches=0)    # , format='eps', dpi=1000)
            print 'Figure saved as ', savefile

        if params.show_plot:
            plt.show()
