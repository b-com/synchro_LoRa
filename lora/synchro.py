# This file is part of the syncro_LoRa distribution (https://github.com/b-com/synchro_LoRa).
# Copyright (c) 2021 Christophe Delacourt.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from . import chirp
import numpy as np


class synchro:
    """
        Lora synchronisation class
    """
    SAMPLING_RATE = 125e3
    GRANULARITY = 64

    def __init__(self, sf):
        self.ref_preamble = np.conj(chirp.gen_short_preamble(sf, 8, 8))
        self.ref_preambles = self.gen_preambles(sf)
        self.win_len = len(self.ref_preamble)
        self.step = self.win_len // self.GRANULARITY

    def gen_preambles(self, sf):
        """ Generate reference preambles """
        preambles = []

        for p in range(8):
            preamble = []

            for k in range(8):
                angle = -2.0 * np.pi * k * p / 8.0
                c = np.cos(angle)
                s = np.sin(angle)

                ch = chirp.gen_chirp(sf, 0, 8)
                for l in range(len(ch)):
                    I = ch[l].real
                    Q = ch[l].imag
                    ch[l] = (I * c - Q * s) + 1j * (I * s + Q * c)

                preamble = np.hstack((preamble, ch))
            preambles.append(np.conj(preamble))

        return preambles

    def max_fft(self, signal, preamble, tau, osr):
        """
            max_fft : TODO see reference paper
        """
        factor = 8 // osr
        dechirp = preamble[0:len(preamble):factor] * \
            signal[tau:tau + len(preamble):factor]
        fft = np.abs(np.fft.fft(dechirp))

        max_i = np.argmax(fft)
        return (max_i, fft[max_i])

    def max_fft8(self, signal, tau, osr):
        """
            max_fft8 : TODO see reference paper
        """
        max_i = 0
        max_v = 0.0
        for p in range(8):
            loc_max_i, loc_max_v = self.max_fft(
                signal, self.ref_preambles[p], tau, osr)
            if loc_max_v > max_v:
                max_v = loc_max_v
                max_i = loc_max_i
        return (max_i, max_v)

    def exhaustive(self, signal_8x):
        """
            exhaustive : TODO see reference paper
        """
        osr = 8
        max_i = 0  # found delay
        max_b = 0  # index of maximum bin
        max_v = 0.0  # maximum value

        # find maximum on the whole window
        for i in range(self.win_len):
            local_max_i, local_max_v = self.max_fft(
                signal_8x, self.ref_preamble, i, 8)
            if local_max_v > max_v:
                max_i = i
                max_v = local_max_v
                max_b = local_max_i

        offset = np.fft.fftfreq(
            self.win_len, 1.0 / (self.SAMPLING_RATE * osr))[max_b]
        return (2 * max_i / osr, offset)

    def phase1(self, signal):
        """
            phase1 : TODO see reference paper
        """
        max_i = 0
        max_v = 0.0
        for p in range(self.GRANULARITY):
            pos = self.step // 2 + p * self.step
            if(pos >= self.win_len):
                continue
            loc_max_i, loc_max_v = self.max_fft8(signal, pos, 2)
            if loc_max_v > max_v:
                max_v = loc_max_v
                max_i = pos
        return (max_i, max_v)

    def phase2(self, signal, phase1_max):
        """
            phase2 : TODO see reference paper
        """

        max_i = 0
        max_v = 0.0
        for p in range(-4, 5):
            pos = phase1_max + p
            if pos < 0 or pos >= self.win_len:
                continue
            loc_max_i, loc_max_v = self.max_fft(
                signal, self.ref_preamble, pos, 2)
            if loc_max_v > max_v:
                max_v = loc_max_v
                max_i = pos
        return (max_i, max_v)

    def phase3(self, signal, phase2_max):
        """
            phase3 : TODO see reference paper
        """

        max_i = 0
        max_v = 0.0
        max_b = 0
        for p in range(-(self.step-4)//16-4, (self.step+5)//16 + 1):
            pos = phase2_max + p * 8
            if pos < 0 or pos >= self.win_len:
                continue
            loc_max_i, loc_max_v = self.max_fft(
                signal, self.ref_preamble, pos, 2)
            if loc_max_v > max_v:
                max_v = loc_max_v
                max_b = loc_max_i
                max_i = pos
        return (max_i, max_b)

    def synchronize(self, signal):
        """
            proposed synchronization algorithm : TODO see reference paper
        """

        p1 = self.phase1(signal)
        p2 = self.phase2(signal, p1[0])
        max_i, max_b = self.phase3(signal, p2[0])
        tau = max_i / 4

        # compute freq
        osr = 2
        factor = 8 // osr
        offset = np.fft.fftfreq(
            len(self.ref_preamble) // factor, 1.0 / (self.SAMPLING_RATE * osr))[max_b]

        return (tau, offset)
