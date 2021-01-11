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

import numpy as np


def _phi_bar(v):
    u = v % 1
    return(u-1) * u


def gen_chirp(SF=7, k=0, osr=2):
    """
    Generate a chirp signal over a symbol period
      Inputs:
        SF: spreading factor
        k: integer defining the modulated chirp symbol (between 0 and 2^SF-1)
        osr: over sampling ratio compared to chip rate
      Output:
        osr x (2^SF) samples representing the requested modulated
        chirp symbol
    """

    Tc = osr
    t = np.arange(0, 2**(SF) * Tc, 1)
    waveform = np.exp(2j * np.pi * 2**(SF-1) * _phi_bar(2**(-SF) * (t / Tc-k)))
    if k != 0:
        ph = np.angle(waveform[0])
        waveform = waveform*np.exp(-1j*ph)
    return(waveform)


def gen_short_preamble(sf, preamble_len=8, osr=2):
    """
    Generate the first part of the LoRa preamble (preamble_len x 0 symbols)
      Inputs:
        sf: spreading factor
        preamble_len: number of 0 symbols, typically 8
        osr: over sampling ratio compared to chip rate
      Output:
        osr x (2^SF) * preable_len samples, containing preamble_len 0 symbols

    """

    c0 = gen_chirp(SF=sf, k=0, osr=osr)
    return np.tile(gen_chirp(SF=sf, k=0, osr=osr), preamble_len)


def gen_preamble(sf, preamble_len=8, sync_word=0x34, osr=2):
    """
    Generate a complete LoRa preamble
      Inputs:
        sf: spreading factor
        preamble_len: number of 0 symbols, typically 8
        sync_word: LoRa 'sync_word' network type idendifier
        osr: over sampling ratio compared to chip rate
      Output:
        osr x (2^SF) * 12.5 samples, representing a LoRa preamble
    """

    # 8 zero symbols
    preamble = gen_short_preamble(sf, preamble_len, osr)

    # sync word
    preamble = np.hstack((preamble,
                          gen_chirp(SF=sf, k=(2 ** sf) - ((sync_word & 0xF0) >> 1), osr=osr)))
    preamble = np.hstack((preamble,
                          gen_chirp(SF=sf, k=(2 ** sf) - ((sync_word & 0xF) << 3), osr=osr)))

    # 2 inverted symbols
    c0 = gen_chirp(SF=sf, k=0, osr=osr)
    preamble = np.hstack((preamble, np.conj(c0)))
    preamble = np.hstack((preamble, np.conj(c0)))

    # guard
    guardtime = np.conj(c0)[0:osr * 2 ** (sf - 2)]

    preamble = np.hstack((preamble, guardtime))

    return preamble


def apply_frequency_offset(s, freq, osr):
    """
    Apply a frequency offet to the input singal
    Inputs:
        s: input signal
        freq: offset to apply
        osr: over sampling ratio compared to chip rate
    Output:
        offseted signal
    """

    sampling_rate = osr * 125e3
    r = np.copy(s)
    for i in range(len(s)):
        angle = 2 * np.pi * i * freq / sampling_rate
        r[i] = s[i] * np.exp(1j*angle)
    return r


def generate_noise(snr, len):
    """
    Generate Gaussian noise
    Inputs:
        snr: level
        len: requested len
    Outputs:
        len samples noise signal
    """

    power_lin = 10 ** (-snr/10)
    return np.sqrt(power_lin/2) * (np.random.normal(0, 1, len) + 1j*np.random.normal(0, 1, len))


def gen_lora_signal(sf, osr, paylod_len, tau, offset, snr):
    """
    Generate a complete LoRa signal
    Inputs:
        sf: spreading factor
        osr: over sampling ratio compared to chip rate
        payload: frame payload
        preamble_len: number of 0 symbols, typically 8
        tau: delay, signal starting position
        offset: frequency offset
        snr: signal noise ratio
    """

    # generate lora frame
    signal = gen_preamble(sf, osr=osr)
    payload_symbols = np.random.randint(0, 2**sf, size=paylod_len)
    for i in range(paylod_len):
        c = gen_chirp(sf, payload_symbols[i], osr)
        signal = np.hstack((signal, c))

    # apply fequency offset
    signal = apply_frequency_offset(signal, offset, osr)

    # generate and add noise
    noise = generate_noise(snr, tau + len(signal))
    noise[tau:-1] = noise[tau:-1] + signal[0:-1]
    return noise


def interpolate(signal, N, M, offset):
    """
    Lagragian interpolation
    Inputs:
        N: oversampling factor
        M: downsampling factor
        offset: input offset
    Outputs:
        resampled signal 
    """

    factor = float(N) / float(M)
    res = np.zeros(int(len(signal) * factor), dtype="complex64")

    xx = offset

    for i in range(len(res)):
        idx = int(xx)
        x = xx - int(xx) + 1
        y0 = signal[idx-1] if idx-1 >= 0 and idx-1 < len(signal) else 0
        y1 = signal[idx] if idx >= 0 and idx < len(signal) else 0
        y2 = signal[idx + 1] if idx+1 >= 0 and idx+1 < len(signal) else 0
        y3 = signal[idx + 2] if idx+2 >= 0 and idx+2 < len(signal) else 0

        res[i] = 1.0*x*y1*(0.5*x - 1.5)*(1.0*x - 2.0) - 0.5*x*y2*(1.0*x - 3.0)*(x - 1) + 0.166666666666666666666666666666666666 * \
            x*y3*(1.0*x - 2.0)*(x - 1) - 1.0*y0 * \
            (0.33333333333333333333333333333333333333333333*x - 1.0) * \
            (0.5*x - 1.0)*(x - 1)

        xx += 1.0 / factor

    return res
