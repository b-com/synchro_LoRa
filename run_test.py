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

import lora.chirp
import lora.synchro

import argparse
import os
from datetime import datetime
import numpy as np
import sys
import random


def is_file_empty(file_path):
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0


def run_tests(sf, snr, count, output_file):

    sync = lora.synchro.synchro(sf)

    osr = 2
    payload_length = 16

    file = open(output_file, "a")
    if is_file_empty(output_file):
        file.write(
            "Date;SF;SNR;Delay;Offset;ExDelay;ExOffset;SynDelay;SynOffset\n")

    for loop in range(count):

        delay = 5 + random.randint(0, 2**sf)
        freq_offset = random.uniform(-10e3, 10e3)

        signal = lora.chirp.gen_lora_signal(
            sf, osr, payload_length, delay, freq_offset, snr)
        signal_8x = lora.chirp.interpolate(signal, 4, 1, 0)

        ex_delay, ex_offset = sync.exhaustive(signal_8x)
        syn_delay, syn_offset = sync.synchronize(signal_8x)

        res = "{};{};{};{};{};{};{};{};{}\n".format(datetime.now(
        ), sf, snr, delay, freq_offset, ex_delay, ex_offset, syn_delay, syn_offset)
        print(res, end='')
        file.write(res)
        file.flush()

    file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='LoRa synchro performance testing tool')
    parser.add_argument('--sf', action='store', default=7)
    parser.add_argument('--snr', action='store', default=0)
    parser.add_argument('--count', action='store', default=1000)
    parser.add_argument('--output-file', action='store', default='results.csv')
    args = parser.parse_args()

    run_tests(int(args.sf), float(args.snr), int(
        args.count), args.output_file)
