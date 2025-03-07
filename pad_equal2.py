#!/usr/bin/env python

import os
import sys
import numpy as np


pattern = sys.argv[1]
total_length = int(sys.argv[2])
sequence_file = sys.argv[3]

sequences = []
rest_of_line = []
with open(sequence_file) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0 ):
            continue
        sp = line.split()
        if ( len(sp) > 1 ):
            rest_of_line.append(" ".join(sp[1:]))
            sequences.append(sp[0])
        else:
            sequences.append(line)

if ( len(rest_of_line) > 0 ):
    assert(len(rest_of_line) == len(sequences))

mul = int(total_length / len(pattern) + 1)
full_pad = pattern * mul


for iseq, sequence in enumerate(sequences):

    missing = total_length - len(sequence)
    assert( missing >= 0 )

    missing_left = int(np.floor(missing/2))
    missing_right = int(np.ceil(missing/2))

    left_pad = ""
    if ( missing_left > 0 ):
        left_pad = full_pad[-missing_left:]
    right_pad = ""
    if ( missing_right > 0 ):
        right_pad = full_pad[:missing_right]

    final = left_pad + sequence + right_pad
    assert(len(final) == total_length)
    if ( len(rest_of_line) > 0):
        print(final, rest_of_line[iseq])
    else:
        print(final)



