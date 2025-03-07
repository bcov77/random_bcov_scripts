#!/usr/bin/env python

import os
import sys
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("max_length", type=int)
parser.add_argument("sequence_file", type=str)
parser.add_argument("--opening_pad", type=str, default="GGS")
parser.add_argument("--only_gs", action='store_true')
parser.add_argument("--no_g", action='store_true')
args = parser.parse_args(sys.argv[1:])


aa_probs = {
    "G":10,
    "S":10,
    "N":0,
    "Q":5,
    "E":1,
    "D":1,
    "A":1,
    "T":1,
}

if ( args.only_gs ):
    for key in aa_probs:
        if ( key not in "GS" ):
            aa_probs[key] = 0

if ( args.no_g ):
    aa_probs['G'] = 0

expanded_list = ""
aa_sum = sum(aa_probs.values())

for letter in aa_probs:
    this_many = int(aa_probs[letter] / aa_sum * 1000 + 0.5)
    expanded_list += letter * this_many

expanded_list = list(expanded_list)

def generate_random_pad(size):
    passing = False
    while not passing:
        rands = np.random.random(size)
        pad = ""
        for rand in rands:
            pad += expanded_list[int(rand*len(expanded_list))]

        passing = True
        for i in range(size//5):
            start = i*5
            end = (i+1)*5
            if ( end > size ):
                break
            sub = pad[start:end]
            num_gs = sub.count("G") + sub.count("S")
            if ( num_gs < 3 ):
                passing = False
                break
        if ( pad.count("E") + pad.count("D") > 2 ):
            passing = False

        if ( pad.count("A") + pad.count("T") > 2 ):
            passing = False

    return pad

total_length = args.max_length
sequence_file = args.sequence_file
opening_pad = args.opening_pad

sequences = []
with open(sequence_file) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0 ):
            continue
        sequences.append(line)

# mul = int(total_length / len(pattern) + 1)
# full_pad = pattern * mul


for sequence in sequences:

    missing = total_length - len(sequence)
    assert( missing >= 0 )

    missing_left = int(np.floor(missing/2))
    missing_right = int(np.ceil(missing/2))



    left_pad = ""
    if ( missing_left > 0 ):
        if ( missing_left > len(opening_pad) ):
            left_pad = str( opening_pad[::-1])

            left_pad = generate_random_pad(missing_left - len(opening_pad)) + left_pad
        else:
            left_pad = str(opening_pad[::-1])[-missing_left:]

    right_pad = ""
    if ( missing_right > 0 ):
        if ( missing_right > len(opening_pad) ):
            right_pad = str( opening_pad)

            right_pad += generate_random_pad(missing_right - len(opening_pad))
        else:
            right_pad = str(opening_pad)[:missing_right]

    final = left_pad + sequence + right_pad
    print(final)
    # print(len(final), total_length)
    assert(len(final) == total_length)



