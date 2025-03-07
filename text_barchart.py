#!/usr/bin/env python

import os
import sys

# input via stdin or sys.argv[1]
# bar values are last column
# Previous columns are labels
# " ".join(labels)


lines = []
skip = False
if ( not skip ):
    for line in sys.stdin.readlines():
       line = line.strip()
       if (len(line) == 0):
           continue
       lines.append(line)


if ( len(sys.argv) > 1):
    with open(sys.argv[1]) as f:
        for line in f:
            line = line.strip()
            if (len(line) == 0):
                continue
            lines.append(line)


all_sp = []

for line in lines:
    line = line.strip()
    if (len(line) == 0):
        continue
    sp = line.split()

    all_sp.append(sp)


max_data = 0
max_length = 0

for sp in all_sp:
    if (len(sp) > 1):
        length = 0
        if (len(sp) > 1):
            length = len(" ".join(sp[:-1]))
        max_length = max(max_length, length)

        max_data = max(float(sp[-1]), max_data)

stars_per_data = 100/max_data


the_format = "%%%is : %%s"%max_length


for sp in all_sp:
    to_print = ""
    if (len(sp) > 1):
        print(the_format%(" ".join(sp[:-1]), "*"*int(float(sp[-1])*stars_per_data)))



