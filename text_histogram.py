#!/usr/bin/env python

import os
import sys
import numpy as np
import subprocess
import random
from seaborn.utils import iqr

# input via stdin or sys.argv[1]
# all values are used to make histogram

def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return the_stuff[0] + the_stuff[1]


lines = []
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

def freedman_diaconis_bins(a):
    """Calculate number of hist bins using Freedman-Diaconis rule."""
    # From https://stats.stackexchange.com/questions/798/
    a = np.asarray(a)
    if len(a) < 2:
        return 1
    h = 2 * iqr(a) / (len(a) ** (1 / 3))
    # fall back to sqrt(a) bins if iqr is 0
    if h == 0:
        return int(np.sqrt(a.size))
    else:
        return int(np.ceil((a.max() - a.min()) / h))




values = []
for line in lines:
    for item in line.split():
        try:
            values.append(float(item))
        except:
            pass
values = np.array(values)

bins = freedman_diaconis_bins(values)


hist, bin_edges = np.histogram(values, bins)

random_name = ("%.8f"%random.random())[2:] + ".tmp"

try:
    f = open(random_name, "w")
except:
    home = os.path.expanduser("~/")
    tempdir = os.path.join(home, "temp")
    if (not os.path.exists(tempdir) ):
        os.mkdir(tempdir)
    random_name = os.path.join(tempdir, random_name)
    f = open(random_name, "w")


for i in range(len(hist)):
    f.write("%.2f %.3f\n"%(bin_edges[i], hist[i]))

f.close()

print(cmd("/home/bcov/util/text_barchart.py %s"%random_name))

os.remove(random_name)














