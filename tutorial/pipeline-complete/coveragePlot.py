#!/usr/bin/env python

import sys, os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use('Agg')

fin = []
title, fout = sys.argv[1:3]
fin = sys.argv[3:]

if type(fin) == str:
    fin = [fin]

data = [np.loadtxt(f, delimiter="\t", dtype=str)[:,-1].astype(int) 
        for f in fin]

coverage = np.array([[100*sum(x[x > k])/sum(x) 
                      for k in range(1,41)] for x in data])

summary = np.mean(coverage, axis=0)

plt.plot(coverage.transpose(), color='b')
if coverage.shape[0] > 1:
    plt.plot(summary, 'r')
plt.grid(1)
plt.xticks(list(range(0,45,5)), [str(x) + 'x' for x in range(0,45,5)])
plt.ylabel('coverage %')
plt.title(title)
plt.savefig(fout)

