#!/usr/bin/env python

import sys,os
import pandas as pd
import numpy as np

fnames = sys.argv[1:]

data = [ np.loadtxt(x, delimiter="\t", dtype=str) for x in fnames]

cnt = np.sum(np.array([x[:,3].astype(int) for x in data]),axis=0)

for i,r in enumerate(data[0]):
    print('\t'.join(list(r[:-1]) + [str(cnt[i])]))

