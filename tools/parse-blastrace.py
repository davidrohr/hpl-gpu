#!/usr/bin/env python

from __future__ import division

import sys
import csv

class SetOfOps:
    def __init__(self):
        self.ops = []
        self.n_ops = 0
        self.time_total = 0.0

    def append(self, op):
        """ Append a single BLAS operation. *op* must be a dict """
        self.ops.append(op)
        self.n_ops += 1
        self.time_total += float(op['TIME'])

    def __str__(self):
        return "%d BLAS operations; %f ns total runtime" % (self.n_ops, self.time_total)

    def bin(self, key):
        bins = {}
        for op in self.ops:
            val = op[key]
            if val not in bins:
                bins[val] = SetOfOps()
            bins[val].append(op)
        return bins

    def filter(self, key, value):
        """ Filter all OPs for given Key/Value pair """
        match_ops = SetOfOps()
        for op in self.ops:
            if key in op:
                if op[key] == value:
                    match_ops.append(op)
        return match_ops



if __name__ == "__main__":
    files = sys.argv[1:]

    if len( files ) < 1:
        print "Please specify some input files generated using BLAS-trace code!"
        exit()

    all_ops = SetOfOps()
    for filename in files:
        for line in csv.reader( open( filename ) ):
            op = {}
            op['FUNC'] = line[0]
            for kv in line[1:]:
                key, val = kv.split("=")
                op[key] = val
            all_ops.append(op)

    # all ops in all_ops
    print "== Total statistics =="
    print "  %d BLAS function calls; %f ms total BLAS runtime" % (all_ops.n_ops, all_ops.time_total/1000)
    for func, f_ops in all_ops.bin('FUNC').items():
        print
        print "== %s ==" % func
        print 
        print "   %d calls to %s; %f ms function runtime (%4.1f %%)" % \
            (f_ops.n_ops, func, f_ops.time_total/1000, 100*f_ops.time_total/all_ops.time_total)
        print
        for n, n_ops in f_ops.bin('N').items():
            print "      with N=%d: %5d calls, %8.1f ms runtime (%4.1f %% from total)" % \
                (int(n), n_ops.n_ops, n_ops.time_total/1000, 100*n_ops.time_total/all_ops.time_total)

