#!/usr/bin/env python
# vim: set fileencoding=utf-8 :

#
# The source code is property of the Frankfurt Institute for Advanced Studies
# (FIAS). None of the material may be copied, reproduced, distributed,
# republished, downloaded, displayed, posted or transmitted in any form or by
# any means, including, but not limited to, electronic, mechanical,
# photocopying, recording, or otherwise, without the prior written permission
# of FIAS.
#
# Authors:
# David Rohr (drohr@jwdt.org)
# Matthias Bach (bach@compeng.uni-frankfurt.de)
# Matthias Kretz (kretz@compeng.uni-frankfurt.de)
# Jörg Bornschein (bornschein@fias.uni-frankfurt.de)
# 

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
                try: # workaround for broken trace files
                    key, val = kv.split("=")
                    op[key] = val
                except ValueError:
                    # handle throughput numbers and everything else that uses ≅.
                    try:
                        key, val = kv.split("≅")
                        op[key] = val
                    except ValueError:
                        print "Trace output for op " + line[0] + " seems to be broken."
                        pass
                    pass
            all_ops.append(op)

    # all ops in all_ops
    print "== Total statistics =="
    print
    print "  %d BLAS function calls; %f ms total BLAS runtime" % (all_ops.n_ops, all_ops.time_total/1000)
    for func, f_ops in all_ops.bin('FUNC').items():
        print
        print "== %s ==" % func
        print 
        print "   %d calls to %s; %f ms function runtime (%4.1f %%)" % \
            (f_ops.n_ops, func, f_ops.time_total/1000, 100*f_ops.time_total/all_ops.time_total)
        print
        try:
            mnk_dict = f_ops.bin('N')       # Bin according to N

            # Sort bins according to aggregate runtime
            sorted_mnk = sorted(mnk_dict, key=lambda k: mnk_dict[k].time_total, reverse=True)
            for mnk in sorted_mnk[:10]:
                ops = mnk_dict[mnk]      # All ops in this bin
                percentage = 100*ops.time_total/all_ops.time_total

                print "      with N=%8s: %6d calls, %8.1f ms runtime (%4.1f %% of total)" % \
                    (mnk, ops.n_ops, ops.time_total/1000, percentage)
        except KeyError:
            # Handle functions that do not have a parameter N
            # In that case we only care about the summary and not about the sub-binning
            pass

