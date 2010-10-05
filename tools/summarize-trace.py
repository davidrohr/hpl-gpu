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

class Op:
    def __init__(self, func, time):
        self.func = func
        self.total_calls = 1
        self.total_time = time

    def append(self, time):
        self.total_calls += 1
        self.total_time += time

class SetOfOps:
    def __init__(self):
        self.ops = {}
        self.n_ops = 0
        self.time_total = 0

    def append(self, func, time):
        """ Append a single Op operation. *op* must be a dict """
        self.n_ops += 1
        self.time_total += time
        try:
            self.ops[ func ].append( time )
        except KeyError:
            self.ops[ func ] = Op( func, time );

    def __str__(self):
        return "%d total invocations; %d ns total runtime" % (self.n_ops, self.time_total)

    def getOps(self):
        return self.ops.values()


if __name__ == "__main__":
    files = sys.argv[1:]

    if len( files ) < 1:
        print "Please specify some input files generated using BLAS-trace code!"
        exit()

    all_ops = SetOfOps()
    for filename in files:
        print "Now parsing " + filename + "."
        for line in csv.reader( open( filename ) ):
            func = line[0]
            time = -1
            for kv in line[1:]:
                try: # workaround for broken trace files
                    key, val = kv.split("=")
                except ValueError:
                    # handle throughput numbers and everything else that uses ≅.
                    try:
                        key, val = kv.split("≅")
                    except ValueError:
                        print "Trace output for function " + func + " seems to be broken."
                        continue
                    continue
                if key == 'TIME':
                    time = int(val)
            if time == -1:
                print "Trace output for function " + func + " does not contain a TIME value."
                continue
            all_ops.append(func,time)

    # all ops in all_ops
    print "== Total statistics =="
    print
    print "  %d  function calls; %f ms total  runtime" % (all_ops.n_ops, float(all_ops.time_total)/1000)
    print
    for op in all_ops.getOps():
        print "   %d calls to %s; %f ms function runtime (%4.1f %%)" % \
            (op.total_calls, op.func, float(op.total_time)/1000, 100*float(op.total_time)/float(all_ops.time_total))

