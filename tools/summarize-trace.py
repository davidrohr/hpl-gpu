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
#

from __future__ import division

import sys
import csv

class Func:
    def __init__(self, line):
        self.name = line[0]
        self.invocations = int(line[1])
        self.walltime = int(line[2])
        self.cputime = int(line[3])

    def add(self, other):
        self.invocations += other.invocations
        self.walltime += other.walltime
        self.cputime += other.cputime

class AllFuncs:
    def __init__(self):
        self.funcs = {}
        self.total_invocations = 0
        self.total_walltime = 0
        self.total_cputime = 0

    def append(self, func):
        """ Add data of one function in one file """
        self.total_invocations += func.invocations;
        self.total_walltime += func.walltime;
        self.total_cputime += func.cputime;
        try:
            self.funcs[ func.name ].add( func )
        except KeyError:
            self.funcs[ func.name ] = func;

    def __str__(self):
        return "%d total invocatioms; %f s total runtime; %f s total cputime; cputime / walltime = %f %%" \
            % (self.total_invocations, self.total_walltime / 1e6, self.total_cputime / 1e6, self.total_cputime / self.total_walltime * 100 )

    def getFuncs(self):
        return self.funcs.values()


if __name__ == "__main__":
    files = sys.argv[1:]

    if len( files ) < 1:
        print "Please specify some input files generated using BLAS-trace code!"
        exit()

    all_funcs = AllFuncs();

    for filename in files:
        reader = csv.reader( open( filename, 'rb' ) )
        reader.next() # skip header
        for line in reader:
            all_funcs.append(Func(line))

    # all ops in all_ops
    print
    print "== Total statistics =="
    print
    print all_funcs
    print
    print "== Function statistics =="
    print
    for func in sorted( all_funcs.getFuncs(), key=lambda entry: entry.walltime, reverse=True ):
        print "%d calls to %s; %f ms function runtime (%4.1f %%); %f ms cputime (%4.1f %%)" % \
            (func.invocations, func.name, func.walltime/1000, 100*func.walltime/all_funcs.total_walltime, \
             func.cputime, 100*func.cputime/func.walltime)

