#!/usr/bin/env python

import sys
import csv

def nflop(M,N,K):
	"Calculate the number of floating point operations required to perfom the DGEMM."""
	return M * N * ( 2 * K + 2 )

def Gflops(M,N,K,time):
	"""Calculate the DGEMM performance in Gflops.

	 Time is expected in mus."""
	return (nflop(M,N,K) / time / 1000.) if time <> 0 else 0

def score(gflops,time):
	"""Calculates the score for this call

	   Calculation attempts to high score high performance calls but don't
	   punish calls that don't take a lot of time as they have a small impact
	   on overall performance."""

	return gflops / time if time <> 0 else 0

files = sys.argv[1:]

if len( files ) < 1:
	print "Please specify some input files generated using Linpack DGEMM trace code!"
	exit()

print "Scanning the following files:"
for filename in files:
	print filename
print

# Format for the files is
# 0: Row / Col Major -> 101 / 102
# 1: TransA -> 111 (non) / 112 (trans)
# 2: TransB -> 111 (non) / 112 (trans)
# 3: M
# 4: N
# 5: K
# 6: LDA
# 7: LDB
# 8: LDC
# 9: time
# added as 10: Gflops
# added as 11: score

# the histogram of all values
# key (TransA,TransB,M,N,K)
# value (score, avg. flops, avg. time)
binned = {}

for filename in files:
	for line in csv.reader( open( filename ) ):
		bin = ( line[1], line[2], int(line[3]), int(line[4]), int(line[5]) )
		time = int( line[9] )
		gflops = Gflops( bin[2],bin[3],bin[4], time )

		try:
			oldval = binned[bin]
			# in theory flops shouldbe the same ...
			oldtime = oldval[2]
			newtime = oldval[2] + time;
			newgflops = ( oldval[1] * oldtime + gflops * time ) / newtime if newtime <> 0 else 0;
			val = ( score( newgflops, newtime ), newgflops, newtime )
		except KeyError:
			val = ( score( gflops, time ), gflops, time )

		binned[bin] = val

# sort bins by score
scored = []
for (bin,val) in binned.items():
	scored.append( (val,bin) )
# sort by score (and flops)
scored = sorted( scored, key=lambda entry: entry[0][0] )

# print (one score per line)
nZeroScores = 0
lowTime = 0
printed = 0
for value in scored:
	if value[0][0] <> 0:
		if value[0][2] >= 100000:
			print value
			printed += 1
		else:
			lowTime += 1
	else:
		nZeroScores += 1

print printed,"values have been displayed"
print nZeroScores,"bins had a score of 0"
print lowTime,"bins took less than one 100 millis in total"

