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

#
# Calculate possible configurations for hpl on the given number of nodes and
# the given blocksize.
#
# A configuation is namely a set of P, Q and N.
#
# It is recommended for the blocksize to be in the area of 1024
#
# A few restrictions are applied:
#
# a) N % NB = 0
# b) NB % 64 = 0 (for caldgemm we currently require 64 - but 8 is the goal)
# c1) (N / NB ) % Q = 0
# c2) (N / NB) % P = 0
# d) isOdd(N / 64 / Q) (should be 8 like in (b) )
# e) isOdd(NB / 64) (should be 8 like in (b) )
#
# Not to forget about the physical limit  N <= sqrt( 0.8 * 64 Gi * P * Q / 8 )
#

# TODO sweep NB
# TODO fuzzy N (try smaller Ns, too)
# TODO project performance (only needed for fuzzy N)

import sys
import math

class Configuration:
	def __init__( self, p, q, n, nb ):
		self.p = p
		self.q = q
		self.n = n
		self.nb = nb

	def __str__( self ):
		return "P x Q = %d x %d -- N = %d -- NB = %d" % ( self.p, self.q, self.n, self.nb )

	def fullfillsRestrictions( self ):
		# a)
		if self.n % self.nb != 0:
			return False

		# b)
		if self.nb % 64 != 0:
			return False

		# c1)
		if ( self.n / self.nb ) % self.q != 0:
			return False

		# c2)
		if ( self.n / self.nb ) % self.p != 0:
			return False

		# d)
		if ( self.n / 64 / self.q ) % 2 == 0:
			return False

		# e)
		if ( self.nb / 64 ) % 2 == 0:
			return False

		return True

def splitIn2Factors( n ):
	factors = []
	for p in range( n < 4 and 1 or 2, int( math.sqrt( n ) ) + 1 ):
		if n % p == 0:
			factors.append( (p, n/p) )
	return factors

def memoryLimitForN( nodes ):
	return int( math.sqrt( .8 * 64 * 1024 * 1024 * 1024 * nodes / 8 ) )

def roundDown( n, nb ):
	return n - n % nb; # we are not using div from future, but just be sure

#
# MAIN
#
if __name__ == "__main__":

	if len( sys.argv ) != 2:
		print 'Usage: suggest-configs.py <nodes>'
		exit(1)

	nodes = int( sys.argv[1] )

	# Get all possible process configurations
	pqs = splitIn2Factors( nodes )

	print "Factors of %d" % nodes
	print pqs

	memoryLimit = memoryLimitForN( nodes )
	print "Memory Limit %d" % memoryLimit

	configs = []

	# Now check which Ns work for these configurations
	for pq in pqs:
		for nb in range( 512, 2048+1, 8 ):
			blockSpecificLimit = roundDown( memoryLimit, nb )
			for n in range( blockSpecificLimit, int( 0.8 * blockSpecificLimit ), -nb ):
				config = Configuration( pq[0], pq[1], n, nb )
				if config.fullfillsRestrictions():
					configs.append( config )

	# sort by N for output
	for config in sorted( configs, key = lambda config: ( config.n, config.p, config.nb ), reverse=True ):
		print config


