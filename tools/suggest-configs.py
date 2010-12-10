#!/usr/bin/env python
# vim: set fileencoding=utf-8 :

#
# Calculate possible configurations for hpl on the given number of nodes.
#
# A configuation is namely a set of P, Q, N and NB.
#
# A few restrictions are applied:
#
# a) N % NB = 0
# c) (N / NB ) % Q = 0
# c2) (N / NB) % P = 0
#
# NB = 1024
#
# Not to forget about the physical limit  N <= sqrt( 0.8 * 64 Gi * P * Q / 8 )
#
#
#
# Copyright 2010:
#  - David Rohr (drohr@jwdt.org)
#  - Matthias Bach (bach@compeng.uni-frankfurt.de)
#  - Matthias Kretz (kretz@compeng.uni-frankfurt.de)
#
# This file is part of HPL-GPU.
#
# HPL-GPU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HPL-GPU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HPL-GPU.  If not, see <http://www.gnu.org/licenses/>.
#
# In addition to the rules layed out by the GNU General Public License
# the following exception is granted:
#
# Use with the Original BSD License.
#
# Notwithstanding any other provision of the GNU General Public License
# Version 3, you have permission to link or combine any covered work with
# a work licensed under the 4-clause BSD license into a single combined
# work, and to convey the resulting work.  The terms of this License will
# continue to apply to the part which is the covered work, but the special
# requirements of the 4-clause BSD license, clause 3, concerning the
# requirement of acknowledgement in advertising materials will apply to
# the combination as such.
#

# TODO project performance (only needed for fuzzy N)
# TODO prefer larger matrices more
# TODO prefer squarer (so not fully square) pxq

import sys
import math
import optparse

class Configuration:
	def __init__( self, p, q, n, nb ):
		self.p = p
		self.q = q
		self.n = n
		self.nb = nb

	def __str__( self ):
		return "Nodes = %d -- P x Q = %d x %d -- N = %d -- NB = %d -- Mem. used per node = %.2f GiB" % \
			( self.p*self.q, self.p, self.q, self.n, self.nb, self.memoryPerNode() )

	def fullfillsRestrictions( self ):
		# a)
		if self.n % self.nb != 0:
			return False

		# c1)
		if ( self.n / self.nb ) % self.q != 0:
			return False

		# c2)
		if ( self.n / self.nb ) % self.p != 0:
			return False

		return True

	def memoryPerNode( self ):
		""" Calculate the memory used per node in GiB """
		return float( self.n**2 ) * 8 / self.p / self.q / 1024**3

	def eta( self, gflopsPerNode ):
		""" Calculate estimated computation time in seconds """
		return ( float(2) / 3 * self.n**3 - float( self.n**2 ) / 2 ) / ( gflopsPerNode * self.p * self.q * 10**9 )

	def memoryPerNode( self ):
		""" Calculate the memory used per node in GiB """
		return float( self.n**2 ) * 8 / self.p / self.q / 1024**3

def splitIn2Factors( n, flat=False ):
	factors = []
	min_p = flat and 2 or int( math.sqrt( n / 2 ) )
	for p in range( n < 4 and 1 or min_p, int( math.sqrt( n ) ) + 1 ):
		if n % p == 0:
			factors.append( (p, n/p) )
	return factors

def memoryLimitForN( nodes, memPerNode ):
	return int( math.sqrt( memPerNode * 1024 * 1024 * 1024 * nodes / 8 ) )

def roundDown( n, nb ):
	return n - n % nb; # we are not using div from future, but just be sure

#
# MAIN
#
if __name__ == "__main__":

	parser = optparse.OptionParser( usage='Usage: %prog [options] nodes', description='Suggests values for PxQ and N for a given number of nodes.' )
	parser.add_option( '-m', '--memory', metavar='memory', default=.8*64, type=float, help='Memory to use per node in GiB' )
	parser.add_option( '-p', '--performance', metavar='performance', default=500, type=int, help='Performance of one node in Gflops' )
	parser.add_option( '-r', '--relaxed', metavar='relaxed', action='store_true', default=False, help='Do not apply tiling restricitons, simply give maximum matrix size and posible process configurations' )
	parser.add_option( '-f', '--flat', metavar='relaxed', action='store_true', default=False, help='Enable flat configurations' )


#	parser.add_option('dirs', metavar='dir', nargs='+', help='The direcotories to scan for HPL*out files' )
	( args, tmp ) = parser.parse_args()

	if len( tmp ) != 1:
		parser.error( 'Need to give a number of nodes' )
		exit -1

	nodes = int( tmp[0] )
	nb = 1024

	memPerNode = args.memory
	perfPerNode = args.performance

	configs = []

	while len( configs ) < 3 and nodes > 0:
		# Get all possible process configurations
		pqs = splitIn2Factors( nodes, args.flat )

		memoryLimit = memoryLimitForN( nodes, memPerNode )
		blockSpecificLimit = roundDown( memoryLimit, nb )

		# Now check which Ns work for these configurations
		for pq in pqs:
			if args.relaxed:
				configs.append( Configuration( pq[0], pq[1], blockSpecificLimit, nb ) )
			else:
				for n in range( blockSpecificLimit, int( 0.8 * blockSpecificLimit ), -nb ):
					config = Configuration( pq[0], pq[1], n, nb )
					if config.fullfillsRestrictions():
						configs.append( config )
						break

		# in case we did not find at least three configs till now we will try
		# again with a smaller number of nodes
		nodes -= 1

	# sort by N for output
	# We prefer, in that order
	#  1) Large matrices
	#  2) With many process rows (more quadratic)
	#  3) That waste little GPU buffer
	#  4) And have large blocks
	for config in sorted( configs, key = lambda config: ( config.p*config.q, config.n, config.p ), reverse=True ):
		eta = int( config.eta( perfPerNode ) );
		etaHour = eta / 3600;
		etaMin = eta % 3600 / 60;
		print "%s -- ETA = %d:%02d" % \
			( config, etaHour, etaMin )

