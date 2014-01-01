#!/usr/bin/env python

#	combine_microarrays.py

#	Transform a Affymetrix Human Exon 1.0 ST array or Human Gene 1.0 array
# into a Human Gene 1.1 array for application to MaLTE using GTEx training data
# It requires several resource files to run. Please 

#	Copyright (C) 2013  Paul K. Korir

#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.

#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.

#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import sys
import argparse
import random

comment = [ "#" ] 		# comments
header = [ "p", "A" ]	# header markers

class ProbeIntensity( object ):
	"""
	probe objects created from each row of an intensity file
	"""
	def __init__( self, probe_id, x, y, probe_type, probeset_id, probeset_type, block, intensities ):
		self.probe_id = probe_id
		self.x = x
		self.y = y
		self.probe_type = probe_type
		self.probeset_id = probeset_id
		self.probeset_type = probeset_type
		self.block = block
		self.intensities = intensities
	
	def __repr__( self ):
		probe_intensity_str = "\t".join([ self.probe_id, self.x, self.y, self.probe_type, self.probeset_id, self.probeset_type, self.block ] + self.intensities )
		return probe_intensity_str
	
	def print_as( self, probe_id, probeset_id ):
		probe_intensity_str = "\t".join([ probe_id, self.x, self.y, self.probe_type, probeset_id, self.probeset_type, self.block ] + self.intensities )
		return probe_intensity_str

class Probe( object ):
	def __init__( self, probe_id, probe_type, gc_count, probe_length, interrogation_position, probe_sequence ):
		self.probe_id = probe_id
		self.probe_type = probe_type
		self.gc_count = gc_count
		self.probe_length = probe_length
		self.interrogation_position = interrogation_position
		self.probe_sequence = probe_sequence
	
	def __eq__( self, Probe ):
		if self.probe_sequence == Probe.probe_sequence:
			return True
		else:
			return False
		
	def __repr__( self ):
		return self.probe_id #+ "(" + self.probe_sequence + ")"
	
class Probeset( object ):
	"""
	probe set objects created from the probe set details file
	"""
	def __init__( self, probeset_id ):
		self.probeset_id = probeset_id
		self.probes = list()
	
	def add_probe( self, probe ):
		self.probes.append( probe )
	
	def __repr__( self ):
		return "Probeset ID %s has %s probes (%s)" % ( self.probeset_id, len( self.probes ), self.probes )
	
class ProbeIntensityFile( object ):
	def __init__( self, name ):
		self.name = name
		self.header = None		
		self.probes = dict()
		
		self.parse()
	
	def parse( self ):
		print >> sys.stderr, "Reading intensities from %s ..." % self.name,
		with open( self.name ) as f:
			for row in f:
				if row[0] in comment: continue
				if row[0] in header: self.header = row.strip()
				L = row.strip().split( "\t" )
				self.probes[ L[0] ] = ProbeIntensity( L[0], L[1], L[2], L[3], L[4], L[5], L[6], L[7:] )
		print >> sys.stderr, "OK"
	
	def __getitem__( self, probe_id ):
		try:
			probes = self.probes[ probe_id ]
		except KeyError:
			probes = None
		return probes
	
class ProbesetDetailFile( object ):
	def __init__( self, name ):
		self.name = name
		self.probesets = dict()
		
		self.parse()
	
	def parse( self ):
		print >> sys.stderr, "Reading probe details from %s ..." % self.name,
		with open( self.name ) as f:
			for row in f:
				if row[0] in comment: continue
				if row[0] in header: continue
				L = row.strip().split( "\t" )
				probe = Probe( L[1], L[2], L[3], L[4], L[5], L[6] )
				probeset = Probeset( L[0] )
				if L[0] not in self.probesets:
					probeset.add_probe( probe )
					self.probesets[ L[0] ] = probeset
				else:
					self.probesets[ L[0] ].add_probe( probe )
		print >> sys.stderr, "OK"
	
	def __getitem__( self, probeset_id ):
		try:
			probes = self.probesets[ probeset_id ].probes
		except KeyError:
			probes = None
		return probes

if __name__ == "__main__":
	parser = argparse.ArgumentParser( description="Script to combine two different microarray types into one" )
	
	parser.add_argument( "-e", "--huex-intensities", help="probe intensities from Affymetrix HuEx_1-0_st arrays" )
	parser.add_argument( "-g", "--huge-intensities", help="probe intensities from Affymetrix HuGe_1-0_st or closely related arrays" )
	parser.add_argument( "-c", "--comparison-map", help="file containing a map between array types" )
	parser.add_argument( "-f", "--huex-details", help="output from 'apt-dump-pgf' with probe set details for Affymetrix HuEx_1-0_st arrays" )
	parser.add_argument( "-i", "--huge-details", help="output from 'apt-dump-pgf' with probe set details for Affymetrix HuGe_1-0_st arrays" )
	parser.add_argument( "-p", "--huex-mps-map", help="metaprobeset-to-probeset map for Affymetrix HuEx_1-0_st arrays" )
	parser.add_argument( "-q", "--huge-mps-map", help="metaprobeset-to-probeset map for Affymetrix HuGe_1-0_st arrays" )
	parser.add_argument( "--huex-out", help="output file for new Affymetrix HuEx_1-0_st probes" )
	parser.add_argument( "--huge-out", help="output file for new Affymetrix HuGe_1-0_st probes" )
	
	args = parser.parse_args()
	
	huex_fn = args.huex_intensities
	huge_fn = args.huge_intensities
	cmp_fn = args.comparison_map
	huexde_fn = args.huex_details
	hugede_fn = args.huge_details
	huex_mps_fn = args.huex_mps_map
	huge_mps_fn = args.huge_mps_map
	huex_out_fn = args.huex_out
	huge_out_fn = args.huge_out
	
	# get the map between probe sets
	print >> sys.stderr, "Building map of HuGe metaprobesets to HuEx metaprobesets ...",
	metaprobeset_map = dict()
	with open( cmp_fn ) as f:
		for row in f:
			if row[0] in header: continue
			L = row.strip().split( "\t" )
			if L[5] not in metaprobeset_map:
				metaprobeset_map[ L[5] ] = [ L[2] ] # huge -> huex
			else:
				metaprobeset_map[ L[5] ] += [ L[2] ]
	print >> sys.stderr, "OK"
	
	# get the map between mps and ps
	print >> sys.stderr, "Building map of HuEx metaprobesets to probesets ...",
	huex_metaprobesets = dict()
	with open( huex_mps_fn ) as f:
		for row in f:
			if row[0] in comment or row[0] in header: continue
			L = row.strip().split( "\t" )
			huex_metaprobesets[ L[0] ] = L[2].split( " " )
	print >> sys.stderr, "OK"
	
	print >> sys.stderr, "Building map of HuGe metaprobesets to probesets ...",
	huge_metaprobesets = dict()
	with open( huge_mps_fn ) as f:
		for row in f:
			if row[0] in comment or row[0] in header: continue
			L = row.strip().split( "\t" )
			huge_metaprobesets[ L[0] ] = L[2].split( " " )
	print >> sys.stderr, "OK"
	
	# probeset map
	print >> sys.stderr, "Building map of HuGe probesets to HuEx probesets ...",
	probeset_map = dict()
	for huge_mps in metaprobeset_map:
		# if one huge_mps is associated with more than one pick one huex_mps at random
		if len( metaprobeset_map[ huge_mps ] ) > 1:
			huex_mps = random.choice( metaprobeset_map[ huge_mps ] )
		else:
			huex_mps = metaprobeset_map[ huge_mps ][0]
	
		# huge: see if there are probesets
		try:
			huge_probesets = huge_metaprobesets[ huge_mps ]
		except KeyError:
			huge_probesets = None
		
		if huge_probesets == None: continue
		
		# huex: ditto
		try:
			huex_probesets = huex_metaprobesets[ huex_mps ]
		except KeyError:
			huex_probesets = None
			
		if huex_probesets == None: continue
		
		for huge_ps in huge_probesets:
			probeset_map[ huge_ps ] = huex_probesets
	print >> sys.stderr, "OK"
	
	huex_int = ProbeIntensityFile( huex_fn )
	huge_int = ProbeIntensityFile( huge_fn )
	
	huex_det = ProbesetDetailFile( huexde_fn )
	huge_det = ProbesetDetailFile( hugede_fn )
	
	of1 = open( huge_out_fn, 'w' )
	of2 = open( huex_out_fn, 'w' )
	of3 = open( "huex_probes.txt", 'w' )
	
	print >> of1, huge_int.header
	print >> of2, huex_int.header
	print >> of3, "probe_id"
	
	print >> sys.stderr, "Writing new probe intensity files ...",
	c = 0
	# get the ps
	for huge_id in probeset_map:
		if c > 10: break
		# get the corresponding probe set ID in huex
		for huex_id in probeset_map[ huge_id ]:
			# get the probes in each from the detail file
			huex_probes = huex_det[ huex_id ]	# have to create __getitem__() method
			huge_probes = huge_det[ huge_id ]
			
			# get exactly matching probes
			if huex_probes == None or huge_probes == None:
				continue
			
			# compare probes
			similar_probes = list()
			for huge_probe in huge_probes:
				for huex_probe in huex_probes:
					if huge_probe == huex_probe: # have to create __eq__() method
						similar_probes += [ ( huge_probe.probe_id, huex_probe.probe_id ) ]
		
			# if there are no probes in common then move on...		
			if len( similar_probes ) == 0: continue

			# now we have the set of probes that are similar between platforms
			for huge_probe, huex_probe in similar_probes:
				# write the intensity of this huge probe to an output file
				print >> of1, huge_int[ huge_probe ] # __getitem__()
			
				# write the intensity of this huex probe to an output file
				# use the huge_probe and huge_id (huge probe set ID)
				print >> of2, huex_int[ huex_probe ].print_as( probe_id=huge_probe, probeset_id=huge_id )
				print >> of3, huex_probe
		c += 0
	
	of1.close()
	of2.close()
	of3.close()
	
	print >> sys.stderr, "OK"
