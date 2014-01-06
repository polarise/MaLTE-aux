#!/usr/bin/env python

#	create_samples_template.py
#
#	Create a template 'samples.txt' file using the pairs of training-
#	test samples provided in the GTEx dataset using the file
#	'training_samples.txt.gz' as input. See '--help' for details.
#
#	Copyright (C) 2014 Paul K. Korir
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import division, print_function
import sys
import argparse
import random

def pick_random_samples( list_samples, no_samples ):
	picked_samples = list()
	while len( picked_samples ) < no_samples:
		picked_samples.append( random.choice( list_samples ))
	
	return picked_samples

if __name__ == "__main__":
	parser = argparse.ArgumentParser( "Script to create a template 'samples.txt' file that optionally adds the test samples." )

	parser.add_argument( '-t', "--train", default="training_samples.txt", help="the file containing a map of training sample names between RNA-Seq and microarray samples [default: 'training_samples.txt' (provided)]" )
	parser.add_argument( '-e', "--test", default=None, help="the file containing a map of test sample names; expected to only contain microarray samples with the first column blank (e.g. '\tTEST_SAMPLE.CEL' as one such test sample row)" )
	parser.add_argument( '-a', "--no-train", type=int, default=150, help="the number of training samples to be used; should not exceed the number of training samples [default: 150]" )
	parser.add_argument( '-b', "--no-test", type=int, default=0, help="the number of test samples; if no test samples are specified then a samples template file will be generated [default: None]" )
	parser.add_argument( '-o', "--outfile", default="samples.txt", help="the name of the output file [default: 'samples.txt']" )

	args = parser.parse_args()
	train_fn = args.train
	test_fn = args.test
	no_train = args.no_train
	no_test = args.no_test
	outfile = args.outfile

	with open( train_fn ) as f:
		all_train_samples = [ pair.strip() for pair in f ]

	#print( all_train_samples )
	
	if test_fn != None:
		with open( test_fn ) as f:
			all_test_samples = [ pair.strip() for pair in f ]
	
	# simple validation
	# validate the number of training/test samples
	try:
		assert no_train <= len( all_train_samples )
		if test_fn != None:
			assert no_test <= len( all_test_samples )
	except:
		raise ValueError( "Invalid number of training/test samples." )

	# other necessary validation

	# get training samples
	training_samples = pick_random_samples( all_train_samples, no_train )

	# get test samples (if required)	
	if test_fn != None:
		test_samples = pick_random_samples( all_test_samples, no_test )

	# create the 'samples.txt' file
	with open( outfile, 'w' ) as of:
		# print the header
		print( "hts", "ma", sep="\t", file=of )
	
		# print the training samples
		for pair in training_samples:
			print( pair, file=of )

		# print the test samples (if any)
		if test_fn != None:
			for pair in test_samples:
				print( pair, file=of )
