#!/usr/bin/env python
from __future__ import division
import sys
import random

def pick_samples( names, no_train, no_test ):
	"""
	This sequence is designed to pick 'no_train' training and 'no_test' test
	samples in such a way that the order in which they are assigned to either set
	is randomised.
	
	There are two main scenarios:
	(i) the sum of the sizes of training and test sets is equal to the number of samples, and
	(ii) the sum of the sizes of training and test sets is less than the number of samples
	
	Consider (i).
	We begin by picking two samples at random without replacement.
	We do not know which of the two will go to the training or test set.
	We then pick a random number in [0,1]: if it is less than 0.5 the first will go
	to the training set; otherwise the first will go to the test set PROVIDED the 
	neither are full.
	In either case we will remove the sample once it is added to either the training
	or test set.
	
	We can further break (i) into two scenarios:
	(i-i) there are fewer training than test samples,
	(i-ii) training and test sets will be of equal size, or
	(i-iii) there are more training than test samples.
	
	Suppose that the first sample picked should go to the training set.
	Consider (i-i).
	If both sets have not been filled then we will:
	- put pick1 into train
	- remove pick1 from names
	- put pick2 into test
	- remove pick2 from names
	
	If only the test set is full then we will:
	- put pick1 into train
	- remove pick1 from names
	- ignore pick2
	On the next iteration we will have another pick2.
	
	If only the training set is full then we will:
	- ignore pick1
	- put pick2 into test
	- remove pick2 from names
	On the next iteration we will have another pick1.
	
	If both sets are full then will:
	- break from the while loop
	
	The exact same scenario but trivially reversed plays out when the first sample
	picked should go to the test set.
	
	It is easy to see that (i-ii) and (i-iii) are already present.
	It is also easy to see that (ii) and its cases will also be represented.
	"""
	train = list()
	test = list()

	while len( train ) < no_train or len( test ) < no_test:
		pick1 = random.choice( names )
		if len( names ) > 1:
			names.remove( pick1 )
		pick2 = random.choice( names )
		names.append( pick1 )
		
		train_first = True
		if random.random() < 0.5:
			train_first = False
		
		if len( train ) < no_train:
			if train_first:
				train.append( pick1 )
				names.remove( pick1 )
			else:
				train.append( pick2 )
				names.remove( pick2 )
		
		if len( test ) < no_test:
			if train_first:
				test.append( pick2 )
				names.remove( pick2 )
			else:
				test.append( pick1 )
				names.remove( pick1 )
		
	return train, test 
	
def make_row( name, names, affy_ids, train=True, test_present=True ):
	if train and test_present:
		print "%s\t%s" % ( name, affy_ids[ names[ name ]] )
	elif not train and test_present:
		print "*%s\t%s" % ( name, affy_ids[ names[ name ]] )
	elif not train and not test_present:
		print "*NA\t%s" % ( affy_ids[ names[ name ]] )
	else:
		raise ValueError( "Hmmm... I'm not designed for this. What should I do, my lord?" )
		return
	
def check_for_duplicates( train, test ):
	
	intsct = len( set( train ).intersection( set( test ) ))
	
	if intsct == 0:
		print >> sys.stderr, "Good: Correctly formed samples file; no duplicates found"
	else:
		print >> sys.stderr, "Warning: Corrupt samples file; duplicates found"
	
	return
	
if __name__ == "__main__":
	try:
		tissue_fn = sys.argv[1]
		affy_ids_fn = sys.argv[2]
		tissue_rnaseq_fn = sys.argv[3]
		no_train = int( sys.argv[4] )
		no_test = int( sys.argv[5] )
	except IndexError:
		print >> sys.stderr, """\
Script to create a tissue-specific 'samples.txt' for MaLTE.
Selects a random subset of samples as training and rest as test.

usage: ./script.py <tissue_fn> <affy_ids_fn> <tissue_rnaseq_gene_fn> <no_train> <no_test>

example:
./script.py mapped_sample_names_by_tissue/Whole_Blood.txt Affymetrix_sample_names.txt Whole_Blood/GTEx_gene_rpkm.txt 20 30 > Whole_Blood/samples.txt"""
		sys.exit( 0 )

	f = open( tissue_fn ) # tissue sample names
	names = dict()
	for row in f:
		L = row.strip().split( '\t' )
		names[ L[0] ] = L[1]
	f.close()

	f = open( affy_ids_fn ) # affymetrix mapping names
	affy_ids = dict()
	for row in f:
		L = row.strip().split( '\t' )
		affy_ids[ L[0] ] = L[1]
	f.close()

	f = open( tissue_rnaseq_fn ) # GTEx rnaseq (hts) names
	gtex_rnaseq_names = f.readline().strip().split( '\t' )[1:]
	f.close()
	
	try:
		assert len( gtex_rnaseq_names ) >= no_train + no_test
	except:
		raise ValueError( "Required number of training and test samples exceeds available samples. No. of RNA-Seq samples: %s" % len( gtex_rnaseq_names ) )
		sys.exit( 1 )

	train, test = pick_samples( gtex_rnaseq_names, no_train, no_test )
	
	check_for_duplicates( train, test )

	print "hts\tma"
	for s in train:
		make_row( s, names, affy_ids )
	for s in test:
		make_row( s, names, affy_ids, train=False )
