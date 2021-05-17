#!/usr/bin/python

import sys
import os
import functools
import pandas as pd
import numpy as np
import getopt

def splitext_(path):
    if len(path.split('.')) > 2:
        return path.split('.')[0],'.'.join(path.split('.')[-2:])
    return splitext(path)

geneFile=''
output=''
study=''
try:
	opts, args = getopt.getopt(sys.argv[1:],"hg:s:o:",["gfile=","sdir=","ofile="])
except getopt.GetoptError:
	print( 'LOF.py -g <gene list> -s <study> -o <output file>' )
	sys.exit(2)

if len(sys.argv) == 1:
	print( 'LOF.py -g <gene list> -s <study> -o <output file>' )
	sys.exit(2)

for opt, arg in opts:
	if opt == '-h':
		print( 'LOF.py -g <gene list> -s <study> -o <output file>' )
		sys.exit(2)
	elif opt in ("-g", "--gfile"):
		geneFile = arg
	elif opt in ("-s", "--sdir"):
		study = arg
	elif opt in ("-o", "--ofile"):
		output = arg
	else:
		print( 'LOF.py -g <gene list> -s <study> -o <output file>' )
		sys.exit()

if output == '':
	output = 'final.tsv'

os.chdir("/home/sbowler/biocore_lts/scott/")
isPath = os.path.exists( study )
if not isPath:
	print( 'Could not locate study directory' )
	sys.exit(2)
else:
	os.chdir( study )
	os.chdir( "varscan/lof/filtered/lof" )
	checkGfile = os.path.isfile( geneFile )
	if not checkGfile:
		print( 'Could not locate gene file.' )
		sys.exit(2)
	tName =[ 'GeneID' ]
	table = pd.read_csv( geneFile, sep="|", header=None, names=tName )
	table.to_csv( 'table.txt', sep='\t', quoting=None, index=False )
	cwd = os.getcwd()
	data = ''
	for entry in os.listdir( cwd ):
		if not os.path.isfile( os.path.join( cwd, entry ) ):
			continue
		elif os.path.splitext( entry )[1] != ".txt":
			continue
		NAME = os.path.splitext( entry )[0]
		fName =[ 'GeneID', '' ]
		fName[1]  = NAME
		data = pd.read_csv( entry, header=None, sep="\t", names=fName )
		table = pd.read_csv( 'table.txt', sep="\t" )
		table_nodups = table.drop_duplicates()
		result = pd.merge(left=table_nodups, right=data, how='left', left_on='GeneID', right_on='GeneID')
		result.to_csv( 'table.txt', sep='\t', quoting=None, index=False )
