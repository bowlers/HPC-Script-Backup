#!/usr/bin/python

import sys
import os
import functools
import pandas as pd
import numpy as np
import glob

count = 1
os.chdir( "results" )
cwd = os.getcwd()
for entry in glob.glob("*.sorted.table"):
	if not os.path.isfile( os.path.join( cwd, entry ) ):
		continue
	if count == 1 or count == 51 or count == 101:
		data = pd.read_csv( "header.lst" )
	table = pd.read_csv( entry, sep="\t" )
	if len( table['GENE'].unique() ) < len( table.index ):
		unique = table['GENE'].unique()
		names = list(table.columns)
		df = pd.DataFrame( columns = names )
		for i in unique:
			if len( table[table['GENE'] == i] ) < 1:
				continue
			tmp = table[table['GENE'] == i]
			total = tmp.sum()
			freq = len( tmp.index )
			new_row = pd.DataFrame( [[i, freq, 1, 1]], columns = names )
			df = df.append( new_row, ignore_index=True )
		df = df.sort_values( by=['GENE'] )
		table = df
	table = table.drop(["DP", "AD"], axis=1)
	data = data.merge( table, on="GENE", how="outer" )

	if count == 50:
		data.fillna(".", inplace=True)
		data.to_csv( "../Final.table1", sep="\t", quoting=None, index=False )
	elif count == 100:
		data.fillna(".", inplace=True)
		data.to_csv( "../Final.table2", sep="\t", quoting=None, index=False )
	count += 1

data.fillna(".", inplace=True)
data.to_csv( "../Final.table3", sep="\t", quoting=None, index=False )

table1 = pd.read_csv( "../Final.table1", sep="\t" )
if os.path.exists("../Final.table2"):
	table2 = pd.read_csv( "../Final.table2", sep="\t" )
table3 = pd.read_csv( "../Final.table3", sep="\t" )

if os.path.exists("../Final.table2"):
	new = pd.merge( table1, table2, on="GENE", how="outer" )
	new = pd.merge( new, table3, on="GENE", how="outer" )
elif os.path.exists("../Final.table3"):
	new = pd.merge( table1, table3, on="GENE", how="outer" )
else:
	new = table1
new.to_csv( "../done.table", sep="\t", quoting=None, index=False )

new = new.replace('.', 0)
new = new[new.astype('bool').mean(axis=1)>=0.25]

new.to_csv( "../done.75.percent.zero.removed.table", sep="\t", quoting=None, index=False )
