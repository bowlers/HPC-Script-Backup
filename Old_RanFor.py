#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

if len(sys.argv) == 1:
	sys.exit("No study number provided")
else:
	ARG1 = sys.argv[1]

DIR = "/home/sbowler/biocore_lts/scott/"
DIR += ARG1

if not os.path.exists(Path(DIR)):
	sys.exit("Could not locate study directory")
else:
	os.chdir(Path(DIR))

FILE = DIR
FILE += "/matrix.tsv"

if not os.path.exists(Path(FILE)):
	sys.exit("Matrix file does not exist within study directory")

tsv = pd.read_csv(Path(FILE), sep='\t', encoding='ISO-8859-1')

#####################################################
# Files are loaded.  Now to setup the Random Forest #
#####################################################
# Labels are values we want to predict (in this case: Overall Survival [OS])
labels = np.array(tsv['OS'])

# Remove the labels from the tsv data structure
# Axis 1 refers to the columns
tsv = tsv.drop('OS', axis = 1)

# Save tsv names for later use
tsv_list = list(tsv.columns)

# Convert to numpy array
tsv = np.array(tsv)

# Split the data into training and testing sets
# Using a random_state = 42 to allow for reproducible results 
train_tsv, test_tsv, train_labels, test_labels = train_test_split(tsv, labels, test_size = 0.25, random_state = 42 )

print('Training Features Shape:', train_tsv.shape)
print('Training Labels Shape:', train_labels.shape)
print('Testing Features Shape:', test_tsv.shape)
print('Testing Labels Shape:', test_labels.shape)

# The baseline predictions are the historical averages
baseline_preds = 813

# Baseline errors, and display average baseline error
baseline_errors = abs(baseline_preds)

print('Average baseline error: ', round(np.mean(baseline_errors), 2))

# Instantiate model with 1000 decision trees
rf = RandomForestRegressor(n_estimators = 1000, random_state = 42)

# Train the model on training data
rf.fit(train_tsv, train_labels);

# User the forest's predict method on the test data
predictions = rf.predict(test_tsv)

# Calculate the absolute errors
errors = abs(predictions - test_labels)

# Print out the mean absolute error (mae)
print('Mean Absolute Error:', round(np.mean(errors), 2), 'days.')

# Calculate mean absolute percentage error
mape = 100 * (errors / test_labels)

# Calculate and display accuracy
accuracy = 100 - np.mean(mape)
print('Accuracy: ', round(accuracy, 20), '%.')

from sklearn.tree import export_graphviz
import pydot

tree = rf.estimators_[5]
export_graphviz(tree, out_file = 'tree.dot', feature_names = tsv_list, rounded = True, precision = 1)
(graph, ) = pydot.graph_from_dot_file('tree.dot')
graph.write_png('tree.png')
