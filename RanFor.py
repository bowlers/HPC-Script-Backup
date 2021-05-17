#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score

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

data = pd.read_csv(Path(FILE), sep='\t', encoding='ISO-8859-1')

#####################################################
# Files are loaded.  Now setup for training/test.   #
#####################################################
data = data.drop('PID', axis =1)
data = data.drop('Study', axis =1)
data = data.drop('Type', axis =1)
y = data['OS']
data = data.drop('OS', axis =1)
featuers = data.columns

X = np.nan_to_num(X)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25)
rfc = RandomForestClassifier(n_estimators=1000, random_state=42, max_depth=None)
rfc.fit(X_train, Y_train)
Y_pred = regressor.predict(X_test)

print("Features sorted by their score: ")
print(sorted(zip(map(lambda x: round(x, 4), rfc.feature_importances_), names), reverse=True))


