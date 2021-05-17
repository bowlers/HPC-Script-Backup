#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeRegressor

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
Y = pd.DataFrame(data.iloc[:,data.columns.get_loc('OS')])
data = data.drop('OS', axis =1)
X = pd.DataFrame(data.iloc[:,:])
names = data.columns

X = np.nan_to_num(X)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.25)
tree = DecisionTreeRegressor()
tree.fit(X, Y)
print('Model Accuracy: ', tree.score(X, Y))

model = RandomForestRegressor(n_estimators=20, random_state=0)
model.fit(X_train, Y_train)

# rf_predictions = model.predict(test)
# rf_probs = model.predict_proba(test)[:,1]

# from sklearn.metrics import roc_auc_score

# Calculate roc auc
# roc_value = roc_auc_score(Y_test, rf_probs)

