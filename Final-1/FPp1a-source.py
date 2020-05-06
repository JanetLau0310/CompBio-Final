#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from collections import Counter
import math

tr_data =pd.read_table('GSE994-train.txt', sep='\s+')
te_data = pd.read_table('GSE994-test.txt', sep='\s+')
y_train = np.array(tr_data.loc['Class'])
for i in range(0,len(y_train)):
    if y_train[i] == 'NeverSmoker':
        y_train[i] = 0
    else:
        y_train[i] =1
y_train = np.array(y_train,dtype=int)

tr_idx_col = tr_data.columns.tolist()
tr_data = tr_data[:-1]

te_idx_col = te_data.columns.tolist()
te_data = te_data[:-1]

x_train = []
x_test = []

for i in range(0,len(tr_idx_col)):
    tmp = tr_data[tr_idx_col[i]].tolist()
    tmp = [ float(x) for x in tmp ]
    x_train.append(tmp)

for i in range(0,len(te_idx_col)):
    tmp = te_data[te_idx_col[i]].tolist()
    tmp = [ float(x) for x in tmp ]
    x_test.append(tmp)

def knn_(x, y, x_train, y_train, k):
    distances = [np.linalg.norm(np.array(x) - np.array(i)) for i in x_train]
    dist_label = list(zip(distances, y_train))
    sorted_label = sorted(dist_label, key=lambda x: x[0])
    k_label = [i[1] for i in sorted_label[:k]]
    counter = Counter(k_label)
    sorted_counter = sorted(counter.items(), key=lambda x: x[1], reverse=True)
    if len(sorted_counter) == 1:
        return sorted_counter[0][0], y

    if sorted_counter[0][1] == sorted_counter[1][1]:
        return knn(x, y, x_train, y_train, k - 1)
    else:
        return sorted_counter[0][0], y

def knn_list(X_te, y_te, X_tr, y_tr, k):
    y_pred = [0]*len(X_te)
    y_act = [0]*len(y_te)
    for i in range(0,len(X_te)):
        y_pred[i], y[i] = knn_(X_te[i], y_te[i], X_tr, y_tr, k)
    
    return y_pred, y_act

def split_data(X, Y, i):
    x = X[i]
    y = Y[i]
    X_tr = np.delete(X, i, 0)
    y_tr = np.delete(Y, i, 0)
    return x, y, X_tr, y_tr

acc = []
for k in range(1, 31):
    result = []
    for i in range(0,30):
        x, y, X_tr, y_tr = split_data(x_train, y_train, i)
        y_pred, y_act = knn_(x, y, X_tr, y_tr, k)
        result.append(y_pred == y_act)
        acc.append([k, sum(result) / len(result)])

#print(sorted(acc, key=lambda x: x[1], reverse=True))
#plt.figure()
#plt.plot([i[0] for i in acc], [i[1] for i in acc])
#plt.title('Accurancy in different k')

def transfer_smoker(y_test):
    res = []
    for i in range(0,len(y_test)):
        if y_test[i] == 0:
            res.append('NeverSmoker')
        else:
            res.append('CurrentSmoker')
    return res

knn=KNeighborsClassifier(n_neighbors=1,metric='euclidean')
knn.fit(x_train,y_train)
y_test = knn.predict(x_test)

y_te_1 = transfer_smoker(y_test)

knn=KNeighborsClassifier(n_neighbors=3,metric='euclidean')
knn.fit(x_train,y_train)
y_test = knn.predict(x_test)
print(y_test)
y_te_3 = transfer_smoker(y_test)

with open('FPp1-1NNoutput.txt','a') as f1:
    for i in range(0,len(y_te_1)):
        f1.write(te_idx_col[i]+'\t'+y_te_1[i]+'\n')

with open('FPp1-3NNoutput.txt','a') as f2:
    for i in range(0,len(y_te_3)):
        f2.write(te_idx_col[i]+'\t'+y_te_3[i]+'\n')


# fold 1: PATIENT1 - PATIENT6<br>
# fold 2: PATIENT7 - PATIENT12<br>
# fold 3: PATIENT13 - PATIENT18<br>
# fold 4: PATIENT19 - PATIENT24<br>
# fold 5: PATIENT25 - PATIENT30

def output_fold(test_idx, y_pred,k):
    y_pred = transfer_smoker(y_pred)
    if k == 1:
        with open('FPp1d-predictions1','a') as f1:
            for i in range(0,len(y_pred)):
                f1.write(tr_idx_col[test_idx[i]]+'\t'+y_pred[i]+'\n')
    elif k ==3:
        with open('FPp1d-predictions3','a') as f2:
            for i in range(0,len(y_pred)):
                f2.write(tr_idx_col[test_idx[i]]+'\t'+y_pred[i]+'\n')
    else:
        print("error k")

folds = KFold(n_splits=5,shuffle = False)
X = np.array(x_train)
y = np.array(y_train)
acc_1 = []
acc_3 = []

for train_index, test_index in folds.split(X,y ):
    print('test_idx = ',test_index)
    X_tr, X_te = X[train_index], X[test_index]
    y_tr, y_te = y[train_index], y[test_index]
    
    y_pred_1, y_act_1 = knn_list(X_te, y_te, X_tr, y_tr, 1)
    acc_1.append(accuracy_score(y_act_1,y_pred_1))
    output_fold(test_index,y_pred_1,1)
    
    y_pred_3, y_act_3 = knn_list(X_te, y_te, X_tr, y_tr, 3)
    acc_3.append(accuracy_score(y_act_3,y_pred_3))
    output_fold(test_index,y_pred_3,3)
