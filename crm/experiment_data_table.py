import pandas as pd
import numpy as np
#import itertools

expdattab = pd.read_table('experimentDataTable.txt')

expdattab[~expdattab.isna()['date']]
expdattab[~expdattab.isna()['date']].duplicated('time').sum()



expdattab[expdattab['expType']=='cr']['subjectNo'].unique().shape
expdattab[expdattab['expType']=='cf']['subjectNo'].unique().shape


set.intersection(set(expdattab[expdattab['expType']=='cr']['subjectNo']),
                 set(expdattab[expdattab['expType']=='bl']['subjectNo']))
set.intersection(set(expdattab[expdattab['expType']=='sr']['subjectNo']),
                 set(expdattab[expdattab['expType']=='bl']['subjectNo']))

set.intersection(set(expdattab[expdattab['expType']=='cr']['subjectNo']),
                 set(expdattab[expdattab['expType']=='cf']['subjectNo']))
# Out[121]: {20, 21, 23, 25, 26, 27, 28, 29}

set.intersection(set(expdattab[expdattab['expType']=='sr']['subjectNo']),
                 set(expdattab[expdattab['expType']=='cf']['subjectNo']),
                 set(expdattab[expdattab['expType']=='ir']['subjectNo']))
# Out[2]: {31, 32, 33, 34, 35, 37, 40, 42, 44, 47, 48, 49}

#expdattab


exp_overlap = {}
for i, exp1 in enumerate(np.unique(expdattab['expType'])[:-1]):
    exp_overlap[exp1] = {}
    for j, exp2 in enumerate(np.unique(expdattab['expType'])[i+1:]):
        exp_overlap[exp1][exp2] = len(set.intersection(
            set(expdattab[expdattab['expType']==exp1]['subjectNo']),
            set(expdattab[expdattab['expType']==exp2]['subjectNo'])))


thresh = 10
for exp1 in exp_overlap:
    if exp1 == 'td':
        continue
    for exp2 in exp_overlap[exp1]:
        if exp2 == 'td':
            continue
        if exp_overlap[exp1][exp2] >= thresh:
            print(exp1, exp2, exp_overlap[exp1][exp2])



# expdattab_old = pd.read_table('experimentDataTable_old.txt')
# indx=5
# adder = 0
# for i, row in enumerate(expdattab.iloc[:,indx]):
#     if i==190:
#         adder = 1
#     if row != expdattab_old.iloc[i+adder,indx]:
#         #adder = 1
#         print(i, row, expdattab_old.iloc[i,indx])


            
# exp_pairs = itertools.product(np.unique(expdattab['expType']), repeat=2)

# counter = 0
# for exp1, exp2 in exp_pairs:
#     if exp1 != exp2:
#         print(exp1, exp2)
#         counter += 1

# len(set.intersection(set(expdattab[expdattab['expType']=='cr']['subjectNo']), set(expdattab[expdattab['expType']=='cf']['subjectNo'])))


# pairs = {}
# tmp = itertools.product(np.arange(4), repeat=2)
# for t1, t2 in tmp:
#     if t1==t2:
#         continue
#     if t1 in pairs:
#         if t2 not in pairs[t1]:
#             pairs[t1][t2] = t1+t2
#     elif t2 in pairs:
#         if t1 not in pairs[t2]:
#             pairs[t2][t1] = t1+t2
#     else:
#         pairs[t1] = {t2: t1+t2}
#     #print(t1,t2)



# pairs2 = {}
# for i, t1 in enumerate(np.arange(4)[:-1]):
#     pairs2[t1] = {}
#     for j, t2 in enumerate(np.arange(4)[i+1:]):
#         pairs2[t1][t2] = t1+t2


