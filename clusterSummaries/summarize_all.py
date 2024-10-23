"""
Script to summarize clusters in all experiments.
"""
import logging

import pandas as pd
from exp_summary import expSummary
from ReberEtAl2019.exp_summary import expSummary as reberExpSummary

logging.basicConfig(level=logging.INFO)

bniSubjectsDir = '/Volumes/BniData/Subjects/'

# crm
# exclude s13e15cr per crm/spike_pattern.py because CSC5 sorting output missing
df = expSummary(bniSubjectsDir, "s*e*cr*/CrmLog_*", ['s23e15cr'])
df['expType'] = 'cr'
allDf = df

# ir
# exclude s42e4ir and s49e3ir per ir/spike_pattern_reps_and_cats.py
# exclude s51e4ir and s46e8ir per ir/create_events_ir.py
df = expSummary(bniSubjectsDir, "s*e*ir*/IRLog_*", ['s42e4ir', 's49e3ir', 's51e4ir', 's46e8ir'])
df['expType'] = 'ir'
allDf = pd.concat([allDf, df])

# oi
df = expSummary(bniSubjectsDir, 's*e*oi/OILog_*', [])
df['expType'] = 'oi'
allDf = pd.concat([allDf, df])

# rp
# exclude s14e2rp per rp/create_events_rp.py
# exclude s28e15rp and s14e26rp per rp/spike_pattern_reps_and_cats.py
df = expSummary(bniSubjectsDir, 's*e*rp/RPLog_*', ['s14e2rp', 's28e15rp', 's14e26rp'])
df['expType'] = 'rp'
allDf = pd.concat([allDf, df])

# sc
df = expSummary(bniSubjectsDir, 's*e*sc/SplitCategoryLog_*', [])
df['expType'] = 'sc'
allDf = pd.concat([allDf, df])

# sr
df = expSummary(bniSubjectsDir, 's*e*sr/SRLog_*', [])
df['expType'] = 'sr'
allDf = pd.concat([allDf, df])

# ReberEtAl2019
df = reberExpSummary('/Users/Shared/Mormann/abstractRepresentationsInMTL/zvals_trials.h5')
df['expType'] = 'rb'
allDf = pd.concat([allDf, df])

allDf.to_csv('../results/analyzedClusters.txt', sep='\t', index=False)