"""
Script to summarize clusters in all experiments.
"""

from exp_summary import expSummary

bniSubjectsDir = '/Volumes/BniData/Subjects/'

df = expSummary(bniSubjectsDir, 'IRLog_*', ['s42e4ir', 's49e3ir', 's51e4ir', 's46e8ir'])

df.to_csv('analyzedClusters.txt', sep='\t', index=False)