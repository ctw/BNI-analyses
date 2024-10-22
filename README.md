# BNI-analyses
Patterns of distributed neural firing code for perceptual features analyses files.

## Dataflow

There are 8 experiments analyzed for the paper: cf, crm, ir, oi, ReberEtAl2019, rp, sc, sr. The 7 
other than ReberEtAl2019 were all recorded at the Barrow Neurological Institute in Phoenix. 

### cf - Conflict Flankers

Events, or trials, in the experiment are extracted from the CFLog_*.txt files in the corresponding 
BniData {subject}/data/{expID} directories and then stored together in the .data/bni_evs_cf.npz 
file. This is done in the cf/create_events_cf.py file and includes all such CFLog_*.txt files in 
experiment sessions.

Spike counts may never have been computed since this experiment was not included in further 
analyses.

### crm - Continuous Recognition Memory

Events in the experiment are extracted from the CrmLog_*.txt files in the corresponding 
BniData/{subject}/data/{expID} directories. Apparently they are stored in a 
'/home/ctw/data/BNI/data/bni_evs_crm.npz/bni_evs_crm.npz' file though it is not clear from the 
extant crm/create_events_crm.py.old file. 

Spike counts are calculated and used in k-means fitting in both spike_pattern.py and 
spike_pattern_firstsecond.py. This occurs for all the .Nse files for all experiments for which 
the events have been isolated. The names of clusters are then combined with 
the brain areas in a 'cluster' numpy array but it is unclear where this is stored or output. 

### ir - invariant recognition

Events in the experiment are extracted from the IRLog_*.txt files in the corresponding
BniData/{subject}/data/{expID} directories by the ir/create_events_ir.py script. The 
results is then stored in the ./data/bni_evs_ir.npz file.

Spike counts are calculated from all .Nse files which have a corresponding .clu.1 file
available by ir/spike_pattern_reps_and_cats.py.

Cluster information is extracted from the corresponding clusterInfo.txt files and only 
clusters which are POTENTIAL or SPIKE are used.

### overall analyses and graphs

all_exps.py iterates over 7 experiments, excluding cf, and gathers z_mean scores by experiment and 
brain area. Each of these is stored in a corresponding all_ps_zmean.pickle file. 