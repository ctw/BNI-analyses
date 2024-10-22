"""
Script to gather information about clusters analyzed for IR experiment.
"""
from glob import glob

bniBaseDir = '/Volumes/BniData/Subjects'

logFiles = sorted(glob(bniBaseDir+'/s???/data/s*e*ir/IRLog_*'))

bad_exps = ['s42e4ir', 's49e3ir', 's51e4ir', 's46e8ir']
for datfile in logFiles:
    subj = datfile.split('/')[-4]
    expstr = datfile.split('/')[-2]
    if expstr in bad_exps:
        continue
    print(f"{datfile}, {subj}, {expstr}")

    # generate list of nse files
    nsefiles = sorted(glob(bniBaseDir + subj + '/analysis/' + expstr + '/KK/CSC?.Nse'))
    nsefiles.extend(sorted(glob(bniBaseDir + subj + '/analysis/' + expstr + '/KK/CSC??.Nse')))
    nsefiles.extend(sorted(glob(bniBaseDir + subj + '/analysis/' + expstr + '/KK/CSC???.Nse')))

    # generate list of clusterfiles
    clusterfiles = sorted(glob(bniBaseDir + subj + '/analysis/' + expstr + '/KK/CSC?.clu*'))
    clusterfiles.extend(sorted(glob(bniBaseDir + subj + '/analysis/' + expstr + '/KK/CSC??.clu*')))
    clusterfiles.extend(sorted(glob(bniBaseDir+ subj + '/analysis/' + expstr + '/KK/CSC???.clu*')))

