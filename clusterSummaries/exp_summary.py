"""
Summarize information for clusters in one experiment type.
"""
from glob import glob
import logging
import os

import numpy as np
from numpy.lib import recfunctions
import pandas as pd

def expSummary(bniSubjectsDir, logFilePattern, bad_exps):
    """
    :param bniSubjectDir: BniData/Subjects directory containing all subject subdirectories
    :param logFilePattern: filename pattern to match log files for experiment
    :param bad_exps: list of expIds of experiments to exclude

    :return pd.DataFrame with information for each cluster
    """

    logFiles = sorted(glob(bniSubjectsDir + '/s???/data/' + logFilePattern))

    # Exclude original _fixed log files as per sr/create_events_sr.py. In these data sets
    # that only matters for sr experiments.
    for lf in logFiles:
        if lf[:-4] + '_fixed.txt' in logFiles:
            logging.warning(f"excluding {lf} as fixed copy exists")
            logFiles.remove(lf)

    df = pd.DataFrame()

    for datfile in logFiles:
        subj = datfile.split('/')[-4]
        expstr = datfile.split('/')[-2]
        if expstr in bad_exps:
            continue
        logging.info(f"{datfile}, {subj}, {expstr}")

        # generate list of nse files
        nsefiles = sorted(glob(bniSubjectsDir + subj + '/analysis/' + expstr + '/KK/CSC?.Nse'))
        nsefiles.extend(sorted(glob(bniSubjectsDir + subj + '/analysis/' + expstr + '/KK/CSC??.Nse')))
        nsefiles.extend(sorted(glob(bniSubjectsDir + subj + '/analysis/' + expstr + '/KK/CSC???.Nse')))

        # generate list of clusterfiles
        clusterfiles = sorted(glob(bniSubjectsDir + subj + '/analysis/' + expstr + '/KK/CSC?.clu*'))
        clusterfiles.extend(sorted(glob(bniSubjectsDir + subj + '/analysis/' + expstr + '/KK/CSC??.clu*')))
        clusterfiles.extend(sorted(glob(bniSubjectsDir + subj + '/analysis/' + expstr + '/KK/CSC???.clu*')))

        # compare nse and .clu.1 files channel numbers
        assert len(nsefiles) == len(clusterfiles), \
            'len(nsefiles) == len(clusterfiles)'
        assert np.all([int(clusterfile.split('CSC')[-1].split('.')[0]) ==
                       int(nsefile.split('CSC')[-1].split('.')[0])
                       for clusterfile, nsefile in
                       zip(clusterfiles, nsefiles)
                       if clusterfile is not None]), 'mismatch in nse & cluster files'

        # get cluster quality, area, etc.
        clusterinfo_file = (bniSubjectsDir + subj + '/analysis/' + expstr + '/clusterInfo.txt')
        if os.path.exists(clusterinfo_file):
            try:
                clusterinfo = np.loadtxt(
                    clusterinfo_file, delimiter='\t', skiprows=1,
                    dtype={'names': ('expstr', 'cluster_id',
                                     'hemisphere', 'area', 'quality'),
                           'formats': ('U16', 'U16', 'U8', 'U8', 'U16')})
            except ValueError:
                clusterinfo = np.loadtxt(
                    clusterinfo_file, delimiter='\t', skiprows=1,
                    dtype={'names': ('cluster_id', 'hemisphere',
                                     'area', 'quality'),
                           'formats': ('U16', 'U8', 'U8', 'U16')})
                clusterinfo = recfunctions.append_fields(
                    clusterinfo, 'expstr', [expstr] * len(clusterinfo),
                    usemask=False, asrecarray=True)
            if len(clusterinfo) > 0:
                assert len(np.unique(clusterinfo['expstr'])) == 1, \
                    'clusterinfo expstr 1'
                assert clusterinfo['expstr'][0] == expstr, 'clusterinfo expstr 2'
        else:
            logging.warning(f'Skipping -- no clusterinfo file {clusterinfo_file}')
            continue

        # consolidate information for each cluster that was analyzed
        for clusterfile, nsefile in zip(
                clusterfiles, nsefiles):
            if clusterfile is None:
                continue
            chan = nsefile.split('.Nse')[0]
            channum = int(chan.split('CSC')[-1])

            sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]

            for cluster_num in np.unique(sess_clusts):
                ci = clusterinfo[clusterinfo[
                                     'cluster_id'] == ('ch' + str(channum) + 'cl' + str(cluster_num))]
                if ci.size == 0:
                    continue
                if ci['quality'] == 'NOISE':
                    continue
                row = pd.DataFrame(ci)
                df = pd.concat([df, row], axis=0, ignore_index=True)

    return df
