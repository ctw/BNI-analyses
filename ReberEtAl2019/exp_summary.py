"""
Generate summary of clusters analyzed in Reber et al. (2019).
"""

from glob import glob
import logging
import os

import h5py
import numpy as np
import pandas as pd

# map from Reber qualities to BML
qualityMap = {'SU': 'SPIKE', 'MU': 'POTENTIAL'}

# map from Reber regions to BML
areaMap = {'AM':'A', 'HC':'H', 'EC':'EC', 'PHC':'PHC', 'other':'NA'}


def expSummary(h5FileName):
    data = h5py.File(h5FileName, 'r')

    res = pd.DataFrame()
    i = 0
    for subj in np.unique(data['cluster_lookup'][:]['subjid']):
        for sess in np.unique(data['cluster_lookup'][:][data['cluster_lookup'][:][
                                                                'subjid'] == subj]['sessid']):
            expstr = 'reber_' + str(subj) + '_' + str(sess)
            logging.info(f"{expstr}")

            filt = ((data['cluster_lookup'][:]['subjid']==subj) &
                    (data['cluster_lookup'][:]['sessid']==sess))

            for cluster in data['cluster_lookup'][filt]:
                channum = cluster['channo']
                cluster_num = cluster['clusid']
                clusterId = f"ch{channum}cl{cluster_num}"
                region = cluster['regionname'].decode('utf-8')
                area = areaMap[cluster['regionname'].decode('utf-8')]
                hemisphere = cluster['hemisphere'].decode('utf-8')
                quality = qualityMap[cluster['clustype'].decode('utf-8')]

                row = pd.DataFrame(data={'expstr':expstr, 'cluster_id':clusterId,
                                         'hemisphere':hemisphere,  'area':area,
                                         'quality':quality},
                                   index=[i])
                res = pd.concat([res, row], axis=0, ignore_index=True)
                i += 1

    return res