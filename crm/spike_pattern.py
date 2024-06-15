
import numpy as np
from glob import glob
#import os
#from scipy import signal
#import bni_helpers
#from pathlib import Path

import h5py
#from numpy.lib import recfunctions
from ptsa.data.timeseries import TimeSeries

from sklearn.cluster import KMeans
# from sklearn.cluster import MeanShift
import itertools


evs = np.load('/home/ctw/data/BNI/data/bni_evs_crm.npz')['evs'].view(
    np.recarray)
#evs = recfunctions.append_fields(evs, 'missing_neuro', [False]*len(evs),
#                                 usemask=False, asrecarray=True)

rhino_root = '/home/ctw/fusemounts/rhino'
nse_dtype = [('timestamp', 'uint64'), ('channel', 'uint32'),
             ('unit_id', 'uint32'), ('params', 'uint32', (8,)),
             ('samples', 'uint16', (32,))]
header_size=16*1024
                


# targvecs_all = []
all_vecs = []
for i, expstr in enumerate(np.unique(evs['expstr'])):
    print(i+1, '/', len(np.unique(evs['expstr'])), expstr)
    ev_sess = evs[evs['expstr']==expstr]
    assert len(np.unique(ev_sess['subject']))==1, 'subj len'
    subject = ev_sess['subject'][0]
    nsefiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                           subject + '/analysis/' + expstr +
                           '/KK/CSC?.Nse'))
    nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                                subject + '/analysis/' + expstr +
                                '/KK/CSC??.Nse')))
    # not really needed, but just in case:
    nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                                subject + '/analysis/' + expstr +
                                '/KK/CSC???.Nse')))
    
    clusterfiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                               subject + '/analysis/' + expstr +
                               '/KK/CSC?.clu*'))
    clusterfiles.extend(sorted(glob(
        rhino_root+'/scratch/josh/BniData/Subjects/' +
        subject + '/analysis/' + expstr + '/KK/CSC??.clu*')))
    # not really needed, but just in case:
    clusterfiles.extend(sorted(glob(
        rhino_root+'/scratch/josh/BniData/Subjects/' +
        subject + '/analysis/' + expstr + '/KK/CSC???.clu*')))
    if expstr == 's23e15cr':
        # CSC5.clu.1 and CSC5.clu.1 are missing, so we're inserting
        # Nones to make sure clusterfiles and nsefiles align (will
        # skip these later):
        clusterfiles.insert(4, None)
        clusterfiles.insert(4, None)
        # there are no NCS files for channels 41-48 (corresponding to
        # indices 31-39) so we're removing them (starting with the
        # last one to be able to remove a range of indices rather than
        # removing the same index multiple times):
        for i in range(39, 31, -1):
            clusterfiles.pop(i)
            nsefiles.pop(i)
    assert len(nsefiles) == len(clusterfiles),\
        'len(nsefiles) == len(clusterfiles)'
    assert np.alltrue([int(clusterfile.split('CSC')[-1].split('.')[0]) ==
                       int(nsefile.split('CSC')[-1].split('.')[0])
                       for clusterfile, nsefile in
                       zip(clusterfiles, nsefiles)
                       if clusterfile is not None]), 'compare nses & clusters'
    clusterinfo_file = (rhino_root+'/scratch/josh/BniData/Subjects/' +
                        subject + '/analysis/' + expstr +
                        '/clusterInfo.txt')
    clusterinfo = None
    if os.path.exists(clusterinfo_file):
        clusterinfo = np.loadtxt(clusterinfo_file, delimiter='\t', skiprows=1,
                                 dtype={'names': (
                                     'expstr', 'cluster_id',
                                     'hemisphere', 'area', 'quality'),
                                        'formats': ('U16', 'U16', 'U8',
                                                    'U8', 'U16')})
        assert len(np.unique(clusterinfo['expstr'])) == 1,\
            'clusterinfo expstr 1'
        assert clusterinfo['expstr'][0] == expstr, 'clusterinfo expstr 2'
        # if expstr != 's21e12cr':s23e11cr
        #     # not sure what's wrong with clusterinfo for 's21e12cr'
        #     assert np.alltrue(np.unique([
        #         int(clid.split('cl')[0].split('ch')[1])
        #         for clid in clusterinfo['cluster_id']]) == np.array([
        #                 int(nsefile.split(
        #                     'CSC')[-1].split('.')[
        #                         0]) for nsefile in nsefiles])), 'chan chk'
    else:
        raise ValueError('No clusterinfo file!')
    # print('******************************************')
    # print('******************************************')
    # print('******************************************')
    # print('********** NO CLUSTERINFO FILE! **********')
    # print('******************************************')
    # print('******************************************')
    # print('******************************************')
    # os.remove('/home/ctw/data/BNI/data/time_series/'+expstr+'_ts.hdf')
    # continue

    sme_filt = ev_sess['first'] &  (ev_sess['paired'] != 'unpaired')
    sme_resp = np.empty(np.sum(sme_filt), dtype = 'U3')
    sme_bins = {'old': [], 'new': []}
    sme_timewindowlength = 1500000  # 1.5 s
    for sme_indx, sme_trial in enumerate(ev_sess[sme_filt]):
        # np.where(ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'],
        #          ev_sess[~sme_filt], [None]*np.sum(~sme_filt))
        #indx = np.flatnonzero(
        #    ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'])
        indices = np.flatnonzero(
            ev_sess['stim_word'] == sme_trial['stim_word'])
        assert len(indices) == 2, 'len(indices) == 2'
        old_resp = np.sum(
            ev_sess[indices[1]]['resp_count_old_1']+
            ev_sess[indices[1]]['resp_count_old_2']+
            ev_sess[indices[1]]['resp_count_old_3'])
        new_resp = np.sum(
            ev_sess[indices[1]]['resp_count_new_1']+
            ev_sess[indices[1]]['resp_count_new_2']+
            ev_sess[indices[1]]['resp_count_new_3'])
        if (old_resp > 0) & (new_resp == 0):
            sme_resp[sme_indx] = 'old'
            sme_bins['old'].extend([ev_sess[indices[0]]['phase1_time'],
                                    ev_sess[indices[0]]['phase1_time'] +
                                    sme_timewindowlength])
        elif (old_resp == 0) & (new_resp > 0):
            sme_resp[sme_indx] = 'new'
            sme_bins['new'].extend([ev_sess[indices[0]]['phase1_time'],
                                    ev_sess[indices[0]]['phase1_time'] +
                                    sme_timewindowlength])
        else:
            sme_resp[sme_indx] = '---'
    if (np.sum(sme_resp=='old') < 1) or (np.sum(sme_resp=='new') < 1):
        continue
    # spikesums = {'old': np.empty(np.sum(sme_resp=='old'), int),
    #             'new': np.empty(np.sum(sme_resp=='new'), int)}
    spikesums = {'old': [], 'new': []}
    spikesums_cont = []
    #binwidth = 100000  # 100 ms
    # binwidth = 200000  # 200 ms
    binwidth = 100000  # 100 ms
    binwidth = 500000  # 500 ms
    time_offset = 5000000  # 5s
    bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
                     ev_sess['phase3_time'][-1]+time_offset, binwidth)
    channels = []
    clusters = []
    clusterIDs = []
    for clusterfile, nsefile in zip(
            clusterfiles, nsefiles):
        if clusterfile is None:
            continue
        chan = nsefile.split('.Nse')[0]
        channum = int(chan.split('CSC')[-1])
        # # clusters[channum] = {}
        # cluster[chan] = {}
        nsetimes = np.memmap(nsefile, dtype=nse_dtype, mode='r',
                             offset=header_size)['timestamp']
        numclusts = np.fromfile(clusterfile, dtype=int, sep='\n')[0]
        sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]
        for cluster_num in np.unique(sess_clusts):
            ci = clusterinfo[clusterinfo[
                'cluster_id'] == ('ch'+str(channum)+'cl'+str(cluster_num))]
            if ci['quality'] == 'NOISE':
                #print(channum, cluster_num)
                continue
            # cluster_num = int(cluster_num)
            clustimes = nsetimes[sess_clusts == cluster_num]
            for oldnew in sme_bins:
                hist, _ = np.histogram(clustimes, sme_bins[oldnew])
                # every other bin is the time between item presentations,
                # so we're skipping these:
                spikesums[oldnew].append(hist[::2])
            channels.append(channum)
            clusters.append(cluster_num)
            clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
            hist, _ = np.histogram(clustimes, bins)
            spikesums_cont.append(hist)
    times_cont = bins[:-1] + binwidth/2
    spikesums_cont = np.array(spikesums_cont)
    spikesums_cont_norm = spikesums_cont/np.atleast_2d(spikesums_cont.max(1)).T
    firing_rate_cont = spikesums_cont/(binwidth/1000000)
    times = {}
    spikesums_norm = {}
    firing_rate = {}
    for oldnew in spikesums:
        spikesums[oldnew] = np.array(spikesums[oldnew])
        spikesums_norm[oldnew] = spikesums[oldnew]/np.atleast_2d(
            spikesums[oldnew].max(1)).T
        spikesums_norm[oldnew][spikesums[oldnew].max(1)==0] = 0
        firing_rate[oldnew] = spikesums[oldnew]/(sme_timewindowlength/1000000)
        times[oldnew] = np.array(sme_bins[oldnew][::2]) + sme_timewindowlength/2
    channels = np.array(channels)
    cluster = np.array(clusters)
    clusterIDs = np.array(clusterIDs)
    # target_vectors = {'expstr': expstr}
    target_vectors = {}
    for oldnew in firing_rate:
        target_vectors[oldnew] = firing_rate[oldnew].mean(1)
        target_vectors[oldnew+'_std'] = firing_rate[oldnew].std(1)
        target_vectors[oldnew+'_n'] = firing_rate[oldnew].shape[1]
        #target_vectors[oldnew+'_norm'] = target_vectors[oldnew]/np.sqrt(
        #    (target_vectors[oldnew]**2).sum())
    #target_vectors['oldnew'] = target_vectors['old']-target_vectors['new']
    #target_vectors['oldnew_norm'] = target_vectors['oldnew']/np.sqrt(
    #    (target_vectors['oldnew']**2).sum())
    # targvecs_all.append(target_vectors)
    all_vecs.append({'expstr': expstr, 'times_cont': times_cont,
                     'spikesums_cont': spikesums_cont,
                     'spikesums_cont_norm': spikesums_cont_norm, 'times': times,
                     'firing_rate_cont': firing_rate_cont, 'clusters': clusters,
                     'clusterIDs': clusterIDs, 'channels': channels,
                     'spikesums_norm': spikesums_norm, 'firing_rate': firing_rate,
                     'target_vectors': target_vectors})



import matplotlib.pyplot as plt
plt.ion()

all_figsaxs = []
for i in range(len(all_vecs)):
    norm_vecs = {}
    for oldnew in ['old', 'new']:
        norm_vecs[oldnew] = all_vecs[i]['target_vectors'][oldnew]/np.sqrt((
            all_vecs[i]['target_vectors'][oldnew]**2).sum())
    norm_vecs['diff'] = (all_vecs[i]['target_vectors']['old']-
                         all_vecs[i]['target_vectors']['new'])/np.sqrt(
                             ((all_vecs[i]['target_vectors']['old']-
                               all_vecs[i]['target_vectors']['new'])**2).sum())
    av_norm = all_vecs[i]['firing_rate_cont']/np.sqrt(
        (all_vecs[i]['firing_rate_cont']**2).sum(axis=0))
    av_norm[~np.isfinite(av_norm)]=0
    # match = np.dot(norm_vecs['diff'], av_norm)
    # match = np.dot(norm_vecs['old'], av_norm)
    match = np.dot(norm_vecs['new'], av_norm)
    # minval = match.min() - np.abs(match.min())*0.05
    # maxval = match.min() - np.abs(match.min())*0.05
    if i%15==0:
        all_figsaxs.append(plt.subplots(15, 1))
    all_figsaxs[-1][1][i%15].plot(match, alpha=0.5)
    ev_sess = evs[evs['expstr']==all_vecs[i]['expstr']]
    tsttimes = all_vecs[i]['times']['new']-ev_sess['phase1_time'][0]+time_offset
    # xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    xs = np.ravel([(x, x+3, x+3) for x in np.round(tsttimes/1000/500)])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='yellow', alpha=0.4)
    tsttimes = all_vecs[i]['times']['old']-ev_sess['phase1_time'][0]+time_offset
    #xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    xs = np.ravel([(x, x+3, x+3) for x in np.round(tsttimes/1000/500)])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='red', alpha=0.4)













import scipy as sp

all_figsaxs = []
for i in range(len(all_vecs)):
    mahalanobis = {}
    mahalanobis['cov'] = np.cov(all_vecs[i]['firing_rate']['old'])
    mahalanobis['cov_inv'] = sp.linalg.inv(mahalanobis['cov'])
    mahalanobis['mean'] = all_vecs[i]['firing_rate']['old'].mean(1)
    mahalanobis['delta'] = (all_vecs[i]['firing_rate_cont']-
                            np.atleast_2d(mahalanobis['mean']).T)
    mahalanobis['distance'] = np.array(np.diag(np.dot(
        np.dot(mahalanobis['delta'].T, mahalanobis['cov_inv']),
        mahalanobis['delta'])))
    mahalanobis['distance'][mahalanobis['distance']<0]=0
    mahalanobis['distance'] = np.sqrt(mahalanobis['distance'])
    if i%15==0:
        all_figsaxs.append(plt.subplots(15, 1))
    all_figsaxs[-1][1][i%15].plot(mahalanobis['distance'], alpha=0.5)
    ev_sess = evs[evs['expstr']==all_vecs[i]['expstr']]
    tsttimes = all_vecs[i]['times']['new']-ev_sess['phase1_time'][0]+time_offset
    # xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    xs = np.ravel([(x, x+3, x+3) for x in np.round(tsttimes/1000/500)])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='yellow', alpha=0.4)
    tsttimes = all_vecs[i]['times']['old']-ev_sess['phase1_time'][0]+time_offset
    #xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    xs = np.ravel([(x, x+3, x+3) for x in np.round(tsttimes/1000/500)])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='red', alpha=0.4)









# # targvecs_all = []
# all_vecs = []
# for i, expstr in enumerate(np.unique(evs['expstr'])):
#     print(i+1, '/', len(np.unique(evs['expstr'])), expstr)
#     ev_sess = evs[evs['expstr']==expstr]
#     assert len(np.unique(ev_sess['subject']))==1, 'subj len'
#     subject = ev_sess['subject'][0]
#     nsefiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
#                            subject + '/analysis/' + expstr +
#                            '/KK/CSC?.Nse'))
#     nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
#                                 subject + '/analysis/' + expstr +
#                                 '/KK/CSC??.Nse')))
#     # not really needed, but just in case:
#     nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
#                                 subject + '/analysis/' + expstr +
#                                 '/KK/CSC???.Nse')))
    
#     clusterfiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
#                                subject + '/analysis/' + expstr +
#                                '/KK/CSC?.clu*'))
#     clusterfiles.extend(sorted(glob(
#         rhino_root+'/scratch/josh/BniData/Subjects/' +
#         subject + '/analysis/' + expstr + '/KK/CSC??.clu*')))
#     # not really needed, but just in case:
#     clusterfiles.extend(sorted(glob(
#         rhino_root+'/scratch/josh/BniData/Subjects/' +
#         subject + '/analysis/' + expstr + '/KK/CSC???.clu*')))
#     if expstr == 's23e15cr':
#         # CSC5.clu.1 and CSC5.clu.1 are missing, so we're inserting
#         # Nones to make sure clusterfiles and nsefiles align (will
#         # skip these later):
#         clusterfiles.insert(4, None)
#         clusterfiles.insert(4, None)
#         # there are no NCS files for channels 41-48 (corresponding to
#         # indices 31-39) so we're removing them (starting with the
#         # last one to be able to remove a range of indices rather than
#         # removing the same index multiple times):
#         for i in range(39, 31, -1):
#             clusterfiles.pop(i)
#             nsefiles.pop(i)
#     assert len(nsefiles) == len(clusterfiles),\
#         'len(nsefiles) == len(clusterfiles)'
#     assert np.alltrue([int(clusterfile.split('CSC')[-1].split('.')[0]) ==
#                        int(nsefile.split('CSC')[-1].split('.')[0])
#                        for clusterfile, nsefile in
#                        zip(clusterfiles, nsefiles)
#                        if clusterfile is not None]), 'compare nses & clusters'
#     clusterinfo_file = (rhino_root+'/scratch/josh/BniData/Subjects/' +
#                         subject + '/analysis/' + expstr +
#                         '/clusterInfo.txt')
#     clusterinfo = None
#     if os.path.exists(clusterinfo_file):
#         clusterinfo = np.loadtxt(clusterinfo_file, delimiter='\t', skiprows=1,
#                                  dtype={'names': (
#                                      'expstr', 'cluster_id',
#                                      'hemisphere', 'area', 'quality'),
#                                         'formats': ('U16', 'U16', 'U8',
#                                                     'U8', 'U16')})
#         assert len(np.unique(clusterinfo['expstr'])) == 1,\
#             'clusterinfo expstr 1'
#         assert clusterinfo['expstr'][0] == expstr, 'clusterinfo expstr 2'
#         # if expstr != 's21e12cr':s23e11cr
#         #     # not sure what's wrong with clusterinfo for 's21e12cr'
#         #     assert np.alltrue(np.unique([
#         #         int(clid.split('cl')[0].split('ch')[1])
#         #         for clid in clusterinfo['cluster_id']]) == np.array([
#         #                 int(nsefile.split(
#         #                     'CSC')[-1].split('.')[
#         #                         0]) for nsefile in nsefiles])), 'chan chk'
#     else:
#         raise ValueError('No clusterinfo file!')
#     # print('******************************************')
#     # print('******************************************')
#     # print('******************************************')
#     # print('********** NO CLUSTERINFO FILE! **********')
#     # print('******************************************')
#     # print('******************************************')
#     # print('******************************************')
#     # os.remove('/home/ctw/data/BNI/data/time_series/'+expstr+'_ts.hdf')
#     # continue

#     sme_filt = ev_sess['first'] &  (ev_sess['paired'] != 'unpaired')
#     sme_resp = np.empty(np.sum(sme_filt), dtype = 'U3')
#     sme_bins = {'old': [], 'new': []}
#     sme_timewindowlength = 1500000  # 1.5 s
#     for sme_indx, sme_trial in enumerate(ev_sess[sme_filt]):
#         # np.where(ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'],
#         #          ev_sess[~sme_filt], [None]*np.sum(~sme_filt))
#         #indx = np.flatnonzero(
#         #    ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'])
#         indices = np.flatnonzero(
#             ev_sess['stim_word'] == sme_trial['stim_word'])
#         assert len(indices) == 2, 'len(indices) == 2'
#         old_resp = np.sum(
#             ev_sess[indices[1]]['resp_count_old_1']+
#             ev_sess[indices[1]]['resp_count_old_2']+
#             ev_sess[indices[1]]['resp_count_old_3'])
#         new_resp = np.sum(
#             ev_sess[indices[1]]['resp_count_new_1']+
#             ev_sess[indices[1]]['resp_count_new_2']+
#             ev_sess[indices[1]]['resp_count_new_3'])
#         if (old_resp > 0) & (new_resp == 0):
#             sme_resp[sme_indx] = 'old'
#             sme_bins['old'].extend([ev_sess[indices[0]]['phase1_time'],
#                                     ev_sess[indices[0]]['phase1_time'] +
#                                     sme_timewindowlength])
#         elif (old_resp == 0) & (new_resp > 0):
#             sme_resp[sme_indx] = 'new'
#             sme_bins['new'].extend([ev_sess[indices[0]]['phase1_time'],
#                                     ev_sess[indices[0]]['phase1_time'] +
#                                     sme_timewindowlength])
#         else:
#             sme_resp[sme_indx] = '---'
#     if (np.sum(sme_resp=='old') < 1) or (np.sum(sme_resp=='new') < 1):
#         continue
#     # spikesums = {'old': np.empty(np.sum(sme_resp=='old'), int),
#     #             'new': np.empty(np.sum(sme_resp=='new'), int)}
#     spikesums = {'old': [], 'new': []}
#     spikesums_cont = []
#     binwidth = 100000  # 100 ms
#     # binwidth = 200000  # 200 ms
#     time_offset = 5000000  # 5s
#     bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
#                      ev_sess['phase3_time'][-1]+time_offset, binwidth)
#     channels = []
#     clusters = []
#     clusterIDs = []
#     for clusterfile, nsefile in zip(
#             clusterfiles, nsefiles):
#         if clusterfile is None:
#             continue
#         chan = nsefile.split('.Nse')[0]
#         channum = int(chan.split('CSC')[-1])
#         # # clusters[channum] = {}
#         # cluster[chan] = {}
#         nsetimes = np.memmap(nsefile, dtype=nse_dtype, mode='r',
#                              offset=header_size)['timestamp']
#         numclusts = np.fromfile(clusterfile, dtype=int, sep='\n')[0]
#         sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]
#         for cluster_num in np.unique(sess_clusts):
#             ci = clusterinfo[clusterinfo[
#                 'cluster_id'] == ('ch'+str(channum)+'cl'+str(cluster_num))]
#             if ci['quality'] == 'NOISE':
#                 #print(channum, cluster_num)
#                 continue
#             # cluster_num = int(cluster_num)
#             clustimes = nsetimes[sess_clusts == cluster_num]
#             for oldnew in sme_bins:
#                 hist, _ = np.histogram(clustimes, sme_bins[oldnew])
#                 # every other bin is the time between item presentations,
#                 # so we're skipping these:
#                 spikesums[oldnew].append(hist[::2])
#             channels.append(channum)
#             clusters.append(cluster_num)
#             clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
#             hist, _ = np.histogram(clustimes, bins)
#             spikesums_cont.append(hist)
#     times_cont = bins[:-1] + binwidth/2
#     spikesums_cont = np.array(spikesums_cont)
#     spikesums_cont_norm = spikesums_cont/np.atleast_2d(spikesums_cont.max(1)).T
#     firing_rate_cont = spikesums_cont/(binwidth/1000000)
#     times = {}
#     spikesums_norm = {}
#     firing_rate = {}
#     for oldnew in spikesums:
#         spikesums[oldnew] = np.array(spikesums[oldnew])
#         spikesums_norm[oldnew] = spikesums[oldnew]/np.atleast_2d(
#             spikesums[oldnew].max(1)).T
#         spikesums_norm[oldnew][spikesums[oldnew].max(1)==0] = 0
#         firing_rate[oldnew] = spikesums[oldnew]/(sme_timewindowlength/1000000)
#         times[oldnew] = np.array(sme_bins[oldnew][::2]) + sme_timewindowlength/2
#     channels = np.array(channels)
#     cluster = np.array(clusters)
#     clusterIDs = np.array(clusterIDs)
#     # target_vectors = {'expstr': expstr}
#     target_vectors = {}
#     for oldnew in firing_rate:
#         target_vectors[oldnew] = firing_rate[oldnew].mean(1)
#         target_vectors[oldnew+'_std'] = firing_rate[oldnew].std(1)
#         target_vectors[oldnew+'_n'] = firing_rate[oldnew].shape[1]
#         #target_vectors[oldnew+'_norm'] = target_vectors[oldnew]/np.sqrt(
#         #    (target_vectors[oldnew]**2).sum())
#     #target_vectors['oldnew'] = target_vectors['old']-target_vectors['new']
#     #target_vectors['oldnew_norm'] = target_vectors['oldnew']/np.sqrt(
#     #    (target_vectors['oldnew']**2).sum())
#     # targvecs_all.append(target_vectors)
#     all_vecs.append({'expstr': expstr, 'times_cont': times_cont,
#                      'spikesums_cont': spikesums_cont,
#                      'spikesums_cont_norm': spikesums_cont_norm, 'times': times,
#                      'firing_rate_cont': firing_rate_cont, 'clusters': clusters,
#                      'clusterIDs': clusterIDs, 'channels': channels,
#                      'spikesums_norm': spikesums_norm, 'firing_rate': firing_rate,
#                      'target_vectors': target_vectors})




import matplotlib.pyplot as plt
plt.ion()

all_figsaxs_len = []
for i in range(len(all_vecs)):
    norm_vecs = {}
    for oldnew in ['old', 'new']:
        norm_vecs[oldnew] = all_vecs[i]['target_vectors'][oldnew]/np.sqrt((
            all_vecs[i]['target_vectors'][oldnew]**2).sum())
    norm_vecs['diff'] = (all_vecs[i]['target_vectors']['old']-
                         all_vecs[i]['target_vectors']['new'])/np.sqrt(
                             ((all_vecs[i]['target_vectors']['old']-
                               all_vecs[i]['target_vectors']['new'])**2).sum())
    av_len = np.sqrt(
        (all_vecs[i]['firing_rate_cont']**2).sum(axis=0))
    if i%15==0:
        all_figsaxs_len.append(plt.subplots(15, 1))
    all_figsaxs_len[-1][1][i%15].plot(av_len, alpha=0.5)
    ev_sess = evs[evs['expstr']==all_vecs[i]['expstr']]
    tsttimes = all_vecs[i]['times']['new']-ev_sess['phase1_time'][0]+time_offset
    xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    y1s = [av_len.max()]*len(xs)
    y2s = [av_len.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs_len[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='yellow', alpha=0.4)
    tsttimes = all_vecs[i]['times']['old']-ev_sess['phase1_time'][0]+time_offset
    xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    y1s = [av_len.max()]*len(xs)
    y2s = [av_len.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs_len[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='red', alpha=0.4)











tst = np.c_[spikesums['old'], spikesums['new']]
tstz = (tst - np.atleast_2d(tst.mean(1)).T) / np.atleast_2d(tst.std(1)).T
tstz -= tstz.mean(0)
tstz /= tstz.std(0)
spikesums_z = {'old': tstz[:, :spikesums['old'].shape[1]],
               'new': tstz[:, spikesums['old'].shape[1]:]}
spikesums_cont_z = ((spikesums_cont - np.atleast_2d(spikesums_cont.mean(1)).T) /
                    np.atleast_2d(spikesums_cont.std(1)).T)
spikesums_cont_z -= spikesums_cont_z.mean(0)
spikesums_cont_z /= spikesums_cont_z.std(0)






# targvecs_all = []
all_vecs_z = []
for i, expstr in enumerate(np.unique(evs['expstr'])):
    print(i+1, '/', len(np.unique(evs['expstr'])), expstr)
    ev_sess = evs[evs['expstr']==expstr]
    assert len(np.unique(ev_sess['subject']))==1, 'subj len'
    subject = ev_sess['subject'][0]
    nsefiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                           subject + '/analysis/' + expstr +
                           '/KK/CSC?.Nse'))
    nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                                subject + '/analysis/' + expstr +
                                '/KK/CSC??.Nse')))
    # not really needed, but just in case:
    nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                                subject + '/analysis/' + expstr +
                                '/KK/CSC???.Nse')))
    
    clusterfiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                               subject + '/analysis/' + expstr +
                               '/KK/CSC?.clu*'))
    clusterfiles.extend(sorted(glob(
        rhino_root+'/scratch/josh/BniData/Subjects/' +
        subject + '/analysis/' + expstr + '/KK/CSC??.clu*')))
    # not really needed, but just in case:
    clusterfiles.extend(sorted(glob(
        rhino_root+'/scratch/josh/BniData/Subjects/' +
        subject + '/analysis/' + expstr + '/KK/CSC???.clu*')))
    if expstr == 's23e15cr':
        # CSC5.clu.1 and CSC5.clu.1 are missing, so we're inserting
        # Nones to make sure clusterfiles and nsefiles align (will
        # skip these later):
        clusterfiles.insert(4, None)
        clusterfiles.insert(4, None)
        # there are no NCS files for channels 41-48 (corresponding to
        # indices 31-39) so we're removing them (starting with the
        # last one to be able to remove a range of indices rather than
        # removing the same index multiple times):
        for i in range(39, 31, -1):
            clusterfiles.pop(i)
            nsefiles.pop(i)
    assert len(nsefiles) == len(clusterfiles),\
        'len(nsefiles) == len(clusterfiles)'
    assert np.alltrue([int(clusterfile.split('CSC')[-1].split('.')[0]) ==
                       int(nsefile.split('CSC')[-1].split('.')[0])
                       for clusterfile, nsefile in
                       zip(clusterfiles, nsefiles)
                       if clusterfile is not None]), 'compare nses & clusters'
    clusterinfo_file = (rhino_root+'/scratch/josh/BniData/Subjects/' +
                        subject + '/analysis/' + expstr +
                        '/clusterInfo.txt')
    clusterinfo = None
    if os.path.exists(clusterinfo_file):
        clusterinfo = np.loadtxt(clusterinfo_file, delimiter='\t', skiprows=1,
                                 dtype={'names': (
                                     'expstr', 'cluster_id',
                                     'hemisphere', 'area', 'quality'),
                                        'formats': ('U16', 'U16', 'U8',
                                                    'U8', 'U16')})
        assert len(np.unique(clusterinfo['expstr'])) == 1,\
            'clusterinfo expstr 1'
        assert clusterinfo['expstr'][0] == expstr, 'clusterinfo expstr 2'
        # if expstr != 's21e12cr':s23e11cr
        #     # not sure what's wrong with clusterinfo for 's21e12cr'
        #     assert np.alltrue(np.unique([
        #         int(clid.split('cl')[0].split('ch')[1])
        #         for clid in clusterinfo['cluster_id']]) == np.array([
        #                 int(nsefile.split(
        #                     'CSC')[-1].split('.')[
        #                         0]) for nsefile in nsefiles])), 'chan chk'
    else:
        raise ValueError('No clusterinfo file!')
    # print('******************************************')
    # print('******************************************')
    # print('******************************************')
    # print('********** NO CLUSTERINFO FILE! **********')
    # print('******************************************')
    # print('******************************************')
    # print('******************************************')
    # os.remove('/home/ctw/data/BNI/data/time_series/'+expstr+'_ts.hdf')
    # continue

    sme_filt = ev_sess['first'] &  (ev_sess['paired'] != 'unpaired')
    sme_resp = np.empty(np.sum(sme_filt), dtype = 'U3')
    sme_bins = {'old': [], 'new': []}
    sme_timewindowlength = 1500000  # 1.5 s
    for sme_indx, sme_trial in enumerate(ev_sess[sme_filt]):
        # np.where(ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'],
        #          ev_sess[~sme_filt], [None]*np.sum(~sme_filt))
        #indx = np.flatnonzero(
        #    ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'])
        indices = np.flatnonzero(
            ev_sess['stim_word'] == sme_trial['stim_word'])
        assert len(indices) == 2, 'len(indices) == 2'
        old_resp = np.sum(
            ev_sess[indices[1]]['resp_count_old_1']+
            ev_sess[indices[1]]['resp_count_old_2']+
            ev_sess[indices[1]]['resp_count_old_3'])
        new_resp = np.sum(
            ev_sess[indices[1]]['resp_count_new_1']+
            ev_sess[indices[1]]['resp_count_new_2']+
            ev_sess[indices[1]]['resp_count_new_3'])
        if (old_resp > 0) & (new_resp == 0):
            sme_resp[sme_indx] = 'old'
            sme_bins['old'].extend([ev_sess[indices[0]]['phase1_time'],
                                    ev_sess[indices[0]]['phase1_time'] +
                                    sme_timewindowlength])
        elif (old_resp == 0) & (new_resp > 0):
            sme_resp[sme_indx] = 'new'
            sme_bins['new'].extend([ev_sess[indices[0]]['phase1_time'],
                                    ev_sess[indices[0]]['phase1_time'] +
                                    sme_timewindowlength])
        else:
            sme_resp[sme_indx] = '---'
    if (np.sum(sme_resp=='old') < 1) or (np.sum(sme_resp=='new') < 1):
        continue
    # spikesums = {'old': np.empty(np.sum(sme_resp=='old'), int),
    #             'new': np.empty(np.sum(sme_resp=='new'), int)}
    spikesums = {'old': [], 'new': []}
    spikesums_cont = []
    binwidth = 100000  # 100 ms
    # binwidth = 200000  # 200 ms
    time_offset = 5000000  # 5s
    bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
                     ev_sess['phase3_time'][-1]+time_offset, binwidth)
    channels = []
    clusters = []
    clusterIDs = []
    for clusterfile, nsefile in zip(
            clusterfiles, nsefiles):
        if clusterfile is None:
            continue
        chan = nsefile.split('.Nse')[0]
        channum = int(chan.split('CSC')[-1])
        # # clusters[channum] = {}
        # cluster[chan] = {}
        nsetimes = np.memmap(nsefile, dtype=nse_dtype, mode='r',
                             offset=header_size)['timestamp']
        numclusts = np.fromfile(clusterfile, dtype=int, sep='\n')[0]
        sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]
        for cluster_num in np.unique(sess_clusts):
            ci = clusterinfo[clusterinfo[
                'cluster_id'] == ('ch'+str(channum)+'cl'+str(cluster_num))]
            if ci['quality'] == 'NOISE':
                #print(channum, cluster_num)
                continue
            # cluster_num = int(cluster_num)
            clustimes = nsetimes[sess_clusts == cluster_num]
            for oldnew in sme_bins:
                hist, _ = np.histogram(clustimes, sme_bins[oldnew])
                # every other bin is the time between item presentations,
                # so we're skipping these:
                spikesums[oldnew].append(hist[::2])
            channels.append(channum)
            clusters.append(cluster_num)
            clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
            hist, _ = np.histogram(clustimes, bins)
            spikesums_cont.append(hist)
    times_cont = bins[:-1] + binwidth/2
    spikesums_cont = np.array(spikesums_cont)
    tst = np.c_[spikesums['old'], spikesums['new']]
    tstz = (tst - np.atleast_2d(tst.mean(1)).T) / np.atleast_2d(tst.std(1)).T
    tstz -= tstz.mean(0)
    tstz /= tstz.std(0)
    spikesums_z = {'old': tstz[:, :np.shape(spikesums['old'])[1]],
                   'new': tstz[:, np.shape(spikesums['old'])[1]:]}
    spikesums_cont_z = ((spikesums_cont - np.atleast_2d(spikesums_cont.mean(1)).T) /
                        np.atleast_2d(spikesums_cont.std(1)).T)
    spikesums_cont_z -= spikesums_cont_z.mean(0)
    spikesums_cont_z /= spikesums_cont_z.std(0)
    spikesums_cont_norm_z = spikesums_cont_z/np.atleast_2d(spikesums_cont_z.max(1)).T
    firing_rate_cont_z = spikesums_cont_z/(binwidth/1000000)
    times = {}
    spikesums_norm_z = {}
    firing_rate_z = {}
    for oldnew in spikesums_z:
        # spikesums_z[oldnew] = np.array(spikesums[oldnew])
        spikesums_norm_z[oldnew] = spikesums_z[oldnew]/np.atleast_2d(
            spikesums_z[oldnew].max(1)).T
        spikesums_norm_z[oldnew][spikesums_z[oldnew].max(1)==0] = 0
        firing_rate_z[oldnew] = spikesums_z[oldnew]/(sme_timewindowlength/1000000)
        times[oldnew] = np.array(sme_bins[oldnew][::2]) + sme_timewindowlength/2
    channels = np.array(channels)
    cluster = np.array(clusters)
    clusterIDs = np.array(clusterIDs)
    # target_vectors = {'expstr': expstr}
    target_vectors_z = {}
    for oldnew in firing_rate_z:
        target_vectors_z[oldnew] = firing_rate_z[oldnew].mean(1)
        target_vectors_z[oldnew+'_std'] = firing_rate_z[oldnew].std(1)
        target_vectors_z[oldnew+'_n'] = firing_rate_z[oldnew].shape[1]
        #target_vectors[oldnew+'_norm'] = target_vectors[oldnew]/np.sqrt(
        #    (target_vectors[oldnew]**2).sum())
    #target_vectors['oldnew'] = target_vectors['old']-target_vectors['new']
    #target_vectors['oldnew_norm'] = target_vectors['oldnew']/np.sqrt(
    #    (target_vectors['oldnew']**2).sum())
    # targvecs_all.append(target_vectors)
    all_vecs_z.append({'expstr': expstr, 'times_cont': times_cont,
                     'spikesums_cont_z': spikesums_cont_z,
                     'spikesums_cont_norm_z': spikesums_cont_norm_z, 'times': times,
                     'firing_rate_cont_z': firing_rate_cont_z, 'clusters': clusters,
                     'clusterIDs': clusterIDs, 'channels': channels,
                     'spikesums_norm_z': spikesums_norm_z,
                     'firing_rate_z': firing_rate_z,
                     'target_vectors_z': target_vectors_z})





all_figsaxs_z = []
for i in range(len(all_vecs)):
    norm_vecs = {}
    for oldnew in ['old', 'new']:
        norm_vecs[oldnew] = all_vecs_z[i]['target_vectors_z'][oldnew]/np.sqrt((
            all_vecs_z[i]['target_vectors_z'][oldnew]**2).sum())
    norm_vecs['diff'] = (all_vecs_z[i]['target_vectors_z']['old']-
                         all_vecs_z[i]['target_vectors_z']['new'])/np.sqrt(
                             ((all_vecs_z[i]['target_vectors_z']['old']-
                               all_vecs_z[i]['target_vectors_z']['new'])**2).sum())
    av_norm = all_vecs_z[i]['firing_rate_cont_z']/np.sqrt(
        (all_vecs_z[i]['firing_rate_cont_z']**2).sum(axis=0))
    av_norm[~np.isfinite(av_norm)]=0
    match = np.dot(norm_vecs['diff'], av_norm)
    # minval = match.min() - np.abs(match.min())*0.05
    # maxval = match.min() - np.abs(match.min())*0.05
    if i%15==0:
        all_figsaxs_z.append(plt.subplots(15, 1))
    all_figsaxs_z[-1][1][i%15].plot(match, alpha=0.5)
    ev_sess = evs[evs['expstr']==all_vecs[i]['expstr']]
    tsttimes = all_vecs_z[i]['times']['new']-ev_sess['phase1_time'][0]+time_offset
    xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs_z[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='yellow', alpha=0.4)
    tsttimes = all_vecs_z[i]['times']['old']-ev_sess['phase1_time'][0]+time_offset
    xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    all_figsaxs_z[-1][1][i%15].fill_between(
        xs, y1s, y2s, where=where, color='red', alpha=0.4)






































    



    
norm_vecs = {}
for oldnew in ['old', 'new']:
    norm_vecs[oldnew] = all_vecs[0]['target_vectors'][oldnew]/np.sqrt((
        all_vecs[0]['target_vectors'][oldnew]**2).sum())
norm_vecs['diff'] = (all_vecs[0]['target_vectors']['old']-
                     all_vecs[0]['target_vectors']['new'])/np.sqrt(
                         ((all_vecs[0]['target_vectors']['old']-
                           all_vecs[0]['target_vectors']['new'])**2).sum())

av_norm = all_vecs[0]['firing_rate_cont']/np.sqrt((all_vecs[0]['firing_rate_cont']**2).sum(axis=0))
av_norm[~np.isfinite(av_norm)]=0

tst1 = np.dot(norm_vecs['diff'], av_norm)
tst2 = np.dot(norm_vecs['old'], av_norm)
tst3 = np.dot(norm_vecs['new'], av_norm)

plt.plot(tst1[2500:3000], alpha=0.5)
plt.plot(tst2[2500:3000], alpha=0.5)
plt.plot(tst3[2500:3000], alpha=0.5)


import matplotlib.pyplot as plt
plt.ion()
tst1 = np.dot(norm_vecs['diff'], av_norm)
plt.plot(tst1, alpha=0.5)
ev_sess = evs[evs['expstr']==all_vecs[0]['expstr']]
tsttimes = all_vecs[0]['times']['new']-ev_sess['phase1_time'][0]+time_offset
xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
y1s = [1]*len(xs)
where = np.ones(xs.shape, dtype=np.bool)
where[2::3] = 0
#plt.plot(np.round(tsttimes/1000/100)
plt.fill_between(xs, y1s, where=where, color='yellow', alpha=0.4)
tsttimes = all_vecs[0]['times']['old']-ev_sess['phase1_time'][0]+time_offset
xs = np.ravel([(x, x+15, x+15) for x in np.round(tsttimes/1000/100)])
y1s = [1]*len(xs)
where = np.ones(xs.shape, dtype=np.bool)
where[2::3] = 0
#plt.plot(np.round(tsttimes/1000/100)
plt.fill_between(xs, y1s, where=where, color='red', alpha=0.4)





tsttimes = all_vecs[0]['times']['old']-ev_sess['phase1_time'][0]+time_offset
tst1 = [x for x in np.round(tsttimes/1000/100)]

tsttimes = all_vecs[0]['times']['new']-ev_sess['phase1_time'][0]+time_offset
tst2 = [x for x in np.round(tsttimes/1000/100)]







for tv in targvecs_all:
    print(np.round(np.dot(tv['old']/np.sqrt((tv['old']**2).sum()),
                          tv['new']/np.sqrt((tv['new']**2).sum())), 2), '\t',
          np.round(np.sqrt(np.power(tv['old'], 2).sum()),2), '\t',
          np.round(np.sqrt(np.power(tv['new'], 2).sum()),2), '\t',
          np.round(np.sqrt(np.power(tv['old']-tv['new'], 2).sum()),2))


for tv in targvecs_all:
    old_t = tv['old']/(tv['old_std']/np.sqrt(tv['old_n']))
    old_t[~np.isfinite(old_t)] = 0
    new_t = tv['new']/(tv['new_std']/np.sqrt(tv['new_n']))
    new_t[~np.isfinite(new_t)] = 0
    print(np.round(np.dot(old_t/np.sqrt((old_t**2).sum()),
                          new_t/np.sqrt((new_t**2).sum())), 2), '\t',
          np.round(np.sqrt(np.power(old_t, 2).sum()),2), '\t',
          np.round(np.sqrt(np.power(new_t, 2).sum()),2), '\t',
          np.round(np.sqrt(np.power(old_t-new_t, 2).sum()),2))


oldnorms = []
newnorms = []
for tv in targvecs_all:
    oldnorms.append(np.sqrt(np.power(tv['old'], 2).sum()))
    newnorms.append(np.sqrt(np.power(tv['new'], 2).sum()))

oldnorms = np.array(oldnorms)
newnorms = np.array(newnorms)
diffnorm = oldnorms-newnorms















for tv in targvecs_all:
    print(np.round(np.dot(tv['old_norm'], tv['new_norm']), 2), '\t',
          np.round(np.sqrt(np.power(tv['old'], 2).sum()),2), '\t',
          np.round(np.sqrt(np.power(tv['new'], 2).sum()),2), '\t',
          np.round(np.sqrt(np.power(tv['old']-tv['new'], 2).sum()),2))


oldnorms = []
newnorms = []
for tv in targvecs_all:
    oldnorms.append(np.sqrt(np.power(tv['old'], 2).sum()))
    newnorms.append(np.sqrt(np.power(tv['new'], 2).sum()))

oldnorms = np.array(oldnorms)
newnorms = np.array(newnorms)
diffnorm = oldnorms-newnorms


tst1 = np.dot(target_vectors['old'],target_vectors['new'])/(
    np.sqrt((target_vectors['old']**2).sum())*
    np.sqrt((target_vectors['new']**2).sum()))
tst2 = np.dot(target_vectors['old_norm'],target_vectors['new_norm'])


old_vec = firing_rate['old'].mean(1)
new_vec = firing_rate['new'].mean(1)
oldnew_vec = old_vec-new_vec
old_vec_norm = np.sqrt(


np.dot(oldnew_vec, firing_rate_cont)/






























    
expstr = evs['expstr'][0]
ev_sess = evs[evs['expstr']==expstr]
assert len(np.unique(ev_sess['subject']))==1, 'subj len'
subject = ev_sess['subject'][0]
nsefiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                       subject + '/analysis/' + expstr +
                       '/KK/CSC?.Nse'))
nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                            subject + '/analysis/' + expstr +
                            '/KK/CSC??.Nse')))
# not really needed, but just in case:
nsefiles.extend(sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                            subject + '/analysis/' + expstr +
                            '/KK/CSC???.Nse')))

clusterfiles = sorted(glob(rhino_root+'/scratch/josh/BniData/Subjects/' +
                           subject + '/analysis/' + expstr +
                           '/KK/CSC?.clu*'))
clusterfiles.extend(sorted(glob(
    rhino_root+'/scratch/josh/BniData/Subjects/' +
    subject + '/analysis/' + expstr + '/KK/CSC??.clu*')))
# not really needed, but just in case:
clusterfiles.extend(sorted(glob(
    rhino_root+'/scratch/josh/BniData/Subjects/' +
    subject + '/analysis/' + expstr + '/KK/CSC???.clu*')))
if expstr == 's23e15cr':
    # CSC5.clu.1 and CSC5.clu.1 are missing, so we're inserting
    # Nones to make sure clusterfiles and nsefiles align (will
    # skip these later):
    clusterfiles.insert(4, None)
    clusterfiles.insert(4, None)
    # there are no NCS files for channels 41-48 (corresponding to
    # indices 31-39) so we're removing them (starting with the
    # last one to be able to remove a range of indices rather than
    # removing the same index multiple times):
    for i in range(39, 31, -1):
        clusterfiles.pop(i)
        nsefiles.pop(i)
assert len(nsefiles) == len(clusterfiles),\
    'len(nsefiles) == len(clusterfiles)'
assert np.alltrue([int(clusterfile.split('CSC')[-1].split('.')[0]) ==
                   int(nsefile.split('CSC')[-1].split('.')[0])
                   for clusterfile, nsefile in
                   zip(clusterfiles, nsefiles)
                   if clusterfile is not None]), 'compare nses & clusters'
clusterinfo_file = (rhino_root+'/scratch/josh/BniData/Subjects/' +
                    subject + '/analysis/' + expstr +
                    '/clusterInfo.txt')
clusterinfo = None
if os.path.exists(clusterinfo_file):
    clusterinfo = np.loadtxt(clusterinfo_file, delimiter='\t', skiprows=1,
                             dtype={'names': (
                                 'expstr', 'cluster_id',
                                 'hemisphere', 'area', 'quality'),
                                    'formats': ('U16', 'U16', 'U8',
                                                'U8', 'U16')})
    assert len(np.unique(clusterinfo['expstr'])) == 1,\
        'clusterinfo expstr 1'
    assert clusterinfo['expstr'][0] == expstr, 'clusterinfo expstr 2'
    assert np.alltrue(np.unique([
        int(clid.split('cl')[0].split('ch')[1])
        for clid in clusterinfo['cluster_id']]) == np.array([int(nsefile.split(
                'CSC')[-1].split('.')[0]) for nsefile in nsefiles])), 'chan chk'
else:
    raise ValueError('No clusterinfo file!')
    # print('******************************************')
    # print('******************************************')
    # print('******************************************')
    # print('********** NO CLUSTERINFO FILE! **********')
    # print('******************************************')
    # print('******************************************')
    # print('******************************************')
    # os.remove('/home/ctw/data/BNI/data/time_series/'+expstr+'_ts.hdf')
    # continue



sme_filt = ev_sess['first'] &  (ev_sess['paired'] != 'unpaired')
sme_resp = np.empty(np.sum(sme_filt), dtype = 'U3')
sme_bins = {'old': [], 'new': []}
sme_timewindowlength = 1500000  # 1.5 s
for sme_indx, sme_trial in enumerate(ev_sess[sme_filt]):
    # np.where(ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'],
    #          ev_sess[~sme_filt], [None]*np.sum(~sme_filt))
    #indx = np.flatnonzero(
    #    ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'])
    indices = np.flatnonzero(
        ev_sess['stim_word'] == sme_trial['stim_word'])
    assert len(indices) == 2, 'len(indices) == 2'
    old_resp = np.sum(
        ev_sess[indices[1]]['resp_count_old_1']+
        ev_sess[indices[1]]['resp_count_old_2']+
        ev_sess[indices[1]]['resp_count_old_3'])
    new_resp = np.sum(
        ev_sess[indices[1]]['resp_count_new_1']+
        ev_sess[indices[1]]['resp_count_new_2']+
        ev_sess[indices[1]]['resp_count_new_3'])
    if (old_resp > 0) & (new_resp == 0):
        sme_resp[sme_indx] = 'old'
        sme_bins['old'].extend([ev_sess[indices[0]]['phase1_time'],
            ev_sess[indices[0]]['phase1_time'] + sme_timewindowlength])
    elif (old_resp == 0) & (new_resp > 0):
        sme_resp[sme_indx] = 'new'
        sme_bins['new'].extend([ev_sess[indices[0]]['phase1_time'],
            ev_sess[indices[0]]['phase1_time'] + sme_timewindowlength])
    else:
        sme_resp[sme_indx] = '---'

clusterfile = clusterfiles[0]
nsefile = nsefiles[0]

# spikesums = {'old': np.empty(np.sum(sme_resp=='old'), int),
#             'new': np.empty(np.sum(sme_resp=='new'), int)}
spikesums = {'old': [], 'new': []}
spikesums_cont = []
binwidth = 100000  # 100 ms
# binwidth = 200000  # 200 ms
time_offset = 5000000  # 5s
bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
                 ev_sess['phase3_time'][-1]+time_offset, binwidth)
channels = []
clusters = []
clusterIDs = []
for clusterfile, nsefile in zip(
        clusterfiles, nsefiles):
    if clusterfile is None:
        continue
    chan = nsefile.split('.Nse')[0]
    channum = int(chan.split('CSC')[-1])
    # # clusters[channum] = {}
    # cluster[chan] = {}
    nsetimes = np.memmap(nsefile, dtype=nse_dtype, mode='r',
                         offset=header_size)['timestamp']
    numclusts = np.fromfile(clusterfile, dtype=int, sep='\n')[0]
    sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]
    for cluster_num in np.unique(sess_clusts):
        # cluster_num = int(cluster_num)
        clustimes = nsetimes[sess_clusts == cluster_num]
        for oldnew in sme_bins:
            hist, _ = np.histogram(clustimes, sme_bins[oldnew])
            # every other bin is the time between item presentations,
            # so we're skipping these:
            spikesums[oldnew].append(hist[::2])
        channels.append(channum)
        clusters.append(cluster_num)
        clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
        hist, _ = np.histogram(clustimes, bins)
        spikesums_cont.append(hist)
times_cont = bins[:-1] + binwidth/2
spikesums_cont = np.array(spikesums_cont)
spikesums_cont_norm = spikesums_cont/np.atleast_2d(spikesums_cont.max(1)).T
firing_rate_cont = spikesums_cont/(binwidth/1000000)

times = {}
spikesums_norm = {}
firing_rate = {}
for oldnew in spikesums:
    spikesums[oldnew] = np.array(spikesums[oldnew])
    spikesums_norm[oldnew] = spikesums[oldnew]/np.atleast_2d(
        spikesums[oldnew].max(1)).T
    spikesums_norm[oldnew][spikesums[oldnew].max(1)==0] = 0
    firing_rate[oldnew] = spikesums[oldnew]/(sme_timewindowlength/1000000)
    times[oldnew] = np.array(sme_bins[oldnew][::2]) + sme_timewindowlength/2
  
channels = np.array(channels)
cluster = np.array(clusters)
clusterIDs = np.array(clusterIDs)











# binwidth = 50000  # 50 ms
binwidth = 100000  # 100 ms
# binwidth = 200000  # 200 ms
time_offset = 5000000  # 5s
bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
                 ev_sess['phase3_time'][-1]+time_offset, binwidth)
spikesums = []
channels = []
clusters = []
clusterIDs = []
for clusterfile, nsefile in zip(
        clusterfiles, nsefiles):
    if clusterfile is None:
        continue
    chan = nsefile.split('.Nse')[0]
    channum = int(chan.split('CSC')[-1])
    # # clusters[channum] = {}
    # cluster[chan] = {}
    nsetimes = np.memmap(nsefile, dtype=nse_dtype, mode='r',
                         offset=header_size)['timestamp']
    numclusts = np.fromfile(clusterfile, dtype=int, sep='\n')[0]
    sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]
    for cluster_num in np.unique(sess_clusts):
        # cluster_num = int(cluster_num)
        clustimes = nsetimes[sess_clusts == cluster_num]
        hist, _ = np.histogram(clustimes, bins)
        spikesums.append(hist)
        channels.append(channum)
        clusters.append(cluster_num)
        clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
spikesums = np.array(spikesums)
spikesums_norm = spikesums/np.atleast_2d(spikesums.max(1)).T
channels = np.array(channels)
cluster = np.array(clusters)
clusterIDs = np.array(clusterIDs)
times = bins[:-1] + np.diff(bins)/2




clusterfile = clusterfiles[0]
nsefile = nsefiles[0]

# binwidth = 50000  # 50 ms
binwidth = 100000  # 100 ms
# binwidth = 200000  # 200 ms
time_offset = 5000000  # 5s
bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
                 ev_sess['phase3_time'][-1]+time_offset, binwidth)
spikesums = []
channels = []
clusters = []
clusterIDs = []
for clusterfile, nsefile in zip(
        clusterfiles, nsefiles):
    if clusterfile is None:
        continue
    chan = nsefile.split('.Nse')[0]
    channum = int(chan.split('CSC')[-1])
    # # clusters[channum] = {}
    # cluster[chan] = {}
    nsetimes = np.memmap(nsefile, dtype=nse_dtype, mode='r',
                         offset=header_size)['timestamp']
    numclusts = np.fromfile(clusterfile, dtype=int, sep='\n')[0]
    sess_clusts = np.fromfile(clusterfile, dtype=int, sep='\n')[1:]
    for cluster_num in np.unique(sess_clusts):
        # cluster_num = int(cluster_num)
        clustimes = nsetimes[sess_clusts == cluster_num]
        hist, _ = np.histogram(clustimes, bins)
        spikesums.append(hist)
        channels.append(channum)
        clusters.append(cluster_num)
        clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
spikesums = np.array(spikesums)
spikesums_norm = spikesums/np.atleast_2d(spikesums.max(1)).T
channels = np.array(channels)
cluster = np.array(clusters)
clusterIDs = np.array(clusterIDs)
times = bins[:-1] + np.diff(bins)/2

kmeans_predict = KMeans(n_clusters=2).fit_predict(spikesums.T)
kmeans_transform = KMeans(n_clusters=2).fit_transform(spikesums.T)
meanshift_predict = MeanShift(bandwidth=30).fit_predict(spikesums.T)
meanshift_predict45 = MeanShift(bandwidth=45).fit_predict(spikesums.T)
meanshift_predict49 = MeanShift(bandwidth=49).fit_predict(spikesums.T)
meanshift_predict50 = MeanShift(bandwidth=50).fit_predict(spikesums.T)

kmeans_predict_norm = KMeans(n_clusters=2).fit_predict(spikesums_norm.T)
#meanshift_predict_norm = MeanShift(bandwidth=30).fit_predict(spikesums_norm.T)

for cluster_id, samples in itertools.groupby(kmeans_predict_norm):
    print(cluster_id, len(list(samples)))

    if np.unique(sess_clusts)[0] == 1:
        assert numclusts == len(np.unique(sess_clusts)),\
            'numclusts == len(np.unique(sess_clusts))'
    else:
        assert numclusts == (len(np.unique(sess_clusts))+1),\
            'numclusts == (len(np.unique(sess_clusts))+1)'
        clustcounter -= 1
    clustcounter += numclusts
        for cluster_num in np.unique(sess_clusts):
            cluster_num = int(cluster_num)
            if clusterinfo is None:
                cluster[chan][cluster_num] = {
                    'quality': {'clusterinfo': 'N/A'},
                    'hemisphere': 'N/A',
                    'area': 'N/A',
                    'time': None}
            else:
                ci = clusterinfo[clusterinfo[
                    'cluster_id'] == ('ch'+str(channum)+'cl'+str(cluster_num))]
                if len(ci) == 0:
                    cluster[chan][cluster_num] = {
                        'quality': {'clusterinfo': 'N/A'},
                        'hemisphere': 'N/A',
                        'area': 'N/A',
                        'time': None}
                else:
                    assert len(ci) == 1, 'len(ci) == 1'
                    cluster[chan][cluster_num] = {
                        'quality': {'clusterinfo': ci['quality'][0]},
                        'hemisphere': ci['hemisphere'][0],
                        'area': ci['area'][0],
                        'time': None}
