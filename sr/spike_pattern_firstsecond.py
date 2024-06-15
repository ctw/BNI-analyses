
import numpy as np
from glob import glob
import os
#from scipy import signal
#import bni_helpers
#from pathlib import Path

#import h5py
#from numpy.lib import recfunctions
# from ptsa.data.timeseries import TimeSeries

# from sklearn.cluster import KMeans
# from sklearn.cluster import MeanShift
# import itertools
from numpy.random import default_rng
import matplotlib.pyplot as plt


evs = np.load(
    '/home/ctw/Christoph/Analyses/BNI/sr/data/bni_evs_sr.npz')['evs'].view(
        np.recarray)
#evs = recfunctions.append_fields(evs, 'missing_neuro', [False]*len(evs),
#                                 usemask=False, asrecarray=True)

rhino_root = '/home/ctw/fusemounts/rhino'
nse_dtype = [('timestamp', 'uint64'), ('channel', 'uint32'),
             ('unit_id', 'uint32'), ('params', 'uint32', (8,)),
             ('samples', 'uint16', (32,))]
header_size=16*1024



bad_expstr = ['s55e4sr']
#one item presented 4 times, issues with timing

# targvecs_all = []
all_vecs = []
for i, expstr in enumerate(np.unique(evs['expstr'])):
    if expstr in bad_expstr:
        continue
    print(i+1, '/', len(np.unique(evs['expstr'])), expstr)
    ev_sess = evs[evs['expstr']==expstr]
    assert len(np.unique(ev_sess['subject'])) == 1, 'subj len'
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
    events = {}
    first_filt = ev_sess['first'] & (ev_sess['paired'] != 'unpaired')
    events['first'] = ev_sess[first_filt]
    # sme_resp = np.empty(np.sum(sme_filt), dtype = 'U3')
    # sme_bins = {'old': [], 'new': []}
    bins = {'first': [], 'second': []}
    second_indices = np.zeros(first_filt.sum(), dtype=np.int32)-1
    # timewindowlength = 1500000  # 1.5 s
    timewindowlength = 1000000  # 1 s
    # prev_time = 0
    for indx, trial in enumerate(ev_sess[first_filt]):
        # np.where(ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'],
        #          ev_sess[~sme_filt], [None]*np.sum(~sme_filt))
        #indx = np.flatnonzero(
        #    ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'])
        indices = np.flatnonzero(
            ev_sess['stim_word'] == trial['stim_word'])
        # if ev_sess[indices[0]]['phase1_time'] < (prev_time+timewindowlength):
        #     print('STOP')
        #     break
        # prev_time = ev_sess[indices[0]]['phase1_time']
        
        # # s55e4sr has 4 presentations of 'nose':
        # assert ((len(indices) == 2) or
        #         (len(indices) == 3) or
        #         (len(indices) == 4)), 'len(indices) == 2 or 3 or 4'
        # s55e4sr has 4 presentations of 'nose':
        assert ((len(indices) == 2) or
                (len(indices) == 3)), 'len(indices) == 2 or 3'
        second_indices[indx] = indices[1]
        bins['first'].extend([ev_sess[indices[0]]['phase1_time'],
                              ev_sess[indices[0]]['phase1_time'] +
                              timewindowlength])
        bins['second'].extend([ev_sess[indices[1]]['phase1_time'],
                               ev_sess[indices[1]]['phase1_time'] +
                               timewindowlength])
        
    assert np.alltrue(second_indices>0), 'second_indices'
    events['second'] = ev_sess[second_indices]
    # spikesums = {'old': np.empty(np.sum(sme_resp=='old'), int),
    #             'new': np.empty(np.sum(sme_resp=='new'), int)}
    spikesums = {}
    # spikesums_cont = []
    #binwidth = 100000  # 100 ms
    # binwidth = 200000  # 200 ms
    # binwidth = 100000  # 100 ms
    # binwidth = 500000  # 500 ms
    # time_offset = 5000000  # 5s
    # bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
    #                  ev_sess['phase3_time'][-1]+time_offset, binwidth)
    channels = {}
    clusters = {}
    clusterIDs = {}
    times = {}
    spikesums_norm = {}
    firing_rate = {}
    all_areas = ['H', 'A', 'AC', 'PF']
    for area in all_areas:
        channels[area] = []
        clusters[area] = []
        clusterIDs[area] = []
        spikesums[area] = {'first': [], 'second': []}
        times[area] = {'first': [], 'second': []}
        spikesums_norm[area] = {'first': [], 'second': []}
        firing_rate[area] = {'first': [], 'second': []}
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
            if ci.size == 0:
                continue
            if ci['quality'] == 'NOISE':
                #print(channum, cluster_num)
                continue
            if ci['area'][0] not in channels:
                print(ci['area'][0], 'not in areas!')
                continue
            # cluster_num = int(cluster_num)
            clustimes = nsetimes[sess_clusts == cluster_num]
            hist, _ = np.histogram(clustimes, bins['first'])
            # every other bin is the time between item presentations,
            # so we're skipping these:
            spikesums[ci['area'][0]]['first'].append(hist[::2])
            # bins in np.histogram need to monotonically increase, but
            # the indices for the second presentations do not, so we
            # need to deal with this manually:
            hist, _ = np.histogram(clustimes, sorted(bins['second']))
            indxsort = np.argsort(bins['second'])
            spikesums[ci['area'][0]]['second'].append(hist[indxsort[::2]])
            channels[ci['area'][0]].append(channum)
            clusters[ci['area'][0]].append(cluster_num)
            clusterIDs[ci['area'][0]].append(
                'CSC'+str(channum)+'cl'+str(cluster_num))
            # hist, _ = np.histogram(clustimes, bins)
            # spikesums_cont.append(hist)
    #times_cont = bins[:-1] + binwidth/2
    #spikesums_cont = np.array(spikesums_cont)
    #spikesums_cont_norm = spikesums_cont/np.atleast_2d(spikesums_cont.max(1)).T
    #firing_rate_cont = spikesums_cont/(binwidth/1000000)
    for area in channels:
        for firstsecond in spikesums[area]:
            spikesums[area][firstsecond] = np.array(
                spikesums[area][firstsecond])
            if len(spikesums[area][firstsecond].shape) == 1:
                spikesums_norm[area][firstsecond] = spikesums[area][firstsecond]
            else:
                spikesums_norm[area][firstsecond] = spikesums[area][
                    firstsecond]/np.atleast_2d(
                        spikesums[area][firstsecond].max(1)).T
                spikesums_norm[area][firstsecond][spikesums[area][
                    firstsecond].max(1)==0] = 0
            firing_rate[area][firstsecond] = spikesums[area][
                firstsecond]/(timewindowlength/1000000)
            times[area][firstsecond] = np.array(
                bins[firstsecond][::2]) + timewindowlength/2
        channels[area] = np.array(channels[area])
        clusters[area] = np.array(clusters[area])
        clusterIDs[area] = np.array(clusterIDs[area])
    all_vecs.append({'expstr': expstr, 'times': times, 'clusters': clusters,
                     'clusterIDs': clusterIDs, 'channels': channels,
                     'spikesums_norm': spikesums_norm, 'spikesums': spikesums,
                     'firing_rate': firing_rate, 'events': events})


expstrs = np.unique([av['expstr'] for av in all_vecs])
subjects = np.unique([av['events']['first']['subject'][0] if (
    np.alltrue(av['events']['first']['subject']==av[
        'events']['first']['subject'][0]) and
    np.alltrue(av['events']['second']['subject']==av[
        'events']['first']['subject'][0])) else np.nan
                      for av in all_vecs])
assert np.alltrue([s is not np.nan for s in subjects])

clusterIDs = np.concatenate(
    [av['clusterIDs'][area]
     for av in all_vecs for area in av['clusterIDs']])
# cluster IDs repeat across subjects, so duplicate IDs are expected

# len(expstrs) == 17
# len(subjects)  == 16
# len(clusterIDs) == 1025










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
#     events = {}
#     first_filt = ev_sess['first'] & (ev_sess['paired'] != 'unpaired')
#     events['first'] = ev_sess[first_filt]
#     # sme_resp = np.empty(np.sum(sme_filt), dtype = 'U3')
#     # sme_bins = {'old': [], 'new': []}
#     bins = {'first': [], 'second': []}
#     second_indices = np.zeros(first_filt.sum(), dtype=np.int32)-1
#     timewindowlength = 1500000  # 1.5 s
#     for indx, trial in enumerate(ev_sess[first_filt]):
#         # np.where(ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'],
#         #          ev_sess[~sme_filt], [None]*np.sum(~sme_filt))
#         #indx = np.flatnonzero(
#         #    ev_sess[~sme_filt]['stim_word'] == sme_trial['stim_word'])
#         indices = np.flatnonzero(
#             ev_sess['stim_word'] == trial['stim_word'])
#         assert len(indices) == 2, 'len(indices) == 2'
#         second_indices[indx] = indices[1]
#         bins['first'].extend([ev_sess[indices[0]]['phase1_time'],
#                               ev_sess[indices[0]]['phase1_time'] +
#                               timewindowlength])
#         bins['second'].extend([ev_sess[indices[1]]['phase1_time'],
#                                ev_sess[indices[1]]['phase1_time'] +
#                                timewindowlength])
        
#     assert np.alltrue(second_indices>0), 'second_indices'
#     events['second'] = ev_sess[second_indices]
#     # spikesums = {'old': np.empty(np.sum(sme_resp=='old'), int),
#     #             'new': np.empty(np.sum(sme_resp=='new'), int)}
#     spikesums = {'first': [], 'second': []}
#     # spikesums_cont = []
#     #binwidth = 100000  # 100 ms
#     # binwidth = 200000  # 200 ms
#     # binwidth = 100000  # 100 ms
#     # binwidth = 500000  # 500 ms
#     # time_offset = 5000000  # 5s
#     # bins = np.arange(ev_sess['phase1_time'][0]-time_offset,
#     #                  ev_sess['phase3_time'][-1]+time_offset, binwidth)
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
#             hist, _ = np.histogram(clustimes, bins['first'])
#             # every other bin is the time between item presentations,
#             # so we're skipping these:
#             spikesums['first'].append(hist[::2])
#             # bins in np.histogram need to monotonically increase, but
#             # the indices for the second presentations do not, so we
#             # need to deal with this manually:
#             hist, _ = np.histogram(clustimes, sorted(bins['second']))
#             indxsort = np.argsort(bins['second'])
#             spikesums['second'].append(hist[indxsort[::2]])
#             channels.append(channum)
#             clusters.append(cluster_num)
#             clusterIDs.append('CSC'+str(channum)+'cl'+str(cluster_num))
#             # hist, _ = np.histogram(clustimes, bins)
#             # spikesums_cont.append(hist)
#     #times_cont = bins[:-1] + binwidth/2
#     #spikesums_cont = np.array(spikesums_cont)
#     #spikesums_cont_norm = spikesums_cont/np.atleast_2d(spikesums_cont.max(1)).T
#     #firing_rate_cont = spikesums_cont/(binwidth/1000000)
#     times = {}
#     spikesums_norm = {}
#     firing_rate = {}
#     for firstsecond in spikesums:
#         spikesums[firstsecond] = np.array(spikesums[firstsecond])
#         spikesums_norm[firstsecond] = spikesums[firstsecond]/np.atleast_2d(
#             spikesums[firstsecond].max(1)).T
#         spikesums_norm[firstsecond][spikesums[firstsecond].max(1)==0] = 0
#         firing_rate[firstsecond] = spikesums[firstsecond]/(timewindowlength/
#                                                            1000000)
#         times[firstsecond] = np.array(
#             bins[firstsecond][::2]) + timewindowlength/2
#     channels = np.array(channels)
#     cluster = np.array(clusters)
#     clusterIDs = np.array(clusterIDs)
#     all_vecs.append({'expstr': expstr, 'times': times, 'clusters': clusters,
#                      'clusterIDs': clusterIDs, 'channels': channels,
#                      'spikesums_norm': spikesums_norm, 'spikesums': spikesums,
#                      'firing_rate': firing_rate, 'events': events})





rng = default_rng()
# permutations = 100000
permutations = 10000
all_ps = {}
for area in all_areas:
    all_ps[area] = {'pval': [], 'pmean': [], 'pstd': [], 'sim_mean': [],
                    'hit': [], 'fa': []}
    for i, av in enumerate(all_vecs):
        print(area, i+1, '/', len(all_vecs), av['expstr'])
        if len(av['firing_rate'][area]['first'].shape) != 2:
            continue
        norm_vecs = {}
        for firstsecond in av['firing_rate'][area]:
            norm_vecs[firstsecond] = av['firing_rate'][area][
                firstsecond]/np.sqrt((av['firing_rate'][area][
                    firstsecond]**2).sum(0))
        
        sim_mean = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
            norm_vecs['first'].T, norm_vecs['second'])))))
        sim_perm_mean = np.ones(permutations)*np.nan
        for p in range(permutations):
            second_perm = rng.permutation(norm_vecs['second'], 1)
            sim_perm_mean[p] = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
                norm_vecs['first'].T, second_perm)))))
        # sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
        # sim_perm_mean = np.nanmean(sim_perm, 1)
        all_ps[area]['pval'].append(np.sum(sim_mean<sim_perm_mean)/permutations)
        all_ps[area]['pmean'].append(sim_perm_mean.mean())
        all_ps[area]['pstd'].append(sim_perm_mean.std())
        all_ps[area]['sim_mean'].append(sim_mean)
        hitfilt = (((av['events']['second']['resp_count_old_1']+
                     av['events']['second']['resp_count_old_2']+
                     av['events']['second']['resp_count_old_3']) > 0) &
                   ((av['events']['second']['resp_count_new_1']+
                     av['events']['second']['resp_count_new_2']+
                     av['events']['second']['resp_count_new_3']) == 0) &
                   ((av['events']['second']['resp_count_other_1']+
                     av['events']['second']['resp_count_other_2']+
                     av['events']['second']['resp_count_other_3']) == 0))
        fafilt = (((av['events']['first']['resp_count_old_1']+
                    av['events']['first']['resp_count_old_2']+
                    av['events']['first']['resp_count_old_3']) > 0) &
                  ((av['events']['first']['resp_count_new_1']+
                    av['events']['first']['resp_count_new_2']+
                    av['events']['first']['resp_count_new_3']) == 0) &
                  ((av['events']['first']['resp_count_other_1']+
                    av['events']['first']['resp_count_other_2']+
                    av['events']['first']['resp_count_other_3']) == 0))
        all_ps[area]['hit'].append(hitfilt.mean())
        all_ps[area]['fa'].append(fafilt.mean())

# area = 'all'
# for area in all_areas:
all_ps['all'] = {'pval': [], 'pmean': [], 'pstd': [], 'sim_mean': [],
                 'hit': [], 'fa': []}
for i, av in enumerate(all_vecs):
    print(i+1, '/', len(all_vecs), av['expstr'])
    norm_vecs = {}
    for firstsecond in av['firing_rate'][all_areas[0]]:
        all_fr = np.concatenate([av['firing_rate'][area][
            firstsecond] for area in all_areas if len(av['firing_rate'][area][
                firstsecond].shape)==2])
        norm_vecs[firstsecond] = all_fr/np.sqrt((all_fr**2).sum(0))

    sim_mean = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
        norm_vecs['first'].T, norm_vecs['second'])))))
    sim_perm_mean = np.ones(permutations)*np.nan
    for p in range(permutations):
        second_perm = rng.permutation(norm_vecs['second'], 1)
        sim_perm_mean[p] = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
            norm_vecs['first'].T, second_perm)))))
    # sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
    # sim_perm_mean = np.nanmean(sim_perm, 1)
    all_ps['all']['pval'].append(np.sum(sim_mean<sim_perm_mean)/permutations)
    all_ps['all']['pmean'].append(sim_perm_mean.mean())
    all_ps['all']['pstd'].append(sim_perm_mean.std())
    all_ps['all']['sim_mean'].append(sim_mean)
    hitfilt = (((av['events']['second']['resp_count_old_1']+
                 av['events']['second']['resp_count_old_2']+
                 av['events']['second']['resp_count_old_3']) > 0) &
               ((av['events']['second']['resp_count_new_1']+
                 av['events']['second']['resp_count_new_2']+
                 av['events']['second']['resp_count_new_3']) == 0) &
               ((av['events']['second']['resp_count_other_1']+
                 av['events']['second']['resp_count_other_2']+
                 av['events']['second']['resp_count_other_3']) == 0))
    fafilt = (((av['events']['first']['resp_count_old_1']+
                av['events']['first']['resp_count_old_2']+
                av['events']['first']['resp_count_old_3']) > 0) &
              ((av['events']['first']['resp_count_new_1']+
                av['events']['first']['resp_count_new_2']+
                av['events']['first']['resp_count_new_3']) == 0) &
              ((av['events']['first']['resp_count_other_1']+
                av['events']['first']['resp_count_other_2']+
                av['events']['first']['resp_count_other_3']) == 0))
    all_ps['all']['hit'].append(hitfilt.mean())
    all_ps['all']['fa'].append(fafilt.mean())

# import pickle
pickle.dump(all_ps, open('./data/all_ps_zmean.pickle', 'xb'), -1)



ds = ((np.array(all_ps['all']['sim_mean'])-np.array(all_ps['all']['pmean']))/
      np.array(all_ps['all']['pstd']))
filt = np.isfinite(ds)
np.corrcoef(ds[filt],(np.array(all_ps['all']['hit'])[filt]-
                      np.array(all_ps['all']['fa'])[filt]))
    


fig, axs = plt.subplots(1, 1, figsize=(3, 3))
hist_range = (-1.5, 11)
axs.hist(ds, range=hist_range)
axs.set_xlabel(r'Neural similarity ($d^\prime$)')
axs.set_ylabel('count')
#axs[0].set_title(r'RVL$_{new}$-RVL$_{old}$')
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.999,
                    wspace=0.08, hspace=0.02)
#fig.savefig('./figs/similarity_ds_nocol.pdf')
axs.hist(ds[np.array(all_ps['all']['pval'])<0.05],
            color='yellow', alpha=0.5, range=hist_range)
#fig.savefig('./figs/similarity_ds.pdf')



fig, axs = plt.subplots(1, 1, figsize=(4, 3))
axs.plot(ds[filt],(np.array(all_ps['all']['hit'])[filt]-
                      np.array(all_ps['all']['fa'])[filt]), 'ko')
fig.subplots_adjust(left=0.19, bottom=0.15, right=0.99, top=0.999,
                    wspace=0.08, hspace=0.02)
axs.set_xlabel(r'Neural similarity ($d^\prime$)')
axs.set_ylabel(r'Recognition performance (hit $-$ FA)')
fig.savefig('./figs/rec_ds_scatterplot.pdf')



















rng = default_rng()
# permutations = 100000
permutations = 10000
all_ps = {'pval': [], 'pmean': [], 'pstd': [], 'sim_mean': [],
          'hit': [], 'fa': []}
for i, av in enumerate(all_vecs):
    print(i+1, '/', len(all_vecs), av['expstr'])
    norm_vecs = {}
    for firstsecond in av['firing_rate']:
        norm_vecs[firstsecond] = av['firing_rate'][firstsecond]/np.sqrt((
            av['firing_rate'][firstsecond]**2).sum(0))
    sim_mean = np.nanmean(np.diag(np.dot(norm_vecs['first'].T,
                                         norm_vecs['second'])))
    sim_perm_mean = np.ones(permutations)*np.nan
    for p in range(permutations):
        second_perm = rng.permutation(norm_vecs['second'], 1)
        sim_perm_mean[p]=np.nanmean(np.diag(np.dot(norm_vecs['first'].T,
                                                   second_perm)))
    # sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
    # sim_perm_mean = np.nanmean(sim_perm, 1)
    all_ps['pval'].append(np.sum(sim_mean<sim_perm_mean)/permutations)
    all_ps['pmean'].append(sim_perm_mean.mean())
    all_ps['pstd'].append(sim_perm_mean.std())
    all_ps['sim_mean'].append(sim_mean)
    hitfilt = (((av['events']['second']['resp_count_old_1']+
                 av['events']['second']['resp_count_old_2']+
                 av['events']['second']['resp_count_old_3']) > 0) &
               ((av['events']['second']['resp_count_new_1']+
                 av['events']['second']['resp_count_new_2']+
                 av['events']['second']['resp_count_new_3']) == 0) &
               ((av['events']['second']['resp_count_other_1']+
                 av['events']['second']['resp_count_other_2']+
                 av['events']['second']['resp_count_other_3']) == 0))
    fafilt = (((av['events']['first']['resp_count_old_1']+
                av['events']['first']['resp_count_old_2']+
                av['events']['first']['resp_count_old_3']) > 0) &
              ((av['events']['first']['resp_count_new_1']+
                av['events']['first']['resp_count_new_2']+
                av['events']['first']['resp_count_new_3']) == 0) &
              ((av['events']['first']['resp_count_other_1']+
                av['events']['first']['resp_count_other_2']+
                av['events']['first']['resp_count_other_3']) == 0))
    all_ps['hit'].append(hitfilt.mean())
    all_ps['fa'].append(fafilt.mean())


ds = ((np.array(all_ps['sim_mean'])-np.array(all_ps['pmean']))/
      np.array(all_ps['pstd']))
filt = np.isfinite(ds)
np.corrcoef(ds[filt],(np.array(all_ps['hit'])[filt]-
                      np.array(all_ps['fa'])[filt]))
    
    
rng = default_rng()
# permutations = 100000
permutations = 10000
all_ps_lag = {'pval': [], 'pmean': [], 'pstd': [], 'sim_mean': [],
              'hit': [], 'fa': []}
for i, av in enumerate(all_vecs):
    print(i+1, '/', len(all_vecs), av['expstr'])
    norm_vecs = {}
    sim_mean = {}
    sim_perm_mean = {}
    pval = {}
    pmean = {}
    pstd = {}
    hit = {}
    fafilt = (((av['events']['first']['resp_count_old_1']+
                av['events']['first']['resp_count_old_2']+
                av['events']['first']['resp_count_old_3']) > 0) &
              ((av['events']['first']['resp_count_new_1']+
                av['events']['first']['resp_count_new_2']+
                av['events']['first']['resp_count_new_3']) == 0) &
              ((av['events']['first']['resp_count_other_1']+
                av['events']['first']['resp_count_other_2']+
                av['events']['first']['resp_count_other_3']) == 0))
    for lag in np.unique(av['events']['first']['lag']):
        lagfilt = av['events']['first']['lag']==lag
        hitfilt = (
            ((av['events']['second'][lagfilt]['resp_count_old_1']+
              av['events']['second'][lagfilt]['resp_count_old_2']+
              av['events']['second'][lagfilt]['resp_count_old_3']) > 0) &
            ((av['events']['second'][lagfilt]['resp_count_new_1']+
              av['events']['second'][lagfilt]['resp_count_new_2']+
              av['events']['second'][lagfilt]['resp_count_new_3']) == 0) &
            ((av['events']['second'][lagfilt]['resp_count_other_1']+
              av['events']['second'][lagfilt]['resp_count_other_2']+
              av['events']['second'][lagfilt]['resp_count_other_3']) == 0))
        hit[lag] = hitfilt.mean()
        norm_vecs[lag] = {}
        for firstsecond in av['firing_rate']:
            norm_vecs[lag][firstsecond] = (
                av['firing_rate'][firstsecond][:, lagfilt]/np.sqrt((
                    av['firing_rate'][firstsecond][:, lagfilt]**2).sum(0)))
        sim_mean[lag] = np.nanmean(np.diag(np.dot(norm_vecs[lag]['first'].T,
                                                  norm_vecs[lag]['second'])))
        sim_perm_mean[lag] = np.ones(permutations)*np.nan
        for p in range(permutations):
            second_perm = rng.permutation(norm_vecs[lag]['second'], 1)
            sim_perm_mean[lag][p]=np.nanmean(
                np.diag(np.dot(norm_vecs[lag]['first'].T, second_perm)))
        pval[lag] = np.sum(sim_mean[lag]<sim_perm_mean[lag])/permutations
        pmean[lag] = sim_perm_mean[lag].mean()
        pstd[lag] = sim_perm_mean[lag].std()
    # sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
    # sim_perm_mean = np.nanmean(sim_perm, 1)
    all_ps_lag['pval'].append(pval)
    all_ps_lag['pmean'].append(pmean)
    all_ps_lag['pstd'].append(pstd)
    all_ps_lag['sim_mean'].append(sim_mean)
    all_ps_lag['hit'].append(hit)
    all_ps_lag['fa'].append(fafilt.mean())








rng = default_rng()
# permutations = 100000
permutations = 10000
all_ps_samediff = {'pval': [], 'pmean': [], 'pstd': [], 'sim_mean': [],
                   'hit': [], 'fa': []}
for i, av in enumerate(all_vecs):
    print(i+1, '/', len(all_vecs), av['expstr'])
    norm_vecs = {}
    sim_mean = {}
    sim_perm_mean = {}
    pval = {}
    pmean = {}
    pstd = {}
    hit = {}
    fafilt = (((av['events']['first']['resp_count_old_1']+
                av['events']['first']['resp_count_old_2']+
                av['events']['first']['resp_count_old_3']) > 0) &
              ((av['events']['first']['resp_count_new_1']+
                av['events']['first']['resp_count_new_2']+
                av['events']['first']['resp_count_new_3']) == 0) &
              ((av['events']['first']['resp_count_other_1']+
                av['events']['first']['resp_count_other_2']+
                av['events']['first']['resp_count_other_3']) == 0))

    for samediff in np.unique(av['events']['first']['paired']):
        samedifffilt = av['events']['second']['paired']==samediff
        hitfilt = (
            ((av['events']['second'][samedifffilt]['resp_count_old_1']+
              av['events']['second'][samedifffilt]['resp_count_old_2']+
              av['events']['second'][samedifffilt][
                  'resp_count_old_3']) > 0) &
            ((av['events']['second'][samedifffilt]['resp_count_new_1']+
              av['events']['second'][samedifffilt]['resp_count_new_2']+
              av['events']['second'][samedifffilt][
                  'resp_count_new_3']) == 0) &
            ((av['events']['second'][samedifffilt]['resp_count_other_1']+
              av['events']['second'][samedifffilt]['resp_count_other_2']+
              av['events']['second'][samedifffilt][
                  'resp_count_other_3']) == 0))
        hit[samediff] = hitfilt.mean()
        norm_vecs[samediff] = {}
        for firstsecond in av['firing_rate']:
            samedifffilt = av['events'][firstsecond]['paired']==samediff
            norm_vecs[samediff][firstsecond] = (
                av['firing_rate'][firstsecond][:, samedifffilt]/np.sqrt((
                    av['firing_rate'][firstsecond][:, samedifffilt]**2).sum(0)))
        sim_mean[samediff] = np.nanmean(
            np.diag(np.dot(norm_vecs[samediff]['first'].T,
                           norm_vecs[samediff]['second'])))
        sim_perm_mean[samediff] = np.ones(permutations)*np.nan
        for p in range(permutations):
            second_perm = rng.permutation(norm_vecs[samediff]['second'], 1)
            sim_perm_mean[samediff][p]=np.nanmean(
                np.diag(np.dot(norm_vecs[samediff]['first'].T, second_perm)))
        pval[samediff] = np.sum(
            sim_mean[samediff]<sim_perm_mean[samediff])/permutations
        pmean[samediff] = sim_perm_mean[samediff].mean()
        pstd[samediff] = sim_perm_mean[samediff].std()
    # sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
    # sim_perm_mean = np.nanmean(sim_perm, 1)
    all_ps_samediff['pval'].append(pval)
    all_ps_samediff['pmean'].append(pmean)
    all_ps_samediff['pstd'].append(pstd)
    all_ps_samediff['sim_mean'].append(sim_mean)
    all_ps_samediff['hit'].append(hit)
    all_ps_samediff['fa'].append(fafilt.mean())








ds = []
hits = []
fas = []
for i, av in enumerate(all_vecs):
    ds.append([((all_ps_lag['sim_mean'][i][lag]-all_ps_lag['pmean'][i][lag])/
                all_ps_lag['pstd'][i][lag]) for lag in np.unique(
                    av['events']['first']['lag'])])
    hits.append([all_ps_lag['hit'][i][lag] for lag in np.unique(
        av['events']['first']['lag'])])
    fas.append([all_ps_lag['fa'][i] for lag in np.unique(
        av['events']['first']['lag'])])

dsr = np.ravel(ds)
filt = np.isfinite(dsr)
hitsr = np.ravel(hits)
fasr = np.ravel(fas)
np.corrcoef(dsr[filt], hitsr[filt]-fasr[filt])







ds = []
hits = []
fas = []
for i, av in enumerate(all_vecs):
    ds.append([((all_ps_samediff['sim_mean'][i][samediff]-
                 all_ps_samediff['pmean'][i][samediff])/
                all_ps_samediff['pstd'][i][samediff])
               for samediff in ['same', 'different']])
    hits.append([all_ps_samediff['hit'][i][samediff]
                 for samediff in ['same', 'different']])
    fas.append([all_ps_samediff['fa'][i]
                for samediff in ['same', 'different']])

ds = np.array(ds)
diff = ds[:,0]-ds[:,1]


dsr = np.ravel(ds)
filt = np.isfinite(dsr)
hitsr = np.ravel(hits)
fasr = np.ravel(fas)
np.corrcoef(dsr[filt], hitsr[filt]-fasr[filt])








import matplotlib.pyplot as plt
plt.ion()
from scipy import stats
fig, axs = plt.subplots(8, 6, figsize=(12,9), sharex=True, sharey=True)
rowindx = 0
colindx = -1
for i, av in enumerate(all_vecs):
    ds = [((all_ps['sim_mean'][i][lag]-all_ps['pmean'][i][lag])/
           all_ps['pstd'][i][lag]) for lag in np.unique(
               av['events']['first']['lag'])]
    ds2 = [stats.norm.ppf(1-all_ps['pval'][i][lag])
           if all_ps['pval'][i][lag]>0 else
           stats.norm.ppf(1-(1/permutations)*0.5)
           for lag in np.unique(av['events']['first']['lag'])]
    rowindx = i%8
    if rowindx == 0:
        colindx += 1
    axs[rowindx, colindx].plot(np.unique(av['events']['first']['lag']),
                               ds, 'ko-', alpha=0.5)
    axs[rowindx, colindx].plot(np.unique(av['events']['first']['lag']),
                               ds2, 'rs-', alpha=0.5)



fig, axs = plt.subplots(8, 6, figsize=(12,9), sharex=True, sharey=True)
rowindx = 0
colindx = -1
for i, av in enumerate(all_vecs):
    ds = [((all_ps_samediff['sim_mean'][i][samediff]-
            all_ps_samediff['pmean'][i][samediff])/
           all_ps_samediff['pstd'][i][samediff])
          for samediff in ['same', 'different']]
    rowindx = i%8
    if rowindx == 0:
        colindx += 1
    axs[rowindx, colindx].plot(ds, 'ko-', alpha=0.5)

    










rng = default_rng()
# permutations = 100000
permutations = 10000
all_ps = {'pval': [], 'pmean': [], 'pstd': [], 'sim_mean': []}
for i, av in enumerate(all_vecs):
    print(i+1, '/', len(all_vecs), av['expstr'])
    norm_vecs = {}
    sim_mean = {}
    sim_perm_mean = {}
    pval = {}
    pmean = {}
    pstd = {}
    for lag in np.unique(av['events']['first']['lag']):
        norm_vecs[lag] = {}
        for firstsecond in av['firing_rate']:
            lagfilt = av['events'][firstsecond]['lag']==lag
            norm_vecs[lag][firstsecond] = (
                av['firing_rate'][firstsecond][:, lagfilt]/np.sqrt((
                    av['firing_rate'][firstsecond][:, lagfilt]**2).sum(0)))
        sim_mean[lag] = np.nanmean(np.diag(np.dot(norm_vecs[lag]['first'].T,
                                                  norm_vecs[lag]['second'])))
        sim_perm_mean[lag] = np.ones(permutations)*np.nan
        for p in range(permutations):
            second_perm = rng.permutation(norm_vecs[lag]['second'], 1)
            sim_perm_mean[lag][p]=np.nanmean(
                np.diag(np.dot(norm_vecs[lag]['first'].T, second_perm)))
        pval[lag] = np.sum(sim_mean[lag]<sim_perm_mean[lag])/permutations
        pmean[lag] = sim_perm_mean[lag].mean()
        pstd[lag] = sim_perm_mean[lag].std()
    # sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
    # sim_perm_mean = np.nanmean(sim_perm, 1)
    all_ps['pval'].append(pval)
    all_ps['pmean'].append(pmean)
    all_ps['pstd'].append(pstd)
    all_ps['sim_mean'].append(sim_mean)
    



all_ps = np.array(all_ps)

i=0
norm_vecs = {}
for firstsecond in all_vecs[i]['firing_rate']:
    norm_vecs[firstsecond] = all_vecs[i]['firing_rate'][firstsecond]/np.sqrt((
        all_vecs[i]['firing_rate'][firstsecond]**2).sum(0))
    # norm_vecs[firstsecond][~np.isfinite(norm_vecs[firstsecond])] = 0

sim = np.diag(np.dot(norm_vecs['first'].T, norm_vecs['second']))
rng = default_rng()
sim_perm = []
permutations = 10000
for p in range(permutations):
    second_perm = rng.permutation(norm_vecs['second'], 1)
    sim_perm.append(np.diag(np.dot(norm_vecs['first'].T, second_perm)))

sim_perm_mean = np.sort(np.nanmean(sim_perm, 1))
    

    
first_norm = all_vecs[i]['firing_rate']['first']/np.sqrt((
    all_vecs[i]['firing_rate']['first']**2).sum(0))
first_norm[~np.isfinite(first_norm)] = 0
first_norm = all_vecs[i]['firing_rate']['first']/np.sqrt((
    all_vecs[i]['firing_rate']['first']**2).sum(0))
first_norm[~np.isfinite(first_norm)] = 0

    
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
