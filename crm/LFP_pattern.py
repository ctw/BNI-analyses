
import numpy as np
from glob import glob
import os
#from scipy import signal
#import bni_helpers
#from pathlib import Path

import h5py
#from numpy.lib import recfunctions
from ptsa.data.timeseries import TimeSeries
from ptsa.data.filters import MorletWaveletFilter
import matplotlib.pyplot as plt

#from sklearn.cluster import KMeans
# from sklearn.cluster import MeanShift
#import itertools
import h5py
import pandas as pd
import xarray

# evs = np.load('/home/ctw/data/BNI/data/bni_evs_crm.npz')['evs'].view(
#     np.recarray)

rhino_root = '/home/ctw/fusemounts/rhino'
freqs = np.geomspace(2, 200, 15)
# for i, expstr in enumerate(np.unique(evs['expstr'])):
#     powfile = (
# rhino_root+'/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
#                expstr+'_tspow.hdf')
#     if os.path.exists(powfile):
#         continue
#     print(i+1, '/', len(np.unique(evs['expstr'])), expstr)
#     ts = TimeSeries.from_hdf(
#         rhino_root+'/scratch/josh/BNI_processed/time_series/'+expstr+'_ts.hdf')
#     ts_pow = MorletWaveletFilter(ts, freqs, output='power').filter()
#     ts_pow[:] = np.log(ts_pow)
#     ts_pow = ts_pow.resampled(50)
#     ts_pow = ts_pow.remove_buffer(1)
#     ts_pow.to_hdf(powfile)


evs = np.load('/home/ctw/data/BNI/data/bni_evs_crm.npz')['evs'].view(
    np.recarray)
expstr = np.unique(evs['expstr'])[3]
ncsfile = rhino_root+'/scratch/josh/BNI_processed/ncs/'+expstr+'_ncs.hdf'

for i, expstr in enumerate(np.unique(evs['expstr'])):
    powfile = (rhino_root+'/scratch/josh/BNI_processed/time_series_pow_all/'+
               expstr+'_tspow.hdf')
    if os.path.exists(powfile):
        continue
    h5py.File(powfile, 'w')
    print(i+1, '/', len(np.unique(evs['expstr'])), expstr)
    ncsfile = rhino_root+'/scratch/josh/BNI_processed/ncs/'+expstr+'_ncs.hdf'
    ncs = h5py.File(ncsfile, 'r')
    ncs_ts = TimeSeries(
        ncs['ncs']['data'][:],
        coords=dict(samplerate=1000, channel=ncs['ncs']['channels'],
                    time=ncs['ncs']['times'][:]/1000000),
        dims=['channel', 'time'])
    # time=(ncs['ncs']['times'][:]-ncs['ncs']['times'][0])/
    # 1000000), dims=['channel', 'time'])

    ts_pow = MorletWaveletFilter(ncs_ts, freqs, output='power').filter()
    ts_pow[:] = np.log(ts_pow)
    ts_pow = ts_pow.resampled(50)
    ts_pow = ts_pow.remove_buffer(1)
    ts_pow.to_hdf(powfile)





tspowfiles = glob(
    rhino_root+'/scratch/josh/BNI_processed/time_series_pow_wordevs/s*e*cr_tspow.hdf')


for i, tspowfile in enumerate(tspowfiles):
    # this file is full of NaNs:
    if ((tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                     's23e15cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's23e5cr_tspow.hdf'))):
        indx_adjust +=1
        continue
    print(tspowfile)
    ts_pow_wordevs = TimeSeries.from_hdf(tspowfile)
    # remove missing data (selecting based on missing_neuro is not enough)
    ts_pow_wordevs = ts_pow_wordevs.isel(
        event=((~ts_pow_wordevs['event'].values['missing_neuro']) &
               np.isfinite(ts_pow_wordevs.mean('frequency', skipna=False).mean(
                   'channel', skipna=False).mean('time', skipna=False))))
    ts_pow_mean = ts_pow_wordevs.isel(time=((ts_pow_wordevs['time']>0) &
                                            (ts_pow_wordevs['time']<1.5))).mean('time')
    ts_pow_mean -= ts_pow_mean.mean('event')
    ts_pow_mean /= ts_pow_mean.std('event')
    ts_pow_mean = ts_pow_mean.isel(
        event=~ts_pow_mean['event'].values['resp_mult'])
    # take the mean over all microwires in an electrode:
    channum = np.array([np.int(c.decode().split('CSC')[-1])
                        for c in ts_pow_mean['channel'].values])
    changroups = xarray.DataArray(data=np.int8(((channum-(channum-1)%8)+7)/8),
                                  name='channel', dims='channel')
    ts_pow_mean = ts_pow_mean.groupby(changroups).mean()
    ts_pow_first = ts_pow_mean.isel(event=ts_pow_mean['event'].values['first'])
    #ts_pow_first = ts_pow_first.isel(
    #    event=~ts_pow_first['event'].values['missing_neuro'])
    subsequent_rec = []
    for ev in ts_pow_first['event'].values:
        tmp = ts_pow_mean.isel(event=ts_pow_mean['event'].values[
            'stim_word'] == ev['stim_word'])
        # assert len(tmp['event']) == 2, 'len event'
        if len(tmp['event']) != 2:
            subsequent_rec.append('---')
            continue
        assert np.alltrue(tmp['event'].values['first']==[True, False])
        subsequent_rec.append(tmp.isel(event=~tmp['event'].values['first'])[
            'event'].values['first_resp'][0])
    subsequent_rec = np.array(subsequent_rec)
    # if 'old' not in subsequent_rec or 'new' not in subsequent_rec:
    #     indx_adjust +=1
    #     continue
    if (np.sum(subsequent_rec=='old')<2) or (np.sum(subsequent_rec=='new')<2):
        print('indxadjust', tspowfile, np.sum(subsequent_rec=='old'),
              np.sum(subsequent_rec=='new'))
        continue
    ts_pow_rec = ts_pow_first.isel(event=subsequent_rec=='old')
    ts_pow_nrec = ts_pow_first.isel(event=subsequent_rec=='new')
    ts_pow_rec.to_hdf(tspowfile.replace('tspow.hdf', 'tspow_rec.hdf'))
    ts_pow_nrec.to_hdf(tspowfile.replace('tspow.hdf', 'tspow_nrec.hdf'))


# tspowfiles_rec = glob(
#     rhino_root+
#     '/scratch/josh/BNI_processed/time_series_pow_wordevs/s*e*cr_tspow_rec.hdf')

# for i, tspowfile in enumerate(tspowfiles_rec):

tspowfiles_all = glob(
    rhino_root+
    '/scratch/josh/BNI_processed/time_series_pow_all/s*e*cr_tspow.hdf')

fig, axs = plt.subplots(41, 1, figsize=(12, 9.5), sharex=True) #, sharey=True)

index_offset = 0
for i, tspowfile_all in enumerate(tspowfiles_all):
    tspowfile_rec = tspowfile_all.replace('all', 'wordevs').replace(
        'tspow', 'tspow_rec')
    if not os.path.exists(tspowfile_rec):
        index_offset +=1
        continue
    print(i+1, '/', len(tspowfiles_all), tspowfile)
    tspowfile_nrec = tspowfile_all.replace('all', 'wordevs').replace(
        'tspow', 'tspow_nrec')
    ts_pow_all = TimeSeries.from_hdf(tspowfile_all)
    ts_pow_all -= ts_pow_all.mean('time')
    ts_pow_all /= ts_pow_all.std('time')
    # take the mean over all microwires in an electrode:
    channum = np.array([np.int(c.decode().split('CSC')[-1])
                        for c in ts_pow_all['channel'].values])
    changroups = xarray.DataArray(data=np.int8(((channum-(channum-1)%8)+7)/8),
                                  name='channel', dims='channel')
    ts_pow_all = ts_pow_all.groupby(changroups).mean()
    ts_pow_rec = TimeSeries.from_hdf(tspowfile_rec)
    ts_pow_nrec = TimeSeries.from_hdf(tspowfile_nrec)
    #target_vec = ts_pow_rec.mean('event').stack(freqchan=('frequency', 'channel'))
    target_vec = (ts_pow_rec.mean('event').stack(freqchan=('frequency', 'channel'))-
                  ts_pow_nrec.mean('event').stack(freqchan=('frequency', 'channel')))
    tv_norm = target_vec/np.sqrt((target_vec**2).sum())
    all_vec = ts_pow_all.stack(freqchan=('frequency', 'channel')).T
    av_norm = all_vec/np.sqrt(all_vec**2).sum(axis=0)
    match = np.dot(tv_norm, av_norm)
    axs[i-index_offset].plot(ts_pow_all['time'][:]-ts_pow_all['time'].values[0],
                             match, alpha=0.5)
    # tsttimes = (ts_pow_rec['event'].values['phase1_time']-
    #             ts_pow_rec['event'].values['phase1_time'][0])/1000000
    tsttimes = (ts_pow_rec['event'].values['phase1_time']/1000000-
                ts_pow_all['time'].values[0])
    xs = np.ravel([(x, x+1.5, x+1.5) for x in tsttimes])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    axs[i-index_offset].fill_between(
        xs, y1s, y2s, where=where, color='red', alpha=0.4)
    tsttimes = (ts_pow_nrec['event'].values['phase1_time']/1000000-
                ts_pow_all['time'].values[0])
    xs = np.ravel([(x, x+1.5, x+1.5) for x in tsttimes])
    y1s = [match.max()]*len(xs)
    y2s = [match.min()]*len(xs)
    where = np.ones(xs.shape, dtype=np.bool)
    where[2::3] = 0
    axs[i-index_offset].fill_between(
        xs, y1s, y2s, where=where, color='yellow', alpha=0.4)

# fig.set_figwidth(12)
# axs[-1].get_xlim()
# Out[42]: (-46.94762426707385, 985.9001096085507)
# In [48]: axs[-1].get_ylim()
# Out[48]: (-0.10393209783416825, 0.10539453012292557)

axs[-1].set_xlim(-2,942)
for ax in axs:
    ax.set_yticks([0])
axs[-1].set_xlabel('Time (s)')
fig.subplots_adjust(left=0.02, bottom=0.05, right=0.999, top=0.999, hspace=0.01)
fig.savefig('./figs/power_sme/lfp_SME_similarity.pdf')






# channum = np.array([np.int(c.decode().split('CSC')[-1]) for c in ts_pow['channel'].values])
# changroups = xarray.DataArray(data=np.int8(((channum-(channum-1)%8)+7)/8),
#                               name='channel', dims='channel')
# ts_pow_chanmean = ts_pow.groupby(changroups).mean()


#fig, axs = plt.subplots(45, 1, figsize=(3, 9.5))

fig, axs = plt.subplots(41, 1, figsize=(1, 9.5), sharex=True, sharey=True)

#fig, axs = plt.subplots(5, 1, figsize=(3, 9.5))
indx_adjust = 0
for i, tspowfile in enumerate(tspowfiles):
    # this file is full of NaNs:
    if ((tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                     's23e15cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's23e5cr_tspow.hdf'))):
        indx_adjust +=1
        continue
    print(tspowfile)
    ts_pow = TimeSeries.from_hdf(tspowfile)
    # remove missing data (selecting based on missing_neuro is not enough)
    ts_pow = ts_pow.isel(event=(~ts_pow['event'].values['missing_neuro']) &
                         np.isfinite(
                             ts_pow.mean('frequency', skipna=False).mean(
                                 'channel', skipna=False).mean(
                                     'time', skipna=False)))
    ts_pow_mean = ts_pow.isel(time=((ts_pow['time']>0) &
                                    (ts_pow['time']<1.5))).mean('time')
    ts_pow_mean -= ts_pow_mean.mean('event')
    ts_pow_mean /= ts_pow_mean.std('event')
    ts_pow_mean = ts_pow_mean.isel(
        event=~ts_pow_mean['event'].values['resp_mult'])
    # take the mean over all microwires in an electrode:
    channum = np.array([np.int(c.decode().split('CSC')[-1])
                        for c in ts_pow_mean['channel'].values])
    changroups = xarray.DataArray(data=np.int8(((channum-(channum-1)%8)+7)/8),
                                  name='channel', dims='channel')
    ts_pow_mean = ts_pow_mean.groupby(changroups).mean()
    ts_pow_first = ts_pow_mean.isel(event=ts_pow_mean['event'].values['first'])
    #ts_pow_first = ts_pow_first.isel(
    #    event=~ts_pow_first['event'].values['missing_neuro'])
    subsequent_rec = []
    for ev in ts_pow_first['event'].values:
        tmp = ts_pow_mean.isel(event=ts_pow_mean['event'].values[
            'stim_word'] == ev['stim_word'])
        # assert len(tmp['event']) == 2, 'len event'
        if len(tmp['event']) != 2:
            subsequent_rec.append('---')
            continue
        assert np.alltrue(tmp['event'].values['first']==[True, False])
        subsequent_rec.append(tmp.isel(event=~tmp['event'].values['first'])[
            'event'].values['first_resp'][0])
    subsequent_rec = np.array(subsequent_rec)
    # if 'old' not in subsequent_rec or 'new' not in subsequent_rec:
    #     indx_adjust +=1
    #     continue
    if (np.sum(subsequent_rec=='old')<2) or (np.sum(subsequent_rec=='new')<2):
        print('indxadjust', tspowfile, np.sum(subsequent_rec=='old'),
              np.sum(subsequent_rec=='new'))
        indx_adjust +=1
        continue
    ts_pow_rec = ts_pow_first.isel(event=subsequent_rec=='old')
    ts_pow_nrec = ts_pow_first.isel(event=subsequent_rec=='new')
    n_rec = len(ts_pow_rec['event'])
    n_nrec = len(ts_pow_nrec['event'])
    sp = np.sqrt(((n_rec-1) * ts_pow_rec.var('event', ddof=1) +
                  (n_nrec-1) * ts_pow_nrec.var('event', ddof=1))/
                 (n_rec+n_nrec+2))
    ts_pow_t = ((ts_pow_rec.mean('event')-ts_pow_nrec.mean('event')) /
                (sp*np.sqrt((1/n_rec)+(1/n_nrec))))
    axs[i-indx_adjust].imshow(ts_pow_t, interpolation='nearest', cmap='bwr',
                              origin='lower', aspect='auto')

fig.subplots_adjust(left=0.475, bottom=0.05, right=0.95, top=0.999, hspace=0.3)
axs[-1].set_xlim((-0.5,8.5))
axs[-1].set_ylim((-1,15))
axs[-1].set_xticks([0, 5])
axs[-1].set_yticks([0, 10])
axs[-2].set_yticklabels(['', ''])
axs[-1].set_yticklabels([np.int(freqs[0]), np.int(np.round(freqs[10]))])
axs[-2].set_ylabel('Frequency')
axs[-1].set_xlabel('Channel')

fig.savefig('./figs/lfp_SME.pdf')






# fig, axs = plt.subplots(41, 1, figsize=(1, 9.5), sharex=True, sharey=True)
fig, axs = plt.subplots(37, 1.25, figsize=(1, 9.5), sharex=True, sharey=True)
ts_pow_t_all = []
#fig, axs = plt.subplots(5, 1, figsize=(3, 9.5))
indx_adjust = 0
prefix = rhino_root+'/scratch/josh/BNI_processed/time_series_pow_wordevs/'
excludedfiles = ['s23e15cr_tspow.hdf','s23e5cr_tspow.hdf','s29e22cr_tspow.hdf',
                 's28e16cr_tspow.hdf','s21e12cr_tspow.hdf','s23e11cr_tspow.hdf',
                 's26e10cr_tspow.hdf','s25e11cr_tspow.hdf','s25e2cr_tspow.hdf',
                 's30e5cr_tspow.hdf']
excludedfiles = [prefix+ef for ef in excludedfiles]
for i, tspowfile in enumerate(tspowfiles):
    # these files are full of NaNs:
    if ((tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                     's23e15cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's23e5cr_tspow.hdf')) or
        # no/incomplete location info for these:
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's29e22cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's28e16cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's21e12cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's23e11cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's26e10cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's25e11cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's25e2cr_tspow.hdf')) or
        (tspowfile == (rhino_root+
                       '/scratch/josh/BNI_processed/time_series_pow_wordevs/'+
                       's30e5cr_tspow.hdf'))):
        indx_adjust +=1
        continue
    print(i+1, '/', len(tspowfiles), tspowfile)
    ts = TimeSeries.from_hdf(tspowfile.replace('_pow_wordevs', '').replace('tspow', 'ts'))
    hemispheres = []
    areas = []
    # chans = ['CSC'+str(i+1) for i in range(len(ts.attrs['cluster'].keys()))]
    for chan in ts.attrs['cluster'].keys():  # chans:
        clusnum = list(ts.attrs['cluster'][chan].keys())[0]
        hemispheres.append(ts.attrs['cluster'][chan][clusnum]['hemisphere'])
        areas.append(ts.attrs['cluster'][chan][clusnum]['area'])
    hemispheres = np.array(hemispheres)
    areas = np.array(areas)
    ts_pow = TimeSeries.from_hdf(tspowfile)
    # remove missing data (selecting based on missing_neuro is not enough)
    ts_pow = ts_pow.isel(event=(~ts_pow['event'].values['missing_neuro']) &
                         np.isfinite(
                             ts_pow.mean('frequency', skipna=False).mean(
                                 'channel', skipna=False).mean(
                                     'time', skipna=False)))
    ts_pow_mean = ts_pow.isel(time=((ts_pow['time']>0) &
                                    (ts_pow['time']<1.5))).mean('time')
    ts_pow_mean -= ts_pow_mean.mean('event')
    ts_pow_mean /= ts_pow_mean.std('event')
    ts_pow_mean = ts_pow_mean.isel(
        event=~ts_pow_mean['event'].values['resp_mult'])
    # # take the mean over all microwires in an electrode:
    # channum = np.array([np.int(c.decode().split('CSC')[-1])
    #                     for c in ts_pow_mean['channel'].values])
    # changroups = xarray.DataArray(data=np.int8(((channum-(channum-1)%8)+7)/8),
    #                               name='channel', dims='channel')
    changroups = np.array(['NA']*len(ts_pow_mean['channel']))
    all_areas = ['H', 'A', 'AC', 'PF']
    for aa in all_areas:
        changroups[areas==aa] = aa
    ts_pow_mean = ts_pow_mean.isel(channel=changroups != 'NA')
    changroups = xarray.DataArray(data=changroups[changroups != 'NA'],
                                  name='channel', dims='channel')
    ts_pow_mean = ts_pow_mean.groupby(changroups).mean()
    ts_pow_1st = ts_pow_mean.isel(event=ts_pow_mean['event'].values['first'])
    ts_pow_2nd = ts_pow_mean.isel(event=~ts_pow_mean['event'].values['first'])
    n_1st = len(ts_pow_first['event'])
    n_2nd = len(ts_pow_second['event'])
    sp = np.sqrt(((n_1st-1) * ts_pow_1st.var('event', ddof=1) +
                  (n_2nd-1) * ts_pow_2nd.var('event', ddof=1))/
                 (n_1st+n_2nd+2))
    ts_pow_t = ((ts_pow_1st.mean('event')-ts_pow_2nd.mean('event')) /
                (sp*np.sqrt((1/n_1st)+(1/n_2nd))))
    ts_pow_t_all.append(ts_pow_t)
    axs[i-indx_adjust].imshow(ts_pow_t, interpolation='nearest', cmap='bwr',
                              origin='lower', aspect='auto')

fig.subplots_adjust(left=0.475, bottom=0.05, right=0.95, top=0.999, hspace=0.3)
axs[-1].set_xlim((-0.5,8.5))
axs[-1].set_ylim((-1,15))
axs[-1].set_xticks([0, 5])
axs[-1].set_yticks([0, 10])
axs[-2].set_yticklabels(['', ''])
axs[-1].set_yticklabels([np.int(freqs[0]), np.int(np.round(freqs[10]))])
axs[-2].set_ylabel('Frequency')
axs[-1].set_xlabel('Channel')

fig.savefig('./figs/lfp/lfp_oldnew.pdf')

expstr_indx = pd.Index([tspf.split('/')[-1].split('_')[0]
                        for tspf in tspowfiles if tspf not in excludedfiles],
                       dtype='U16', name='expstr')
ts_pow_t_all = xarray.concat(ts_pow_t_all, expstr_indx)
fig, axs = plt.subplots(37, 1, figsize=(1.25, 9.5), sharex=True, sharey=True)
for i, expstr in enumerate(ts_pow_t_all['expstr']):
    axs[i].imshow(ts_pow_t_all.sel(expstr=expstr), interpolation='nearest', cmap='bwr',
                  origin='lower', aspect='auto')
axs[-1].set_xlim((-0.5,3.5))
axs[-1].set_ylim((-1,15))
axs[-1].set_xticks(range(4))
axs[-1].set_xticklabels(ts_pow_t_all['channel'].values)
axs[-1].set_yticks([0, 10])
#axs[-2].set_yticklabels(['', ''])
axs[-1].set_yticklabels([np.int(freqs[0]), np.int(np.round(freqs[10]))])
axs[-2].set_ylabel('Frequency')
axs[-1].set_xlabel('Area')
#fig.subplots_adjust(left=0.475, bottom=0.05, right=0.95, top=0.999, hspace=0.3)
fig.subplots_adjust(left=0.38, bottom=0.045, right=0.97, top=0.999, hspace=0.3)
#fig.set_size_inches((1.25, 9.5))
fig.savefig('./figs/lfp/lfp_oldnew_area.pdf')

fig, axs = plt.subplots(1, 1, figsize=(1.25, 0.5), sharex=True, sharey=True)
axs.imshow(ts_pow_t_all.mean('expstr'), interpolation='nearest', cmap='bwr',
           origin='lower', aspect='auto')
axs.set_xlim((-0.5,3.5))
axs.set_ylim((-1,15))
axs.set_xticks(range(4))
axs.set_xticklabels(ts_pow_t_all['channel'].values)
axs.set_yticks([0, 10])
#axs[-2].set_yticklabels(['', ''])
axs.set_yticklabels([np.int(freqs[0]), np.int(np.round(freqs[10]))])
axs.set_ylabel('Freq.')
axs.set_xlabel('Area')
#fig.subplots_adjust(left=0.475, bottom=0.05, right=0.95, top=0.999, hspace=0.3)
fig.subplots_adjust(left=0.38, bottom=0.45, right=0.97, top=0.999, hspace=0.3)
fig.savefig('./figs/lfp/lfp_oldnew_area_mean.pdf')




ts_files = sorted(glob('/home/ctw/data/BNI/data/time_series/s*_ts.hdf'))
# ts_file = ts_files[0]

for ts_file in ts_files:
    ts = TimeSeries.from_hdf(ts_file)
    hemispheres = []
    areas = []
    chans = ['CSC'+str(i+1) for i in range(len(ts.attrs['cluster'].keys()))]
    for chan in chans:
        hemispheres.append(ts.attrs['cluster'][chan]['1']['hemisphere'])
        areas.append(ts.attrs['cluster'][chan]['1']['area'])
    
