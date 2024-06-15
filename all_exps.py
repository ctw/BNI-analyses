
import numpy as np
from glob import glob
#import os
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
import pickle
from scipy.stats import ks_2samp
from scipy import stats


exp_2cat = ['oi', 'ir']
exp_1cat = ['sc', 'rp', 'ReberEtAl2019']
exp_0cat = ['crm', 'sr']



# exps = ['oi', 'ir', 'sc', 'rp', 'ReberEtAl2019']


all_ps = {}
comparisons = {}

for exp in exp_0cat:
    all_ps[exp] = pickle.load(open('./'+exp+'/data/all_ps_zmean.pickle', 'rb'))
    comparisons[exp] = {}
    for area in all_ps[exp]:
        comparisons[exp][area] = {'desc': {}}
        for comp in comparisons[exp][area]:
            comparisons[exp][area][comp]['ks'] = []
            comparisons[exp][area][comp]['diff'] = []
            comparisons[exp][area][comp]['d_prime'] = []
            comparisons[exp][area][comp]['t'] = []
            comparisons[exp][area][comp]['t_diff'] = []
            comparisons[exp][area][comp]['n'] = []
        for comp in ['desc']:
            comparisons[exp][area][comp]['mean'] = []
            comparisons[exp][area][comp]['mean_perm'] = []
            comparisons[exp][area][comp]['std_perm'] = []
        sm = np.array(all_ps[exp][area]['sim_mean'])
        spm = np.array(all_ps[exp][area]['pmean'])      
        comparisons[exp][area]['desc']['ks'].append(
            ks_2samp(sm, spm))
        comparisons[exp][area]['desc']['mean'].append(np.nanmean(sm))
        comparisons[exp][area]['desc']['mean_perm'].append(np.nanmean(spm))
        comparisons[exp][area]['desc']['std_perm'].append(np.nanstd(spm))
        comparisons[exp][area]['desc']['diff'].append(
            np.nanmean(sm)-np.nanmean(spm))
        comparisons[exp][area]['desc']['d_prime'].append(
            (np.nanmean(sm)-np.nanmean(spm))/np.nanstd(spm))
        comparisons[exp][area]['desc']['t'].append(
            np.nanmean(sm-spm)/(np.nanstd(sm-spm)/np.sqrt(len(sm))))
        comparisons[exp][area]['desc']['t_diff'].append(
            (np.array(sm)-np.array(spm))/(np.nanstd(sm-spm)/np.sqrt(len(sm))))
        comparisons[exp][area]['desc']['n'].append(len(sm))

for exp in exp_1cat:
    all_ps[exp] = pickle.load(open('./'+exp+'/data/all_ps_zmean.pickle', 'rb'))
    comparisons[exp] = {}
    for area in all_ps[exp]:
        comparisons[exp][area] = {'desc': {}, 'cat': {}, 'desc-cat': {}, 'permcomp': {}}
        for comp in comparisons[exp][area]:
            comparisons[exp][area][comp]['ks'] = []
            comparisons[exp][area][comp]['diff'] = []
            comparisons[exp][area][comp]['d_prime'] = []
            comparisons[exp][area][comp]['t'] = []
            comparisons[exp][area][comp]['t_diff'] = []
            comparisons[exp][area][comp]['n'] = []
        for comp in ['desc', 'cat']:
            comparisons[exp][area][comp]['mean'] = []
            comparisons[exp][area][comp]['mean_perm'] = []
            comparisons[exp][area][comp]['std_perm'] = []
            comparisons[exp][area]['desc-cat']['mean_desc'] = []
        comparisons[exp][area]['desc-cat']['mean_cat'] = []
        comparisons[exp][area]['permcomp']['mean_desc_perm'] = []
        comparisons[exp][area]['permcomp']['mean_cat_perm'] = []
        for sm, spm in zip(all_ps[exp][area]['sim_means'],
                           all_ps[exp][area]['sim_perm_means']):
            # for sm, spm in zip(all_ps[exp][area]['sim_mean'],
            #                    all_ps[exp][area]['pmean']):
            comparisons[exp][area]['desc']['ks'].append(
                ks_2samp(sm['description'], spm['description']))
            comparisons[exp][area]['desc']['mean'].append(sm['description'].mean())
            comparisons[exp][area]['desc']['mean_perm'].append(spm['description'].mean())
            comparisons[exp][area]['desc']['std_perm'].append(spm['description'].std())
            comparisons[exp][area]['desc']['diff'].append(
                sm['description'].mean()-spm['description'].mean())
            comparisons[exp][area]['desc']['d_prime'].append(
                (sm['description'].mean()-spm['description'].mean())/spm['description'].std())
            comparisons[exp][area]['desc']['t'].append(
                np.mean(sm['description']-spm['description'])/
                (np.std(sm['description']-spm['description'])/
                 np.sqrt(len(sm['description']))))
            comparisons[exp][area]['desc']['t_diff'].append(
                (sm['description']-spm['description'])/
                (np.std(sm['description']-spm['description'])/
                 np.sqrt(len(sm['description']))))
            comparisons[exp][area]['desc']['n'].append(len(sm['description']))
            comparisons[exp][area]['cat']['ks'].append(ks_2samp(sm['category'], spm['category']))
            comparisons[exp][area]['cat']['mean'].append(sm['category'].mean())
            comparisons[exp][area]['cat']['mean_perm'].append(spm['category'].mean())
            comparisons[exp][area]['cat']['std_perm'].append(spm['category'].std())
            comparisons[exp][area]['cat']['diff'].append(
                sm['category'].mean()-spm['category'].mean())
            comparisons[exp][area]['cat']['d_prime'].append(
                (sm['category'].mean()-spm['category'].mean())/spm['category'].std())
            comparisons[exp][area]['cat']['t'].append(
                np.mean(sm['category']-spm['category'])/
                (np.std(sm['category']-spm['category'])/
                 np.sqrt(len(sm['category']))))
            comparisons[exp][area]['cat']['t_diff'].append(
                (sm['category']-spm['category'])/
                (np.std(sm['category']-spm['category'])/
                 np.sqrt(len(sm['category']))))
            comparisons[exp][area]['cat']['n'].append(len(sm['category']))
            comparisons[exp][area]['desc-cat']['ks'].append(
            ks_2samp(sm['description'], sm['category']))
            comparisons[exp][area]['desc-cat']['mean_desc'].append(sm['description'].mean())
            comparisons[exp][area]['desc-cat']['mean_cat'].append(sm['category'].mean())
            comparisons[exp][area]['desc-cat']['diff'].append(
                sm['description'].mean()-sm['category'].mean())
            comparisons[exp][area]['permcomp']['ks'].append(
                ks_2samp(spm['description'], spm['category']))
            comparisons[exp][area]['permcomp']['mean_desc_perm'].append(spm['description'].mean())
            comparisons[exp][area]['permcomp']['mean_cat_perm'].append(spm['category'].mean())
            comparisons[exp][area]['permcomp']['diff'].append(
                spm['description'].mean()-sm['category'].mean())


for exp in exp_2cat:
    all_ps[exp] = pickle.load(open('./'+exp+'/data/all_ps_zmean_cat2.pickle', 'rb'))
    comparisons[exp] = {}
    for area in all_ps[exp]:
        comparisons[exp][area] = {'desc': {}, 'cat': {}, 'cat2': {}, 'cat2img': {},
                             'cat2txt': {}, 'cat2imgtxt': {}, 'desc-cat': {},
                             'permcomp': {}}
        for comp in comparisons[exp][area]:
            comparisons[exp][area][comp]['ks'] = []
            comparisons[exp][area][comp]['diff'] = []
            comparisons[exp][area][comp]['d_prime'] = []
            comparisons[exp][area][comp]['t'] = []
            comparisons[exp][area][comp]['t_diff'] = []
        for comp in ['desc', 'cat', 'cat2', 'cat2img', 'cat2txt', 'cat2imgtxt']:
            comparisons[exp][area][comp]['mean'] = []
            comparisons[exp][area][comp]['mean_perm'] = []
            comparisons[exp][area][comp]['std_perm'] = []
        comparisons[exp][area]['desc-cat']['mean_desc'] = []
        comparisons[exp][area]['desc-cat']['mean_cat'] = []
        comparisons[exp][area]['permcomp']['mean_desc_perm'] = []
        comparisons[exp][area]['permcomp']['mean_cat_perm'] = []

        for sm, spm in zip(all_ps[exp][area]['sim_means'],
                           all_ps[exp][area]['sim_perm_means']):
            comparisons[exp][area]['desc']['ks'].append(
                ks_2samp(sm['description'], spm['description']))
            comparisons[exp][area]['desc']['mean'].append(sm['description'].mean())
            comparisons[exp][area]['desc']['mean_perm'].append(spm['description'].mean())
            comparisons[exp][area]['desc']['std_perm'].append(spm['description'].std())
            comparisons[exp][area]['desc']['diff'].append(
                sm['description'].mean()-spm['description'].mean())
            comparisons[exp][area]['desc']['d_prime'].append(
                (sm['description'].mean()-spm['description'].mean())/spm['description'].std())
            comparisons[exp][area]['desc']['t'].append(
                np.mean(sm['description']-spm['description'])/(np.std(
                    sm['description']-spm['description'])/np.sqrt(len(sm['description']))))
            comparisons[exp][area]['desc']['t_diff'].append(
                (sm['description']-spm['description'])/(np.std(
                    sm['description']-spm['description'])/np.sqrt(len(sm['description']))))
            #
            comparisons[exp][area]['cat']['ks'].append(ks_2samp(sm['category'], spm['category']))
            comparisons[exp][area]['cat']['mean'].append(sm['category'].mean())
            comparisons[exp][area]['cat']['mean_perm'].append(spm['category'].mean())
            comparisons[exp][area]['cat']['std_perm'].append(spm['category'].std())
            comparisons[exp][area]['cat']['diff'].append(
                sm['category'].mean()-spm['category'].mean())
            comparisons[exp][area]['cat']['d_prime'].append(
                (sm['category'].mean()-spm['category'].mean())/spm['category'].std())
            comparisons[exp][area]['cat']['t'].append(
                np.mean(sm['category']-spm['category'])/(np.std(
                    sm['category']-spm['category'])/np.sqrt(len(sm['category']))))
            comparisons[exp][area]['cat']['t_diff'].append(
                (sm['category']-spm['category'])/(np.std(
                    sm['category']-spm['category'])/np.sqrt(len(sm['category']))))
            #
            comparisons[exp][area]['cat2']['ks'].append(ks_2samp(sm['category2'], spm['category2']))
            comparisons[exp][area]['cat2']['mean'].append(sm['category2'].mean())
            comparisons[exp][area]['cat2']['mean_perm'].append(spm['category2'].mean())
            comparisons[exp][area]['cat2']['std_perm'].append(spm['category2'].std())
            comparisons[exp][area]['cat2']['diff'].append(
                sm['category2'].mean()-spm['category2'].mean())
            comparisons[exp][area]['cat2']['d_prime'].append(
                (sm['category2'].mean()-spm['category2'].mean())/spm['category2'].std())
            comparisons[exp][area]['cat2']['t'].append(
                np.mean(sm['category2']-spm['category2'])/(np.std(
                    sm['category2']-spm['category2'])/np.sqrt(len(sm['category2']))))
            comparisons[exp][area]['cat2']['t_diff'].append(
                (sm['category2']-spm['category2'])/(np.std(
                    sm['category2']-spm['category2'])/np.sqrt(len(sm['category2']))))
            #
            comparisons[exp][area]['cat2img']['ks'].append(ks_2samp(
                sm['category2_img'], spm['category2_img']))
            comparisons[exp][area]['cat2img']['mean'].append(
                sm['category2_img'].mean())
            comparisons[exp][area]['cat2img']['mean_perm'].append(
                spm['category2_img'].mean())
            comparisons[exp][area]['cat2img']['std_perm'].append(
                spm['category2_img'].std())
            comparisons[exp][area]['cat2img']['diff'].append(
                sm['category2_img'].mean()-spm['category2_img'].mean())
            comparisons[exp][area]['cat2img']['d_prime'].append(
                (sm['category2_img'].mean()-spm['category2_img'].mean()) /
                spm['category2_img'].std())
            comparisons[exp][area]['cat2img']['t'].append(
                np.mean(sm['category2_img']-spm['category2_img'])/(np.std(
                    sm['category2_img']-spm['category2_img'])/np.sqrt(len(sm['category2_img']))))
            comparisons[exp][area]['cat2img']['t_diff'].append(
                (sm['category2_img']-spm['category2_img'])/(np.std(
                    sm['category2_img']-spm['category2_img'])/np.sqrt(len(sm['category2_img']))))
            #
            comparisons[exp][area]['cat2txt']['ks'].append(ks_2samp(
                sm['category2_txt'], spm['category2_txt']))
            comparisons[exp][area]['cat2txt']['mean'].append(
                sm['category2_txt'].mean())
            comparisons[exp][area]['cat2txt']['mean_perm'].append(
                spm['category2_txt'].mean())
            comparisons[exp][area]['cat2txt']['std_perm'].append(
                spm['category2_txt'].std())
            comparisons[exp][area]['cat2txt']['diff'].append(
                sm['category2_txt'].mean()-spm['category2_txt'].mean())
            comparisons[exp][area]['cat2txt']['d_prime'].append(
                (sm['category2_txt'].mean()-spm['category2_txt'].mean()) /
                spm['category2_txt'].std())
            comparisons[exp][area]['cat2txt']['t'].append(
                np.mean(sm['category2_txt']-spm['category2_txt'])/(np.std(
                    sm['category2_txt']-spm['category2_txt'])/np.sqrt(len(sm['category2_txt']))))
            comparisons[exp][area]['cat2txt']['t_diff'].append(
                (sm['category2_txt']-spm['category2_txt'])/(np.std(
                    sm['category2_txt']-spm['category2_txt'])/np.sqrt(len(sm['category2_txt']))))
            #
            comparisons[exp][area]['cat2imgtxt']['ks'].append(ks_2samp(
                sm['category2_imgtxt'], spm['category2_imgtxt']))
            comparisons[exp][area]['cat2imgtxt']['mean'].append(
                sm['category2_imgtxt'].mean())
            comparisons[exp][area]['cat2imgtxt']['mean_perm'].append(
                spm['category2_imgtxt'].mean())
            comparisons[exp][area]['cat2imgtxt']['std_perm'].append(
                spm['category2_imgtxt'].std())
            comparisons[exp][area]['cat2imgtxt']['diff'].append(
                sm['category2_imgtxt'].mean()-spm['category2_imgtxt'].mean())
            comparisons[exp][area]['cat2imgtxt']['d_prime'].append(
                (sm['category2_imgtxt'].mean()-spm['category2_imgtxt'].mean()) /
                spm['category2_imgtxt'].std())
            comparisons[exp][area]['cat2imgtxt']['t'].append(
                np.mean(sm['category2_imgtxt']-spm['category2_imgtxt'])/(np.std(
                    sm['category2_imgtxt']-spm['category2_imgtxt'])/np.sqrt(len(sm['category2_imgtxt']))))
            comparisons[exp][area]['cat2imgtxt']['t_diff'].append(
                (sm['category2_imgtxt']-spm['category2_imgtxt'])/(np.std(
                    sm['category2_imgtxt']-spm['category2_imgtxt'])/np.sqrt(len(sm['category2_imgtxt']))))
            #
            comparisons[exp][area]['desc-cat']['ks'].append(
            ks_2samp(sm['description'], sm['category']))
            comparisons[exp][area]['desc-cat']['mean_desc'].append(sm['description'].mean())
            comparisons[exp][area]['desc-cat']['mean_cat'].append(sm['category'].mean())
            comparisons[exp][area]['desc-cat']['diff'].append(
                sm['description'].mean()-sm['category'].mean())
            comparisons[exp][area]['permcomp']['ks'].append(
                ks_2samp(spm['description'], spm['category']))
            comparisons[exp][area]['permcomp']['mean_desc_perm'].append(spm['description'].mean())
            comparisons[exp][area]['permcomp']['mean_cat_perm'].append(spm['category'].mean())
            comparisons[exp][area]['permcomp']['diff'].append(
                spm['description'].mean()-sm['category'].mean())



# exp_2cat = ['oi', 'ir']
# exp_1cat = ['sc', 'rp', 'ReberEtAl2019']
# exp_0cat = ['crm', 'sr']



area = 'all'
fig, axs = plt.subplots(1, 1, figsize=(5,5))
for i, exp in enumerate(exp_0cat):
    axs.plot(i*2-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t_diff'][0]))-0.5)*0.3,
             comparisons[exp][area]['desc']['t_diff'][0], 'ko', alpha=0.5)
    axs.boxplot(np.array(comparisons[exp][area]['desc']['t_diff'])[
        np.isfinite(comparisons[exp][area]['desc']['t_diff'])], positions=[i*2-0.35],
                notch=True, bootstrap=10000, widths=0.65, showfliers=False,
                # notch=False, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))
axs.set_yscale('symlog')
axs.set_xticklabels(exp_0cat)
axs.set_ylabel('Normalized Similarity ($t$)')
axs.set_xlabel('Recognition Experiment')
axs.set_xticklabels(['Visual (CRM)', 'Auditory (SR)'])
#axs.set_xticklabels(['Continuous Recognition Memory', 'Spoken Recognition'])
axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
fig.subplots_adjust(left=0.15, right=0.99, bottom=0.09, top=0.99)
fig.savefig('./figs/exp_0cat.pdf')
# fig.savefig('./figs/exp_0cat_nonotch.pdf')



for i, exp in enumerate(exp_0cat):
    t = np.nanmean(comparisons[exp][area]['desc']['t_diff'][0])/(
        np.nanstd(comparisons[exp][area]['desc']['t_diff'][0])/
        np.sqrt(np.sum(np.isfinite(comparisons[exp][area]['desc']['t_diff'][0]))))
    df = np.sum(np.isfinite(comparisons[exp][area]['desc']['t_diff'][0]))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# crm 2.799956824088979 45 0.00750355824490666
# sr 1.8454971496354349 15 0.08479711409074246

exp_txt = ['Visual (CRM)', 'Auditory (SR)']
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(6, 3))
for i, exp in enumerate(exp_0cat):
    filt = np.isfinite(comparisons[exp][area]['desc']['t_diff'][0])
    axs[i].plot((np.array(all_ps[exp][area]['hit'])[filt]-
                 np.array(all_ps[exp][area]['fa'])[filt]),
                comparisons[exp][area]['desc']['t_diff'][0][filt],
                'ko', alpha=0.5)
    axs[i].set_yscale('symlog')
    axs[i].set_box_aspect(1)
    axs[i].set_xlabel('Recognition Performance')
    axs[i].text(0.05, 0.1, exp_txt[i],
                horizontalalignment='left', size='large',
                weight=1000, verticalalignment='center',
                transform=axs[i].transAxes)
    axs[i].text(0.2, 0.2, r'$r=$'+str(
        np.round(np.corrcoef(
            (np.array(all_ps[exp][area]['hit'])[filt]-
             np.array(all_ps[exp][area]['fa'])[filt]),
            comparisons[exp][area]['desc']['t_diff'][0][filt])[0,1], 2)),
                horizontalalignment='center', verticalalignment='center',
                transform=axs[i].transAxes)
axs[0].set_ylabel('Normalized Similarity ($t$)')
#axs[0].set_xlabel('Experiment')
#axs[0].set_xticklabels(['Visual (CRM)', 'Auditory (SR)'])
fig.subplots_adjust(left=0.12, right=0.99, bottom=0.1, top=0.999)
fig.savefig('./figs/exp_0cat_recogperf.pdf')

for i, exp in enumerate(exp_0cat):
    filt = np.isfinite(comparisons[exp][area]['desc']['t_diff'][0])
    r = np.corrcoef(
        (np.array(all_ps[exp][area]['hit'])[filt]-
         np.array(all_ps[exp][area]['fa'])[filt]),
        comparisons[exp][area]['desc']['t_diff'][0][filt])[0, 1]
    df = np.sum(np.isfinite(comparisons[exp][area]['desc']['t_diff'][0]))-2
    t = np.abs((r*np.sqrt(df))/np.sqrt(1-r**2))
    pval = stats.t.sf(t, df)*2
    print(exp, r, t, df, pval)
# crm -0.3629903026295787 2.584056209008891 44 0.013161106406128442
# sr -0.1734882907895884 0.6591288189832707 14 0.5205078988514106

# area = 'all'
# fig, axs = plt.subplots(1, 1, figsize=(8,5))
# for i, exp in enumerate(exp_1cat):
#     axs.plot(i*2-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t']))-0.5)*0.3,
#              comparisons[exp][area]['desc']['t'], 'ko', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['desc']['t'])[
#         np.isfinite(comparisons[exp][area]['desc']['t'])], positions=[i*2-0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     axs.plot(i*2+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat']['t'], 'bo', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat']['t'])[
#         np.isfinite(comparisons[exp][area]['cat']['t'])], positions=[i*2+0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
# axs.set_yscale('symlog')
# axs.set_xticks(range(0, len(exp_1cat)+2, 2))
# exp1_cat_labels = exp_1cat.copy()
# exp1_cat_labels[-1] = 'Reber et al. (2019)'
# axs.set_xticklabels(exp1_cat_labels)
# axs.set_ylabel('Normalized Similarity')
# axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
# fig.subplots_adjust(left=0.09, right=0.99, bottom=0.05, top=0.99)
# fig.savefig('./figs/exp_1cat.pdf')


area = 'all'
fig, axs = plt.subplots(1, 1, figsize=(7,5))
for i, exp in enumerate(exp_1cat):
    axs.plot(i*2-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t']))-0.5)*0.3,
             comparisons[exp][area]['desc']['t'], 'ko', alpha=0.5)
    axs.boxplot(np.array(comparisons[exp][area]['desc']['t'])[
        np.isfinite(comparisons[exp][area]['desc']['t'])], positions=[i*2-0.35],
                notch=True, bootstrap=10000, widths=0.65, showfliers=False,
                # notch=False, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))
    axs.plot(i*2+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat']['t']))-0.5)*0.3,
             comparisons[exp][area]['cat']['t'], 'bo', alpha=0.5)
    axs.boxplot(np.array(comparisons[exp][area]['cat']['t'])[
        np.isfinite(comparisons[exp][area]['cat']['t'])], positions=[i*2+0.35],
                notch=True, bootstrap=10000, widths=0.65, showfliers=False,
                # notch=False, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))

axs.set_yscale('symlog')
assert axs.get_ylim()[0] > -200, 'min lim'
assert axs.get_ylim()[1] < 4000, 'max lim'
axs.set_ylim((-200, 4000))

for i, exp in enumerate(exp_1cat):
    # axs.text(i*2-0.35, 0.7*axs.get_ylim()[1], 'repetition', verticalalignment='center',
    #          horizontalalignment='center')
    # axs.text(i*2.+0.35, 0.7*axs.get_ylim()[1], 'broad', verticalalignment='center',
    #          horizontalalignment='center')
    axs.text(i*2-0.35, 0.7*axs.get_ylim()[1], 'repetition', verticalalignment='center',
             horizontalalignment='center')
    axs.text(i*2.+0.35, 0.7*axs.get_ylim()[1], 'category', verticalalignment='center',
             horizontalalignment='center')

#axs.set_yscale('symlog')
axs.set_xticks(range(0, len(exp_1cat)+2, 2))
exp1_cat_labels = exp_1cat.copy()
exp1_cat_labels[-1] = 'Reber et al. (2019)'
axs.set_xticklabels(exp1_cat_labels)
axs.set_ylabel('Normalized Similarity ($t$)')
axs.set_xlabel('Discrimination Experiment')
axs.set_xticklabels(['Split Category (SC)', 'Rossion Pourtois (RP)', 'Reber et al. (2019)'])
axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
fig.subplots_adjust(left=0.12, right=0.99, bottom=0.09, top=0.99)
fig.savefig('./figs/exp_1cat.pdf')
# fig.savefig('./figs/exp_1cat_nonotch.pdf')


for i, exp in enumerate(exp_1cat):
    t = np.nanmean(comparisons[exp][area]['desc']['t'])/(
        np.nanstd(comparisons[exp][area]['desc']['t'])/
        np.sqrt(np.sum(np.isfinite(comparisons[exp][area]['desc']['t']))))
    df = np.sum(np.isfinite(comparisons[exp][area]['desc']['t']))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# sc 4.373112016057351 18 0.0003666387712854639
# rp 5.1232156438664616 34 1.1889125559601344e-05
# ReberEtAl2019 5.833793185761428 58 2.557507255815159e-07
for i, exp in enumerate(exp_1cat):
    t = np.nanmean(comparisons[exp][area]['cat']['t'])/(
        np.nanstd(comparisons[exp][area]['cat']['t'])/
        np.sqrt(np.sum(np.isfinite(comparisons[exp][area]['cat']['t']))))
    df = np.sum(np.isfinite(comparisons[exp][area]['cat']['t']))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# sc 2.7221350095341754 18 0.013978158434426433
# rp 3.4396131322015537 34 0.0015587150114037972
# ReberEtAl2019 8.286580505286153 58 2.0301186317080976e-11
for i, exp in enumerate(exp_1cat):
    desc = np.array(comparisons[exp][area]['desc']['t'])
    cat = np.array(comparisons[exp][area]['cat']['t'])
    t = np.nanmean(desc-cat)/(np.nanstd(desc-cat)/
                              np.sqrt(np.sum(np.isfinite(desc-cat))))
    df = np.sum(np.isfinite(desc-cat))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# sc 4.3185821872421855 18 0.00041370141808557114
# rp 4.553550201928133 34 6.463132550918978e-05
# ReberEtAl2019 3.8713761818751813 58 0.00027714920703981457




# area = 'all'
# fig, axs = plt.subplots(1, 1, figsize=(12,5))
# for i, exp in enumerate(exp_2cat):
#     axs.plot(i*5.5-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t']))-0.5)*0.3,
#              comparisons[exp][area]['desc']['t'], 'ko', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['desc']['t'])[
#         np.isfinite(comparisons[exp][area]['desc']['t'])], positions=[i*5.5-0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     axs.plot(i*5.5+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat']['t'], 'bo', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat']['t'])[
#         np.isfinite(comparisons[exp][area]['cat']['t'])], positions=[i*5.5+0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*5.5+0.35*3+(np.random.uniform(size=len(comparisons[exp][area]['cat2']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2']['t'])], positions=[i*5.5+0.35*3],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
    
#     axs.plot(i*5.5+0.35*5+(np.random.uniform(size=len(comparisons[exp][area]['cat2img']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2img']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2img']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2img']['t'])], positions=[i*5.5+0.35*5],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*5.5+0.35*7+(np.random.uniform(size=len(comparisons[exp][area]['cat2txt']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2txt']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2txt']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2txt']['t'])], positions=[i*5.5+0.35*7],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*5.5+0.35*9+(np.random.uniform(size=len(comparisons[exp][area]['cat2imgtxt']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2imgtxt']['t'], 'go', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2imgtxt']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2imgtxt']['t'])], positions=[i*5.5+0.35*9],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))

# for i, exp in enumerate(exp_2cat):
#     axs.text(i*5.5+0.35, 0.6*axs.get_ylim()[1], 'broad', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*3, 0.6*axs.get_ylim()[1], 'narrow', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*5, 0.6*axs.get_ylim()[1], 'img only', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*7, 0.6*axs.get_ylim()[1], 'txt only', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*9, 0.6*axs.get_ylim()[1], 'img - txt', verticalalignment='center',
#              horizontalalignment='center')
# axs.set_yscale('symlog')
# axs.set_xticks(np.arange(1.5, len(exp_2cat)+6, 5.5))
# axs.set_xticklabels(exp_2cat)
# axs.set_ylabel('Normalized Similarity')
# axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
# fig.subplots_adjust(left=0.06, right=0.99, bottom=0.05, top=0.99)
# fig.savefig('./figs/exp_2cat_old.pdf')


# area = 'all'
# fig, axs = plt.subplots(1, 1, figsize=(12,5))
# for i, exp in enumerate(exp_2cat):
#     axs.plot(i*5.5-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t']))-0.5)*0.3,
#              comparisons[exp][area]['desc']['t'], 'ko', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['desc']['t'])[
#         np.isfinite(comparisons[exp][area]['desc']['t'])], positions=[i*5.5-0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     axs.plot(i*5.5+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat']['t'], 'bo', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat']['t'])[
#         np.isfinite(comparisons[exp][area]['cat']['t'])], positions=[i*5.5+0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*5.5+0.35*3+(np.random.uniform(size=len(comparisons[exp][area]['cat2']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2']['t'])], positions=[i*5.5+0.35*3],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
    
#     axs.plot(i*5.5+0.35*5+(np.random.uniform(size=len(comparisons[exp][area]['cat2img']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2img']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2img']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2img']['t'])], positions=[i*5.5+0.35*5],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*5.5+0.35*7+(np.random.uniform(size=len(comparisons[exp][area]['cat2txt']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2txt']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2txt']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2txt']['t'])], positions=[i*5.5+0.35*7],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*5.5+0.35*9+(np.random.uniform(size=len(comparisons[exp][area]['cat2imgtxt']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2imgtxt']['t'], 'go', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2imgtxt']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2imgtxt']['t'])], positions=[i*5.5+0.35*9],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))

# for i, exp in enumerate(exp_2cat):
#     axs.text(i*5.5+0.35, 0.6*axs.get_ylim()[1], 'broad', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*3, 0.6*axs.get_ylim()[1], 'narrow', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*5, 0.6*axs.get_ylim()[1], 'img only', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*7, 0.6*axs.get_ylim()[1], 'txt only', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*5.5+0.35*9, 0.6*axs.get_ylim()[1], 'img - txt', verticalalignment='center',
#              horizontalalignment='center')
# axs.set_yscale('symlog')
# axs.set_xticks(np.arange(1.5, len(exp_2cat)+6, 5.5))
# axs.set_xticklabels(exp_2cat)
# axs.set_ylabel('Normalized Similarity')
# axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
# fig.subplots_adjust(left=0.06, right=0.99, bottom=0.05, top=0.99)
# fig.savefig('./figs/exp_2cat_detailed.pdf')



# area = 'all'
# fig, axs = plt.subplots(1, 1, figsize=(6,5))
# for i, exp in enumerate(exp_2cat):
#     axs.plot(i*2.5-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t']))-0.5)*0.3,
#              comparisons[exp][area]['desc']['t'], 'ko', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['desc']['t'])[
#         np.isfinite(comparisons[exp][area]['desc']['t'])], positions=[i*2.5-0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     axs.plot(i*2.5+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat']['t'], 'bo', alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat']['t'])[
#         np.isfinite(comparisons[exp][area]['cat']['t'])], positions=[i*2.5+0.35],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))
#     #
#     axs.plot(i*2.5+0.35*3+(np.random.uniform(size=len(comparisons[exp][area]['cat2']['t']))-0.5)*0.3,
#              comparisons[exp][area]['cat2']['t'], 'go',  alpha=0.5)
#     axs.boxplot(np.array(comparisons[exp][area]['cat2']['t'])[
#         np.isfinite(comparisons[exp][area]['cat2']['t'])], positions=[i*2.5+0.35*3],
#                 notch=True, bootstrap=10000, widths=0.65, showfliers=False,
#                 whiskerprops=dict(lw=0), capprops=dict(lw=0))

# axs.set_yscale('symlog')
# assert axs.get_ylim()[0] > -80, 'min lim'
# assert axs.get_ylim()[1] < 1700, 'max lim'
# axs.set_ylim((-100, 2000))

# for i, exp in enumerate(exp_2cat):
#     axs.text(i*2.5-0.35, 0.7*axs.get_ylim()[1], 'repetition', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*2.5+0.35, 0.7*axs.get_ylim()[1], 'broad', verticalalignment='center',
#              horizontalalignment='center')
#     axs.text(i*2.5+0.35*3, 0.7*axs.get_ylim()[1], 'narrow', verticalalignment='center',
#              horizontalalignment='center')
# #axs.set_xticks(np.arange(1.5, len(exp_2cat)+6, 3.5))
# axs.set_xticks(np.arange(0.35, len(exp_2cat)+3, 2.5))
# axs.set_xticklabels(exp_2cat)
# axs.set_ylabel('Normalized Similarity')
# axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
# fig.subplots_adjust(left=0.12, right=0.99, bottom=0.05, top=0.99)
# fig.savefig('./figs/exp_2cat_oldorder.pdf')







area = 'all'
fig, axs = plt.subplots(1, 1, figsize=(7,5))
notch = True
# notch = False
for i, exp in enumerate(exp_2cat):
    axs.plot(i*2.5-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t']))-0.5)*0.3,
             comparisons[exp][area]['desc']['t'], 'ko', alpha=0.5)
    axs.boxplot(np.array(comparisons[exp][area]['desc']['t'])[
        np.isfinite(comparisons[exp][area]['desc']['t'])], positions=[i*2.5-0.35],
                notch=notch, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))
    #
    axs.plot(i*2.5+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat2']['t']))-0.5)*0.3,
             comparisons[exp][area]['cat2']['t'], 'go',  alpha=0.5)
    axs.boxplot(np.array(comparisons[exp][area]['cat2']['t'])[
        np.isfinite(comparisons[exp][area]['cat2']['t'])], positions=[i*2.5+0.35],
                notch=notch, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))
    #
    axs.plot(i*2.5+0.35*3+(np.random.uniform(size=len(comparisons[exp][area]['cat']['t']))-0.5)*0.3,
             comparisons[exp][area]['cat']['t'], 'bo', alpha=0.5)
    axs.boxplot(np.array(comparisons[exp][area]['cat']['t'])[
        np.isfinite(comparisons[exp][area]['cat']['t'])], positions=[i*2.5+0.35*3],
                notch=notch, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))

axs.set_yscale('symlog')
assert axs.get_ylim()[0] > -80, 'min lim'
assert axs.get_ylim()[1] < 1700, 'max lim'
axs.set_ylim((-100, 2000))

for i, exp in enumerate(exp_2cat):
    # axs.text(i*2.5-0.35, 0.7*axs.get_ylim()[1], 'repetition', verticalalignment='center',
    #          horizontalalignment='center')
    # axs.text(i*2.5+0.35, 0.7*axs.get_ylim()[1], 'narrow', verticalalignment='center',
    #          horizontalalignment='center')
    # axs.text(i*2.5+0.35*3, 0.7*axs.get_ylim()[1], 'broad', verticalalignment='center',
    #          horizontalalignment='center')
    axs.text(i*2.5-0.35, 0.7*axs.get_ylim()[1], 'repetition', verticalalignment='center',
             horizontalalignment='center')
    axs.text(i*2.5+0.35, 0.7*axs.get_ylim()[1], 'object', verticalalignment='center',
             horizontalalignment='center')
    axs.text(i*2.5+0.35*3, 0.7*axs.get_ylim()[1], 'category', verticalalignment='center',
             horizontalalignment='center')
#axs.set_xticks(np.arange(1.5, len(exp_2cat)+6, 3.5))
axs.set_xticks(np.arange(0.35, len(exp_2cat)+3, 2.5))
axs.set_xticklabels(exp_2cat)
axs.set_ylabel('Normalized Similarity ($t$)')
axs.set_xlabel('Discrimination Experiment')
axs.set_xticklabels(['Object Invariance (OI)', 'Invariant Representation (IR)'])
axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
fig.subplots_adjust(left=0.1, right=0.99, bottom=0.09, top=0.99)
if notch:
    fig.savefig('./figs/exp_2cat.pdf')
else:
    fig.savefig('./figs/exp_2cat_nonotch.pdf')

for i, exp in enumerate(exp_2cat):
    t = np.nanmean(comparisons[exp][area]['desc']['t'])/(
        np.nanstd(comparisons[exp][area]['desc']['t'])/
        np.sqrt(np.sum(np.isfinite(comparisons[exp][area]['desc']['t']))))
    df = np.sum(np.isfinite(comparisons[exp][area]['desc']['t']))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# oi 3.5369794277485553 29 0.0013830429411979337
# ir 4.849547871920626 35 2.5252368723756333e-05
for i, exp in enumerate(exp_2cat):
    t = np.nanmean(comparisons[exp][area]['cat']['t'])/(
        np.nanstd(comparisons[exp][area]['cat']['t'])/
        np.sqrt(np.sum(np.isfinite(comparisons[exp][area]['cat']['t']))))
    df = np.sum(np.isfinite(comparisons[exp][area]['cat']['t']))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# oi 4.607107389029861 29 7.549626353023227e-05
# ir 3.3659893614907843 35 0.00186341613842946
for i, exp in enumerate(exp_2cat):
    t = np.nanmean(comparisons[exp][area]['cat2']['t'])/(
        np.nanstd(comparisons[exp][area]['cat2']['t'])/
        np.sqrt(np.sum(np.isfinite(comparisons[exp][area]['cat2']['t']))))
    df = np.sum(np.isfinite(comparisons[exp][area]['cat2']['t']))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# oi 4.530905830018866 29 9.322030633228235e-05
# ir 4.329554981600348 35 0.00011913261163996858

for i, exp in enumerate(exp_2cat):
    desc = np.array(comparisons[exp][area]['desc']['t'])
    cat = np.array(comparisons[exp][area]['cat']['t'])
    t = np.nanmean(desc-cat)/(np.nanstd(desc-cat)/
                              np.sqrt(np.sum(np.isfinite(desc-cat))))
    df = np.sum(np.isfinite(desc-cat))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# oi 3.17820148588885 29 0.0035082344354681926
# ir 4.089516224046632 35 0.00024077968340996934
for i, exp in enumerate(exp_2cat):
    desc = np.array(comparisons[exp][area]['desc']['t'])
    cat2 = np.array(comparisons[exp][area]['cat2']['t'])
    t = np.nanmean(desc-cat2)/(np.nanstd(desc-cat2)/
                              np.sqrt(np.sum(np.isfinite(desc-cat2))))
    df = np.sum(np.isfinite(desc-cat2))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# oi 2.8438698214007383 29 0.008086497320978594
# ir 3.1545617723413155 35 0.0032949890736895943
for i, exp in enumerate(exp_2cat):
    cat = np.array(comparisons[exp][area]['cat']['t'])
    cat2 = np.array(comparisons[exp][area]['cat2']['t'])
    t = np.nanmean(cat2-cat)/(np.nanstd(cat2-cat)/
                              np.sqrt(np.sum(np.isfinite(cat2-cat))))
    df = np.sum(np.isfinite(cat2-cat))-1
    pval = stats.t.sf(t, df)*2
    print(exp, t, df, pval)
# oi 3.5809056131419914 29 0.0012315573434561672
# ir 3.577777111963037 35 0.0010380713104313863


















i=1
exp='sr'
axs.plot(i*2-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['t_diff'][0]))-0.5)*0.3,
         comparisons[exp][area]['desc']['t_diff'][0], 'ko', alpha=0.5)
axs.boxplot(np.array(comparisons[exp][area]['desc']['t_diff'])[
    np.isfinite(comparisons[exp][area]['desc']['t_diff'])], positions=[i*2-0.35],
            notch=True, bootstrap=10000, widths=0.65, showfliers=False,
            whiskerprops=dict(lw=0), capprops=dict(lw=0))
axs.set_yscale('symlog')
axs.set_xticklabels(['CRM', 'SR'])
axs.set_ylabel('Normalized Similarity')
axs.plot([-0.75, 2.1], [0, 0], 'b-')
fig.subplots_adjust(left=0.15, right=0.99, bottom=0.05, top=0.99)




fig, axs = plt.subplots(1, 1, figsize=(8,2))
# for i, area in enumerate(areas):
area = 'all'
i=0
axs.plot(i*2-0.35+(np.random.uniform(size=len(comparisons[exp][area]['desc']['d_prime']))-0.5)*0.3,
         comparisons[exp][area]['desc']['d_prime'], 'ko', alpha=0.5)
axs.boxplot(np.array(comparisons[exp][area]['desc']['d_prime'])[
    np.isfinite(comparisons[exp][area]['desc']['d_prime'])], positions=[i*2-0.35],
            notch=True, bootstrap=10000, widths=0.65, showfliers=False,
            whiskerprops=dict(lw=0), capprops=dict(lw=0))
axs.plot(i*2+0.35+(np.random.uniform(size=len(comparisons[exp][area]['cat']['d_prime']))-0.5)*0.3,
         comparisons[exp][area]['cat']['d_prime'], 'bo', alpha=0.5)
axs.boxplot(np.array(comparisons[exp][area]['cat']['d_prime'])[
    np.isfinite(comparisons[exp][area]['cat']['d_prime'])], positions=[i*2+0.35],
            notch=True, bootstrap=10000, widths=0.65, showfliers=False,
            whiskerprops=dict(lw=0), capprops=dict(lw=0))























all_ps = pickle.load(open('./crm/data/all_ps_zmean.pickle', 'rb'))

# ds = ((np.array(all_ps['all']['sim_mean'])-np.array(all_ps['all']['pmean']))/
#       np.array(all_ps['all']['pstd']))
# filt = np.isfinite(ds)
# np.corrcoef(ds[filt],(np.array(all_ps['all']['hit'])[filt]-
#                       np.array(all_ps['all']['fa'])[filt]))
    

# fig, axs = plt.subplots(1, 1, figsize=(3, 3))
# hist_range = (-1.5, 11)
# axs.hist(ds, range=hist_range)
# axs.set_xlabel(r'Neural similarity ($d^\prime$)')
# axs.set_ylabel('count')
# #axs[0].set_title(r'RVL$_{new}$-RVL$_{old}$')
# fig.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.999,
#                     wspace=0.08, hspace=0.02)
# fig.savefig('./figs/similarity_ds_nocol.pdf')
# axs.hist(ds[np.array(all_ps['all']['pval'])<0.05],
#             color='yellow', alpha=0.5, range=hist_range)
# fig.savefig('./figs/similarity_ds.pdf')

# fig, axs = plt.subplots(1, 1, figsize=(4, 3))
# axs.plot(ds[filt],(np.array(all_ps['all']['hit'])[filt]-
#                       np.array(all_ps['all']['fa'])[filt]), 'ko')
# fig.subplots_adjust(left=0.19, bottom=0.15, right=0.99, top=0.999,
#                     wspace=0.08, hspace=0.02)
# axs.set_xlabel(r'Neural similarity ($d^\prime$)')
# axs.set_ylabel(r'Recognition performance (hit $-$ FA)')
# fig.savefig('./figs/rec_ds_scatterplot.pdf')


comparisons = {}

for area in all_ps:
    comparisons[area] = {'desc': {}, 'cat': {}, 'desc-cat': {}, 'permcomp': {}}
    for comp in comparisons[area]:
        comparisons[area][comp]['ks'] = []
        comparisons[area][comp]['diff'] = []
        comparisons[area][comp]['d_prime'] = []
        comparisons[area][comp]['t'] = []
        comparisons[area][comp]['n'] = []
    for comp in ['desc', 'cat']:
        comparisons[area][comp]['mean'] = []
        comparisons[area][comp]['mean_perm'] = []
        comparisons[area][comp]['std_perm'] = []
        comparisons[area]['desc-cat']['mean_desc'] = []
    comparisons[area]['desc-cat']['mean_cat'] = []
    comparisons[area]['permcomp']['mean_desc_perm'] = []
    comparisons[area]['permcomp']['mean_cat_perm'] = []
    # for sm, spm in zip(all_ps[area]['sim_means'],
    #                    all_ps[area]['sim_perm_means']):
    for sm, spm in zip(all_ps[area]['sim_mean'],
                       all_ps[area]['pmean']):
        comparisons[area]['desc']['ks'].append(
            ks_2samp(sm['description'], spm['description']))
        comparisons[area]['desc']['mean'].append(sm['description'].mean())
        comparisons[area]['desc']['mean_perm'].append(spm['description'].mean())
        comparisons[area]['desc']['std_perm'].append(spm['description'].std())
        comparisons[area]['desc']['diff'].append(
            sm['description'].mean()-spm['description'].mean())
        comparisons[area]['desc']['d_prime'].append(
            (sm['description'].mean()-spm['description'].mean())/spm['description'].std())
        comparisons[area]['desc']['t'].append(
            np.mean(sm['description']-spm['description'])/
            (np.std(sm['description']-spm['description'])/
             np.sqrt(len(sm['description']))))
        comparisons[area]['desc']['n'].append(len(sm['description']))
        comparisons[area]['cat']['ks'].append(ks_2samp(sm['category'], spm['category']))
        comparisons[area]['cat']['mean'].append(sm['category'].mean())
        comparisons[area]['cat']['mean_perm'].append(spm['category'].mean())
        comparisons[area]['cat']['std_perm'].append(spm['category'].std())
        comparisons[area]['cat']['diff'].append(
            sm['category'].mean()-spm['category'].mean())
        comparisons[area]['cat']['d_prime'].append(
            (sm['category'].mean()-spm['category'].mean())/spm['category'].std())
        comparisons[area]['cat']['t'].append(
            np.mean(sm['category']-spm['category'])/
            (np.std(sm['category']-spm['category'])/
             np.sqrt(len(sm['category']))))
        comparisons[area]['cat']['n'].append(len(sm['category']))
        comparisons[area]['desc-cat']['ks'].append(
        ks_2samp(sm['description'], sm['category']))
        comparisons[area]['desc-cat']['mean_desc'].append(sm['description'].mean())
        comparisons[area]['desc-cat']['mean_cat'].append(sm['category'].mean())
        comparisons[area]['desc-cat']['diff'].append(
            sm['description'].mean()-sm['category'].mean())
        comparisons[area]['permcomp']['ks'].append(
            ks_2samp(spm['description'], spm['category']))
        comparisons[area]['permcomp']['mean_desc_perm'].append(spm['description'].mean())
        comparisons[area]['permcomp']['mean_cat_perm'].append(spm['category'].mean())
        comparisons[area]['permcomp']['diff'].append(
            spm['description'].mean()-sm['category'].mean())




