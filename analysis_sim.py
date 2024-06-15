
import numpy as np
import matplotlib.pyplot as plt


rng = np.random.default_rng()


n = 40
units = 2
trials = 20

fr = {}
norm = {}
for fs in ['first', 'second']:
    fr[fs] = rng.poisson(np.array([[2, 5]]).T, (units, trials))
    norm[fs] = fr[fs]/np.sqrt((fr[fs]**2).sum(0))

sim_mean = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
    norm['first'].T, norm['second'])))))





fr = {}
norm = {}
for fs in ['first', 'second']:
    #fr[fs] = rng.normal(np.array([[0, 1]]).T, scale=np.array([[0.1, 0.1]]).T, size=(units, trials))
    # fr[fs] = rng.normal(np.array([[0, 1]]).T, scale=0.1, size=(units, trials))
    fr[fs] = np.c_[rng.normal(np.array([[0, 1]]).T, scale=1, size=(units, trials)),
                   rng.normal(np.array([[1, 0]]).T, scale=1, size=(units, trials))]
    norm[fs] = fr[fs]/np.sqrt((fr[fs]**2).sum(0))

sim_mean = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
    norm['first'].T, norm['second'])))))

permutations = 10000
sm = []
spm = []
frs = []
norms = []
for s in range(n):
    print(s+1, '/', n)
    fr = {}
    norm = {}
    for fs in ['first', 'second']:
        #fr[fs] = rng.normal(np.array([[0, 1]]).T, scale=np.array([[0.1, 0.1]]).T, size=(units, trials))
        # fr[fs] = rng.normal(np.array([[0, 1]]).T, scale=0.1, size=(units, trials))
        fr[fs] = np.c_[rng.normal(np.array([[0, 1]]).T, scale=1, size=(units, trials)),
                       rng.normal(np.array([[1, 0]]).T, scale=1, size=(units, trials))]
        norm[fs] = fr[fs]/np.sqrt((fr[fs]**2).sum(0))
    frs.append(fr)
    norms.append(norm)
    sim_mean = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
        norm['first'].T, norm['second'])))))
    sm.append(sim_mean)
    sim_perm_mean = np.ones(permutations)*np.nan
    for p in range(permutations):
        second_perm = rng.permutation(norm['second'], 1)
        sim_perm_mean[p] = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
            norm['first'].T, second_perm)))))
    spm.append(sim_perm_mean.mean())

sm = np.array(sm)
spm = np.array(spm)
t_diff = (sm-spm)/(np.nanstd(sm-spm)/np.sqrt(len(sm)))


sm45 = []
spm45 = []
frs45 = []
norms45 = []
for s in range(n):
    print(s+1, '/', n)
    fr = {}
    norm = {}
    for fs in ['first', 'second']:
        #fr[fs] = rng.normal(np.array([[0, 1]]).T, scale=np.array([[0.1, 0.1]]).T, size=(units, trials))
        # fr[fs] = rng.normal(np.array([[0, 1]]).T, scale=0.1, size=(units, trials))
        fr[fs] = np.c_[rng.normal(np.array([[0, 1]]).T, scale=1, size=(units, trials)),
                       rng.normal(np.array([[0, -1]]).T, scale=1, size=(units, trials))]
        norm[fs] = fr[fs]/np.sqrt((fr[fs]**2).sum(0))
    frs45.append(fr)
    norms45.append(norm)
    sim_mean = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
        norm['first'].T, norm['second'])))))
    sm45.append(sim_mean)
    sim_perm_mean = np.ones(permutations)*np.nan
    for p in range(permutations):
        second_perm = rng.permutation(norm['second'], 1)
        sim_perm_mean[p] = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
            norm['first'].T, second_perm)))))
    spm45.append(sim_perm_mean.mean())

sm45 = np.array(sm45)
spm45 = np.array(spm45)
t_diff45 = (sm45-spm45)/(np.nanstd(sm45-spm45)/np.sqrt(len(sm45)))

fig, axs = plt.subplots(1, 1, figsize=(5,5))
for i in [0]:  # enumerate(exp_0cat):
    axs.plot(i*2-0.35+(np.random.uniform(size=len(t_diff))-0.5)*0.3,
             t_diff, 'ko', alpha=0.5)
    axs.boxplot(np.array(t_diff), positions=[i*2-0.35],
                notch=True, bootstrap=10000, widths=0.65, showfliers=False,
                # notch=False, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))
for i in [1]:  # enumerate(exp_0cat):
    axs.plot(i*2-0.35+(np.random.uniform(size=len(t_diff45))-0.5)*0.3,
             t_diff45, 'ko', alpha=0.5)
    axs.boxplot(np.array(t_diff45), positions=[i*2-0.35],
                notch=True, bootstrap=10000, widths=0.65, showfliers=False,
                # notch=False, bootstrap=10000, widths=0.65, showfliers=False,
                whiskerprops=dict(lw=0), capprops=dict(lw=0))
axs.set_yscale('symlog')
#axs.set_xticklabels(exp_0cat)
axs.set_ylabel('Normalized Similarity')
axs.plot(axs.get_xlim(), [0, 0], color='red', ls='-', lw=2)
fig.subplots_adjust(left=0.15, right=0.99, bottom=0.05, top=0.99)
#fig.savefig('./figs/exp_0cat.pdf')


    
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
    axs[i].text(0.15, 0.1, exp,
                horizontalalignment='center', size='large',
                weight=1000, verticalalignment='center',
                transform=axs[i].transAxes)
    axs[i].text(0.2, 0.2, r'$r=$'+str(
        np.round(np.corrcoef(
            (np.array(all_ps[exp][area]['hit'])[filt]-
             np.array(all_ps[exp][area]['fa'])[filt]),
            comparisons[exp][area]['desc']['t_diff'][0][filt])[0,1], 2)),
                horizontalalignment='center', verticalalignment='center',
                transform=axs[i].transAxes)
axs[0].set_ylabel('Normalized Similarity')
fig.subplots_adjust(left=0.15, right=0.99, bottom=0.05, top=0.99)
fig.savefig('./figs/exp_0cat_recogperf.pdf')








plt.figure()
plt.plot([np.zeros(len(norm['first'][0])), norm['first'][0]],
         [np.zeros(len(norm['first'][1])), norm['first'][1]], 'r-', alpha=0.2)
plt.plot([np.zeros(len(norm['second'][0])), norm['second'][0]],
         [np.zeros(len(norm['second'][1])), norm['second'][1]], 'b-', alpha=0.2)

plt.plot([0, np.mean(norm['first'][0][:trials])],
         [0, np.mean(norm['first'][1][:trials])], 'r-', lw=5)
plt.plot([0, np.mean(norm['first'][0][trials:])],
         [0, np.mean(norm['first'][1][trials:])], 'r-', lw=5)
plt.plot([0, np.mean(norm['second'][0][:trials])],
         [0, np.mean(norm['second'][1][:trials])], 'b-', lw=5)
plt.plot([0, np.mean(norm['second'][0][trials:])],
         [0, np.mean(norm['second'][1][trials:])], 'b-', lw=5)


permutations = 10000

sim_perm_mean = np.ones(permutations)*np.nan
for p in range(permutations):
    second_perm = rng.permutation(norm['second'], 1)
    sim_perm_mean[p] = np.tanh(np.nanmean(np.arctanh(np.diag(np.dot(
        norm['first'].T, second_perm)))))









norm1 = fr1/np.sqrt((fr1**2).sum(0))
# (norm**2).sum(0)
norm2 = fr2/np.sqrt((fr2**2).sum(0))





# tst = rng.poisson(np.atleast_2d(np.arange(50)).T, (50, trials))

# tstnorm = tst/np.sqrt((tst**2).sum(0))
