import numpy as np
import gzip
from glob import glob
from numpy.lib.recfunctions import append_fields
from scipy.stats import linregress
from collections import defaultdict
import sys
# sys.path.append('/home/ctw/src')
# import CELEX
# import ezodf
import pandas as pd



rhino_root = '/home/ctw/fusemounts/rhino'
data_dir = rhino_root+'/scratch/josh/BniData/Subjects'
open_funcs = {'gz': gzip.open, 'txt': open}



dat_files = {}
dat_files['cf'] = sorted(glob(data_dir+'/s???/data/s*e*cf/CFLog_*'))

exp_dict = {}
ev_raw_dtypes = [('time', int), ('code', int), ('param', int)]


stim_tab = pd.read_table('./data/catmapExt_CF.csv')

evs = None
bad_exps = []
for datfile in dat_files['cf']:
    subj = datfile.split('/')[-4]
    expstr = datfile.split('/')[-2]
    if expstr in bad_exps:
        continue
    print(datfile)
    if (evs is not None) and (expstr in evs['expstr']):
        # in some cases log files are saved both in compressed and
        # uncompressed formats, we only want to read it once:
        continue
    if int(subj[1:]) != int(expstr.split('e')[0][1:]):
        raise ValueError('Unsupported subj: '+datfile)
    exp, subj_check, sess = datfile.split('/')[-1].split('.')[0].split('_')
    if int(subj[1:]) != int(subj_check[1:]):
        raise ValueError('Unsupported subj: '+datfile)
    if exp.upper() == 'CFLOG':
        exp = 'cf'
    else:
        raise ValueError('Unsupported exp: '+exp)
    expnum = int(expstr.split('e')[1][:len(expstr.split('e')[1])-len(exp)])
    sess = int(sess)
    analysis_dir = os.path.dirname(datfile).replace(
        'data','analysis')
    apropsfile = analysis_dir + '/analysisProps.txt'
    # if not os.path.exists(apropsfile):
    #     print(apropsfile, "DOES NOT EXIST!")
    #     continue
    keys_reversed = False
    if os.path.exists(apropsfile):
        with open(apropsfile, 'rt') as f:
            for line in f:
                if line == 'jKeyIsAnimate = false':
                    keys_reversed = True
    game_dict_parse = False
    game_dict_used = False # to ignore empty lines right after section header
    ev_parse = False
    ev_parse_used = False
    stim_dict = {}
    key_dict = {}
    sessevs_raw = []
    with open_funcs[datfile.split('.')[-1]](datfile, 'rt') as f:
        for line in f:
            if 'The Game Dictionary' in line:
                game_dict_parse = True
                continue
            elif 'LOCALTIME EVENT PARAMETER' in line:
                ev_parse = True
                continue
            if line.strip() == '':
                if game_dict_used:
                    game_dict_parse = False
                if ev_parse_used:
                    ev_parse = False
                continue
            if game_dict_parse:
                game_dict_used = True
                if 'Image' in line:
                    stim_params = line.split()
                    code = np.int(stim_params[0])
                    tmp_dict = {}
                    tmp_dict['target_file'] = stim_params[2]
                    stimnum = int(stim_params[2].split('.')[0])
                    tmp_dict['target_num'] = stimnum
                    tmp_dict['target_name'] = stim_tab.loc[
                        stim_tab['Number']==tmp_dict['target_num']][
                            'Name'].values[0]
                    tmp_dict['target_cat'] = stim_tab.loc[
                        stim_tab['Number']==tmp_dict['target_num']][
                            'Category'].values[0]
                    tmp_dict['target_animate'] = stim_tab.loc[
                        stim_tab['Number']==tmp_dict['target_num']][
                            'Is It Living? '].values[0] == 'T'
                    if code in stim_dict:
                        assert stim_dict[code] == tmp_dict, \
                            'stim_dict[code] == tmp_dict'
                    else:
                        stim_dict[code] = tmp_dict.copy()
                elif 'Flanker ' in line:
                    # space in 'Flanker ' is important to disambiguate
                    # from 'Conflict Flankers'
                    stim_params = line.split()
                    code = np.int(stim_params[0])
                    tmp_dict = {}
                    tmp_dict['flanker_file'] = stim_params[2]
                    stimnum = int(stim_params[2].split('.')[0])
                    tmp_dict['flanker_num'] = stimnum
                    tmp_dict['flanker_name'] = stim_tab.loc[
                        stim_tab['Number']==tmp_dict['flanker_num']][
                            'Name'].values[0]
                    tmp_dict['flanker_cat'] = stim_tab.loc[
                        stim_tab['Number']==tmp_dict['flanker_num']][
                            'Category'].values[0]
                    tmp_dict['flanker_animate'] = stim_tab.loc[
                        stim_tab['Number']==tmp_dict['flanker_num']][
                            'Is It Living? '].values[0] == 'T'
                    if code in stim_dict:
                        assert stim_dict[code] == tmp_dict, \
                            'stim_dict[code] == tmp_dict'
                    else:
                        stim_dict[code] = tmp_dict.copy()
                elif 'key' in line.lower():
                    key_params = line.split()
                    code = np.int(key_params[0])
                    tmp_dict = {}
                    if 'anim' in line:
                        if keys_reversed:
                            tmp_dict['response'] = 'inanimate'
                        else:
                            tmp_dict['response'] = 'animate'
                    elif 'inan' in line:
                        if keys_reversed:
                            tmp_dict['response'] = 'animate'
                        else:
                            tmp_dict['response'] = 'inanimate'
                    elif 'Other' in line:
                        tmp_dict['response'] = 'other'
                    else:
                        raise ValueError('Unknown key parameter: '+line)
                    if code in key_dict:
                        assert key_dict[code] == tmp_dict, \
                            'key_dict[code] == tmp_dict'
                    else:
                        key_dict[code] = tmp_dict.copy()
                elif 'Synchronize Server' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {}
                    if 'Start' in line:
                        tmp_dict['evtype'] = 'synch_start'
                    elif 'End' in line:
                        tmp_dict['evtype'] = 'synch_end'
                    else:
                        raise ValueError('Unknown exp parameter: '+line)
                    if code in exp_dict:
                        assert exp_dict[code] == tmp_dict, \
                            'exp_dict[code] == tmp_dict'
                    else:
                        exp_dict[code] = tmp_dict.copy()
                elif 'Phase' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {}
                    if 'Start' in line:
                        start_phase_code = code
                        tmp_dict['evtype'] = 'phase_start'
                    else:
                        raise ValueError('Unknown exp parameter: '+line)
                    if code in exp_dict:
                        assert exp_dict[code] == tmp_dict, \
                            'exp_dict[code] == tmp_dict'
                    else:
                        exp_dict[code] = tmp_dict.copy()
                elif 'Conflict' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {}
                    if 'Started' in line:
                        tmp_dict['evtype'] = 'cf_start'
                    elif 'Stopped' in line:
                        tmp_dict['evtype'] = 'cf_end'
                    else:
                        raise ValueError('Unknown exp parameter: '+line)
                    if code in exp_dict:
                        assert exp_dict[code] == tmp_dict, \
                            'exp_dict[code] == tmp_dict'
                    else:
                        exp_dict[code] = tmp_dict.copy()
                elif 'fixation' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {'evtype': 'fixatin'}
                    if code in exp_dict:
                        assert exp_dict[code] == tmp_dict, \
                            'exp_dict[code] == tmp_dict'
                    else:
                        exp_dict[code] = tmp_dict.copy()
                elif 'Initial' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {'evtype': 'init_wait'}
                    if code in exp_dict:
                        assert exp_dict[code] == tmp_dict, \
                            'exp_dict[code] == tmp_dict'
                    else:
                        exp_dict[code] = tmp_dict.copy()
                else:
                    raise ValueError('Unknown Game Dictionary parameter: ' +
                                     line)
            #
            #
            # parsing events:
            #
            elif ev_parse:
                ev_parse_used = True
                exptime, expcode, exppar = line.split()
                exptime = np.int(exptime)
                expcode = np.int(expcode)
                exppar = np.int(exppar)
                sessevs_raw.append((exptime, expcode, exppar))
    sessevs_raw_rec = np.array(sessevs_raw, dtype=ev_raw_dtypes)
    times = []
    time_thresh = 5
    for i,row in enumerate(sessevs_raw_rec):
        if row['code'] != 1:
            continue
        assert sessevs_raw_rec[i-1]['code'] == 0, \
            'sessevs_raw_rec[i-1][\'code\'] == 0'
        if (row['time'] - sessevs_raw_rec[i-1]['time']) > time_thresh:
            continue
        # take mean of send and receive timestamp as exptime to
        # convert to neurotime:
        times.append((np.mean([sessevs_raw_rec[i-1]['time'], row['time']]),
                      row['param']))
    times_rec = np.array(times, dtype=[('exptime', float), ('neurotime', int)])
    if np.std(times_rec['neurotime']) == 0:
        # there are cases where the clock seems to not have worked and
        # all neurotimes are zero, preventing alignment
        continue
    times_mean = []
    for i,timediff in enumerate(np.diff(times_rec['exptime'])):
        if timediff < time_thresh:
            times_mean.append(
                (np.mean([times_rec['exptime'][i], times_rec['exptime'][i+1]]),
                 np.mean([times_rec['neurotime'][i],
                          times_rec['neurotime'][i+1]])))
    times_mean_rec = np.array(times_mean, dtype=[('exptime', float),
                                                 ('neurotime', float)])
    # there should be a pair of synch pulses at the beginning and end:
    if len(times_mean_rec) < 2:
        # one or both synch pulses are missing & we can't align
        continue
    assert len(times_mean_rec) == 2, 'len(times_mean_rec) == 2'
    slope, intercept, r, p, se = linregress(times_mean_rec['exptime'],
                                            times_mean_rec['neurotime'])
    assert se == 0, 'se == 0'
    assert r > 0.999999999999999, 'r > 0.999999999999999'
    sessevs_raw_rec = append_fields(
        sessevs_raw_rec, 'neurotime', sessevs_raw_rec['time']*slope+intercept,
        dtypes=float, usemask=False, fill_value=np.nan, asrecarray=True)
    sessevs = []
    indx = 0
    if sessevs_raw_rec[indx]['code'] not in exp_dict:
        raise ValueError('Invalid start: ' + str(sessevs_raw_rec[indx]))
    if exp_dict[sessevs_raw_rec[indx]['code']]['evtype'] != 'cf_start':
        raise ValueError('cf_start code expected: ' +
                         str(sessevs_raw_rec[indx]) + ' ' +
                         exp_dict[sessevs_raw_rec[indx]['code']]['evtype'])
    indx += 1
    while indx < len(sessevs_raw_rec):
        while ((indx < len(sessevs_raw_rec)) and
               (sessevs_raw_rec[indx]['code'] != start_phase_code)):
            indx += 1
        if indx == len(sessevs_raw_rec):
            break
        if sessevs_raw_rec[indx]['param'] == 0:
            indx += 1
            continue
        if sessevs_raw_rec[indx]['param'] == 1:
            phase = 1
            ev_dict = {'phase1_time': sessevs_raw_rec[indx]['neurotime'],
                       'resp_early': False, 'resp_late': False,
                       'resp_mult': False, 'first_resp': '',
                       'first_resp_time': np.nan,
                       'keys_reversed': keys_reversed,
                       'analysis_dir': analysis_dir[len(rhino_root):]}
            indx += 1
            # sanity check to make sure we have expected number of events:
            targ_counter = 0
            flank_counter = 0
            key_count = {1: defaultdict(int), 2: defaultdict(int),
                         3: defaultdict(int)}
            key_first_time = {1: defaultdict(float), 2: defaultdict(float),
                              3: defaultdict(float)}
            key_last_time = {1: defaultdict(float), 2: defaultdict(float),
                             3: defaultdict(float)}
            while ((indx < (len(sessevs_raw_rec)-1)) and
                   (sessevs_raw_rec[indx]['code'] != start_phase_code)):
                if sessevs_raw_rec[indx]['code'] in stim_dict:
                    if targ_counter+flank_counter > 1:
                        raise ValueError('targ+flank counter > 1: ' +
                                         str(targ_counter) + ' ' +
                                         str(flank_counter))
                    # flanker should be second, so this should be true
                    # also (a bit redundant, but a sanity check that
                    # the order is as expected):
                    if flank_counter > 0:
                        raise ValueError('flank counter > 0: ' +
                                         str(flank_counter))
                    ev_dict.update(stim_dict[sessevs_raw_rec[indx]['code']])
                    if 'target_file' in stim_dict[
                            sessevs_raw_rec[indx]['code']]:
                        ev_dict['target_time'] = sessevs_raw_rec[indx][
                            'neurotime']
                        targ_counter += 1
                    elif 'flanker_file' in stim_dict[
                            sessevs_raw_rec[indx]['code']]:
                        ev_dict['flanker_time'] = sessevs_raw_rec[indx][
                            'neurotime']
                        flank_counter += 1
                    else:
                        raise ValueError('No target or flanker time')
                elif sessevs_raw_rec[indx]['code'] in key_dict:
                    ev_dict['resp_early'] = True
                    key_count[phase][key_dict[sessevs_raw_rec[indx]['code']][
                        'response']] += 1
                    if ((key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] == 0) or
                        (key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] >
                         sessevs_raw_rec[indx]['neurotime'])):
                        key_first_time[phase][
                            key_dict[sessevs_raw_rec[indx]['code']][
                                'response']] = sessevs_raw_rec[indx][
                                    'neurotime']
                        if ev_dict['first_resp'] == '':
                            ev_dict['first_resp'] = key_dict[
                                sessevs_raw_rec[indx][
                                'code']]['response']
                            ev_dict['first_resp_time'] = sessevs_raw_rec[
                                indx]['neurotime']
                    if ((key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] == 0) or
                        (key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] <
                         sessevs_raw_rec[indx]['neurotime'])):
                        key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] = sessevs_raw_rec[indx][
                                'neurotime']
                elif sessevs_raw_rec[indx]['code'] in exp_dict:
                    raise ValueError(
                        'exp code when stim code was expected -- 1: ' +
                        str(sessevs_raw_rec[indx]))
                else:
                    raise ValueError('unkown code: ' +
                                     str(sessevs_raw_rec[indx]))
                indx += 1
                if indx == len(sessevs_raw_rec):
                    raise ValueError('Incomplete trial -- 1')
            assert sessevs_raw_rec[indx]['code'] == start_phase_code, \
                'sessevs_raw_rec[indx][\'code\'] == start_phase_code'
            phase2_present = True
            if sessevs_raw_rec[indx]['param'] != 2:
                phase2_present = False
                # raise ValueError('Phase 2 expected: ' +
                #                  str(sessevs_raw_rec[indx]))
                print('Phase 2 missing for', datfile)
                ev_dict['phase2_time'] = np.nan
            if phase2_present:
                phase = 2
                ev_dict['phase2_time'] = sessevs_raw_rec[indx]['neurotime']
                indx += 1
                while ((indx < (len(sessevs_raw_rec)-1)) and
                       (sessevs_raw_rec[indx]['code'] != start_phase_code)):
                    if sessevs_raw_rec[indx]['code'] in stim_dict:
                        raise ValueError('Stim in phase 2: ' +
                                         str(sessevs_raw_rec[indx]))
                    if sessevs_raw_rec[indx]['code'] in key_dict:
                        key_count[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] += 1
                        if ((key_first_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] == 0) or
                            (key_first_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] >
                             sessevs_raw_rec[indx]['neurotime'])):
                            if ev_dict['first_resp'] == '':
                                ev_dict['first_resp'] = key_dict[
                                    sessevs_raw_rec[indx]['code']]['response']
                                ev_dict['first_resp_time'] = sessevs_raw_rec[
                                    indx]['neurotime']
                            key_first_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] = sessevs_raw_rec[
                                    indx]['neurotime']
                        if ((key_last_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] == 0) or
                            (key_last_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] <
                             sessevs_raw_rec[indx]['neurotime'])):
                            key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                                'code']]['response']] = sessevs_raw_rec[
                                    indx]['neurotime']
                    elif sessevs_raw_rec[indx]['code'] in exp_dict:
                        print('exp code when key code was expected -- 2',
                              sessevs_raw_rec[indx])
                        # raise ValueError(
                        #     'exp code when key code was expected -- 2: ' +
                        #     str(sessevs_raw_rec[indx]))
                    else:
                        raise ValueError('unkown code: ' +
                                         str(sessevs_raw_rec[indx]))
                    indx += 1
                if indx == len(sessevs_raw_rec):
                    raise ValueError('Incomplete trial -- 2')
            if sessevs_raw_rec[indx]['code'] != start_phase_code:
                print('Start phase code expected', sessevs_raw_rec[indx])
            # assert sessevs_raw_rec[indx]['code'] == start_phase_code, \
            #     'sessevs_raw_rec[indx][\'code\'] == start_phase_code'
            phase3_present = True
            if sessevs_raw_rec[indx]['param'] != 3:
                # raise ValueError('Phase 3 expected: ' +
                #                 str(sessevs_raw_rec[indx]))
                print('Phase 3 expected', sessevs_raw_rec[indx])
                phase3_present = False
                ev_dict['phase3_time'] = np.nan
            if phase3_present:
                phase = 3
                ev_dict['phase3_time'] = sessevs_raw_rec[indx]['neurotime']
                indx += 1
                while ((indx < (len(sessevs_raw_rec)-1)) and
                       (sessevs_raw_rec[indx]['code'] != start_phase_code)):
                    if sessevs_raw_rec[indx]['code'] in stim_dict:
                        raise ValueError('Stim in phase 3: ' +
                                         str(sessevs_raw_rec[indx]))
                    elif sessevs_raw_rec[indx]['code'] in key_dict:
                        key_count[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] += 1
                        ev_dict['resp_late'] = True
                        if ((key_first_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] == 0) or
                            (key_first_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] > sessevs_raw_rec[
                                    indx]['neurotime'])):
                            if ev_dict['first_resp'] == '':
                                ev_dict['first_resp'] = key_dict[
                                    sessevs_raw_rec[indx][
                                    'code']]['response']
                                ev_dict['first_resp_time'] = sessevs_raw_rec[
                                    indx]['neurotime']
                            key_first_time[phase][key_dict[
                                sessevs_raw_rec[indx]['code']][
                                    'response']] = sessevs_raw_rec[indx][
                                        'neurotime']
                        if ((key_last_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] == 0) or
                            (key_last_time[phase][key_dict[sessevs_raw_rec[
                                indx]['code']]['response']] <
                             sessevs_raw_rec[indx]['neurotime'])):
                            key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                                'code']]['response']] = sessevs_raw_rec[
                                    indx]['neurotime']
                        indx += 1
                    elif sessevs_raw_rec[indx]['code'] in exp_dict:
                        if sessevs_raw_rec[indx]['code'] != start_phase_code:
                            indx += 1
                    else:
                        raise ValueError('unkown code: ' +
                                         str(sessevs_raw_rec[indx]))
            if not 'target_time' in ev_dict:
                continue
            if not 'flanker_time' in ev_dict:
                continue
            resp_count_all = 0
            for phase in key_count:
                for resp in [key_dict[r]['response'] for r in key_dict]:
                    if phase <= 2:
                        resp_count_all += key_count[phase][resp]
                    ev_dict['resp_count_'+resp+'_'+str(phase)] = key_count[
                        phase][resp]
                    firsttime = key_first_time[phase][resp]
                    if firsttime == 0:
                        firsttime = np.nan
                    lasttime = key_last_time[phase][resp]
                    if lasttime == 0:
                        lasttime = np.nan
                    ev_dict['resp_firsttime_'+resp+'_'+str(phase)] = firsttime
                    ev_dict['resp_lasttime_'+resp+'_'+str(phase)] = lasttime
            if resp_count_all > 1:
                ev_dict['resp_mult'] = True
            sessevs.append(
                (exp, subj, expnum, sess, expstr, ev_dict['target_time'],
                 ev_dict['flanker_time'], ev_dict['phase1_time'],
                 ev_dict['phase2_time'], ev_dict['phase3_time'], 'cf_trial',
                 ev_dict['target_file'], ev_dict['target_num'],
                 ev_dict['target_name'], ev_dict['target_cat'],
                 ev_dict['target_animate'],
                 ev_dict['flanker_file'], ev_dict['flanker_num'],
                 ev_dict['flanker_name'], ev_dict['flanker_cat'],
                 ev_dict['flanker_animate'],
                 ev_dict['flanker_animate'] == ev_dict['target_animate'],
                 ev_dict['resp_early'], ev_dict['resp_late'],
                 ev_dict['resp_mult'], ev_dict['first_resp'],
                 ev_dict['first_resp_time'], ev_dict['keys_reversed'],
                 ev_dict['analysis_dir'], ev_dict['resp_count_animate_1'],
                 ev_dict['resp_count_animate_2'],
                 ev_dict['resp_count_animate_3'],
                 ev_dict['resp_count_inanimate_1'],
                 ev_dict['resp_count_inanimate_2'],
                 ev_dict['resp_count_inanimate_3'],
                 ev_dict['resp_count_other_1'], ev_dict['resp_count_other_2'],
                 ev_dict['resp_count_other_3'],
                 ev_dict['resp_firsttime_animate_1'],
                 ev_dict['resp_firsttime_animate_2'],
                 ev_dict['resp_firsttime_animate_3'],
                 ev_dict['resp_firsttime_inanimate_1'],
                 ev_dict['resp_firsttime_inanimate_2'],
                 ev_dict['resp_firsttime_inanimate_3'],
                 ev_dict['resp_firsttime_other_1'],
                 ev_dict['resp_firsttime_other_2'],
                 ev_dict['resp_firsttime_other_3'],
                 ev_dict['resp_lasttime_animate_1'],
                 ev_dict['resp_lasttime_animate_2'],
                 ev_dict['resp_lasttime_animate_3'],
                 ev_dict['resp_lasttime_inanimate_1'],
                 ev_dict['resp_lasttime_inanimate_2'],
                 ev_dict['resp_lasttime_inanimate_3'],
                 ev_dict['resp_lasttime_other_1'],
                 ev_dict['resp_lasttime_other_2'],
                 ev_dict['resp_lasttime_other_3']))
    if indx < len(sessevs_raw_rec):
        raise ValueError('Incomplete trial -- end')
    ev_dtypes = [('exp', 'U8'), ('subject', 'U8'), ('expnum', int),
                 ('session', int), ('expstr', 'U16'), ('target_time', float),
                 ('flanker_time', float), ('phase1_time', float),
                 ('phase2_time', float), ('phase3_time', float),
                 ('evtype', 'U16'), ('target_file', 'U32'), ('target_num', int),
                 ('target_name', 'U32'), ('target_cat', 'U32'),
                 ('target_animate', bool), ('flanker_file', 'U32'),
                 ('flanker_num', int), ('flanker_name', 'U32'),
                 ('flanker_cat', 'U32'), ('flanker_animate', bool), 
                 ('flanker_match', bool), ('resp_early', bool),
                 ('resp_late', bool), ('resp_mult', bool),
                 ('first_resp', 'U16'),
                 ('first_resp_time', float), ('keys_reversed', bool),
                 ('analysis_dir', 'U128'), ('resp_count_animate_1', int),
                 ('resp_count_animate_2', int), ('resp_count_animate_3', int),
                 ('resp_count_inanimate_1', int),
                 ('resp_count_inanimate_2', int),
                 ('resp_count_inanimate_3', int), ('resp_count_other_1', int),
                 ('resp_count_other_2', int), ('resp_count_other_3', int),
                 ('resp_firsttime_animate_1', float),
                 ('resp_firsttime_animate_2', float),
                 ('resp_firsttime_animate_3', float),
                 ('resp_firsttime_inanimate_1', float),
                 ('resp_firsttime_inanimate_2', float),
                 ('resp_firsttime_inanimate_3', float),
                 ('resp_firsttime_other_1', float),
                 ('resp_firsttime_other_2', float),
                 ('resp_firsttime_other_3', float),
                 ('resp_lasttime_animate_1', float),
                 ('resp_lasttime_animate_2', float),
                 ('resp_lasttime_animate_3', float),
                 ('resp_lasttime_inanimate_1', float),
                 ('resp_lasttime_inanimate_2', float),
                 ('resp_lasttime_inanimate_3', float),
                 ('resp_lasttime_other_1', float),
                 ('resp_lasttime_other_2', float),
                 ('resp_lasttime_other_3', float)]
    sessevs_rec = np.array(sessevs, ev_dtypes)
    if evs is None:
        evs = sessevs_rec
    else:
        evs = np.r_[evs, sessevs_rec]



np.savez_compressed('./data/bni_evs_cf.npz', evs=evs)
# np.save('./data/bni_evs_crm2.npy', evs, allow_pickle=False)



# tst1 = pd.read_csv('/home/ctw/fusemounts/rhino/scratch/josh/BniData/Subjects/s030/analysis/s30e5cr/sliceInfo.txt', sep='\t')
# tstevs1 = evs[evs['expstr'] == 's30e5cr']
# tst2 = pd.read_csv('/home/ctw/fusemounts/rhino/scratch/josh/BniData/Subjects/s012/analysis/s12e10cr/crmTrials.txt', sep='\t')
# tstevs2 = evs[evs['expstr'] == 's12e10cr']





# tmp1 = slice_info['firstResp']/10000
# tmp2 = evs_sess['first_resp_time']-evs_sess['stim_time']
# # tmp2 = evs_sess['first_resp_time']-evs_sess['phase1_time']
# tmp2[~np.isfinite(tmp2)] = -1000
# assert np.allclose(tmp1[(~slice_info['keyEarly']) &
#                         (slice_info['firstResp'] > 0)],
#                    tmp2[(~slice_info['keyEarly']) &
#                         (slice_info['firstResp'] > 0)],
#                    rtol=0.005, atol=0.0005), \
#                    """np.allclose(tmp1[(~slice_info['keyEarly']) &
#                    (slice_info['firstResp'] > 0)],
#                    tmp2[(~slice_info['keyEarly']) &
#                    (slice_info['firstResp'] > 0)],
#                    rtol=0.005, atol=0.0005)"""




# tsta1 = slice_info['firstResp'][(~slice_info['keyEarly']) &
#                         (slice_info['firstResp'] > 0)]-slice_info['endTime'][(~slice_info['keyEarly']) &
#                         (slice_info['firstResp'] > 0)]
# tsta2 = evs_sess['first_resp_time'][(~slice_info['keyEarly']) &
#                         (slice_info['firstResp'] > 0)]-evs_sess['stim_time'][(~slice_info['keyEarly']) &
#                         (slice_info['firstResp'] > 0)]


# tsta1 = slice_info['firstResp']-slice_info['stTime']
# tsta1 = slice_info['endTime']-slice_info['stTime']
# tsta1 = slice_info['firstResp']-slice_info['endTime']
# tsta2 = evs_sess['first_resp_time']-evs_sess['stim_time']


