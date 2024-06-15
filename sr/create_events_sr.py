import numpy as np
import gzip
from glob import glob
from numpy.lib.recfunctions import append_fields
from scipy.stats import linregress
from collections import defaultdict
import sys
sys.path.append('/home/ctw/src')
import CELEX
import ezodf
import pandas as pd



rhino_root = '/home/ctw/fusemounts/rhino'
data_dir = rhino_root+'/scratch/josh/BniData/Subjects'
open_funcs = {'gz': gzip.open, 'txt': open}


stimwords_dict = {}


CELEX.locl = '/home/ctw/src/wordpools/trunk/stats/celex2/'


# Phases:
# 0 wait for user
# 1 prep
# 2 sound
# 3 response
# 4 delay



def add_stimword(word):
    if word in stimwords_dict:
        return
    clx = CELEX.CELEX('english')
    celex_wf = CELEX.Wordforms(clx)
    totfreq = 0
    for wf in clx.ls2ws(word):
        totfreq += celex_wf.ws2wf(wf)  
    stimwords_dict[word] = {'celex_wf': totfreq}


dat_files = {}
dat_files['sr'] = sorted(glob(data_dir+'/s???/data/s*e*sr/SRLog_*'))

for df in dat_files['sr']:
    if df[:-4]+'_fixed.txt' in dat_files['sr']:
        dat_files['sr'].remove(df)

exp_dict = {}
ev_raw_dtypes = [('time', int), ('code', int), ('param', int)]

evs = None
bad_exps = []  #['s14e17cr']
for datfile in dat_files['sr']:
    subj = datfile.split('/')[-4]
    expstr = datfile.split('/')[-2]
    if expstr in bad_exps:
        continue
    print(datfile)
    if (evs is not None) and (expstr in evs['expstr']):
        # in some cases log files are saved both in compressed and
        # uncompressed formats, we only want to read it once:
        if expstr != 's55e4sr':
            # in 's55e4sr' two sessions were run back to back using
            # the same expstr ... we don't want to exclude that 2nd
            # session:
            continue
    if int(subj[1:]) != int(expstr.split('e')[0][1:]):
        raise ValueError('Unsupported subj: '+datfile)
    if 'fixed' in datfile:
        exp, subj_check, sess, _ = datfile.split('/')[-1].split(
            '.')[0].split('_')
    else:
        exp, subj_check, sess = datfile.split('/')[-1].split('.')[0].split('_')
    if int(subj[1:]) != int(subj_check[1:]):
        raise ValueError('Unsupported subj: '+datfile)
    if exp.upper() == 'SRLOG':
        exp = 'sr'
    else:
        raise ValueError('Unsupported exp: '+exp)
    expnum = int(expstr.split('e')[1][:len(expstr.split('e')[1])-len(exp)])
    sess = int(sess)
    analysis_dir = os.path.dirname(datfile).replace(
        'data', 'analysis')
    apropsfile = analysis_dir + '/analysisProps.txt'
    # if not os.path.exists(apropsfile):
    #     print(apropsfile, "DOES NOT EXIST!")
    #     continue
    keys_reversed = False
    if os.path.exists(apropsfile):
        with open(apropsfile, 'rt') as f:
            for line in f:
                if line == 'jKeyIsNew = false':
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
            # print(line)
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
                if 'Sound' in line:
                    stim_params = line.split()
                    code = np.int(stim_params[0])
                    tmp_dict = {}
                    tmp_dict['stim_file'] = stim_params[2]
                    tmp_dict['stim_word'] = stim_params[2].split('_')[0][:-1]
                    tmp_dict['stim_version'] = stim_params[2].split('_')[0][-1]
                    tmp_dict['stim_mf'] = stim_params[2].split('_')[1][0]
                    tmp_dict['stim_mf_version'] = stim_params[2].split('_')[
                        1][1]
                    add_stimword(tmp_dict['stim_word'])
                    tmp_dict['stim_freq'] = stimwords_dict[
                        tmp_dict['stim_word']]['celex_wf']
                    if code in stim_dict:
                        assert stim_dict[code] == tmp_dict, \
                            'stim_dict[code] == tmp_dict'
                    else:
                        stim_dict[code] = tmp_dict.copy()
                elif 'key' in line:
                    key_params = line.split()
                    code = np.int(key_params[0])
                    tmp_dict = {}
                    if 'j' in line:
                        if keys_reversed:
                            tmp_dict['response'] = 'new'
                        else:
                            tmp_dict['response'] = 'old'
                    elif 'f' in line:
                        if keys_reversed:
                            tmp_dict['response'] = 'old'
                        else:
                            tmp_dict['response'] = 'new'
                    else:
                        raise ValueError('Unknown key parameter: '+line)
                    if code in key_dict:
                        assert key_dict[code] == tmp_dict, \
                            'key_dict[code] == tmp_dict'
                    else:
                        key_dict[code] = tmp_dict.copy()
                elif 'Key' in line:
                    key_params = line.split()
                    code = np.int(key_params[0])
                    tmp_dict = {}
                    if 'Other' in line:
                        tmp_dict['response'] = 'other'
                        if code in key_dict:
                            assert key_dict[code] == tmp_dict, \
                                'key_dict[code] == tmp_dict'
                        else:
                            key_dict[code] = tmp_dict.copy()
                    elif 'Prompt' in line:
                        tmp_dict['evtype'] = 'prompt'
                        if code in exp_dict:
                            assert exp_dict[code] == tmp_dict, \
                                'exp_dict[code] == tmp_dict'
                        else:
                            exp_dict[code] = tmp_dict.copy()
                    else:
                        raise ValueError('Unknown key parameter: '+line)
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
                elif 'Spoken' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {}
                    if 'Started' in line:
                        tmp_dict['evtype'] = 'recog_start'
                    elif 'Stopped' in line:
                        tmp_dict['evtype'] = 'recog_end'
                    else:
                        raise ValueError('Unknown exp parameter: '+line)
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
                elif 'Fixation' in line:
                    exp_params = line.split()
                    code = np.int(exp_params[0])
                    tmp_dict = {'evtype': 'fixation'}
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
    assert exp_dict[1]['evtype'] == 'synch_end', 'synch_end code'
    for i,row in enumerate(sessevs_raw_rec):
        if row['code'] != 1:  # synch_end code as tested by above assertion
            continue
        assert sessevs_raw_rec[i-1]['code'] == 0, \
            'sessevs_raw_rec[i-1][\'code\'] == 0'
        assert exp_dict[sessevs_raw_rec[i-1]['code']][
            'evtype'] == 'synch_start', 'synch_start before end'
        if (row['time'] - sessevs_raw_rec[i-1]['time']) >= time_thresh:
            continue
        # take mean of send and receive timestamp as exptime to
        # convert to neurotime:
        times.append((np.mean([sessevs_raw_rec[i-1]['time'], row['time']]),
                      row['param']))
    times_rec = np.array(times, dtype=[('exptime', float), ('neurotime', int)])
    times_mean = []
    for i,timediff in enumerate(np.diff(times_rec['exptime'])):
        if timediff < time_thresh:
            times_mean.append(
                (np.mean([times_rec['exptime'][i], times_rec['exptime'][i+1]]),
                 np.mean([times_rec['neurotime'][i],
                          times_rec['neurotime'][i+1]])))
    times_mean_rec = np.array(times_mean, dtype=[('exptime', float),
                                                 ('neurotime', float)])
    assert len(times_mean_rec) == 2, 'len(times_mean_rec) == 2'
    slope, intercept, r, p, se = linregress(times_mean_rec['exptime'],
                                            times_mean_rec['neurotime'])
    assert se == 0, 'se == 0'
    assert r > 0.999999, 'r > 0.999999'
    sessevs_raw_rec = append_fields(
        sessevs_raw_rec, 'neurotime', sessevs_raw_rec['time']*slope+intercept,
        dtypes=float, usemask=False, fill_value=np.nan, asrecarray=True)
    sessevs = []
    indx = 0
    if sessevs_raw_rec[indx]['code'] not in exp_dict:
        raise ValueError('Invalid start: ' + str(sessevs_raw_rec[indx]))
    if exp_dict[sessevs_raw_rec[indx]['code']]['evtype'] != 'recog_start':
        raise ValueError('recognition started code expected: ' +
                         str(sessevs_raw_rec[indx]))
    indx += 1
    while indx < len(sessevs_raw_rec):
        while ((indx < len(sessevs_raw_rec)) and
               (sessevs_raw_rec[indx]['code'] != start_phase_code)):
            indx += 1
        if indx == len(sessevs_raw_rec):
            break
        if sessevs_raw_rec[indx]['param'] == 0:
            # if sessevs_raw_rec[indx]['param'] != 1:
            indx += 1
            continue
        #if sessevs_raw_rec[indx]['param'] == 2:   
        # if ((sessevs_raw_rec[indx]['param'] == 1) or
        #     (sessevs_raw_rec[indx]['param'] == 2)) :
        stim_counter = 0
        key_count = {1: defaultdict(int), 2: defaultdict(int),
                     3: defaultdict(int), 4: defaultdict(int)}
        key_first_time = {1: defaultdict(float), 2: defaultdict(float),
                          3: defaultdict(float), 4: defaultdict(float)}
        key_last_time = {1: defaultdict(float), 2: defaultdict(float),
                         3: defaultdict(float), 4: defaultdict(float)}
        ev_dict = {'phase1_time': np.nan, 'phase2_time': np.nan,
                   'phase3_time': np.nan, 'phase4_time': np.nan,
                   'resp_early': False, 'resp_late': False,
                   'resp_mult': False, 'first_resp': '',
                   'first_resp_time': np.nan,
                   'keys_reversed': keys_reversed,
                   'analysis_dir': analysis_dir[len(rhino_root):]}
 
        if sessevs_raw_rec[indx]['param'] == 1:
            phase = 1
            indx += 1
            ev_dict['phase1_time'] = sessevs_raw_rec[indx]['neurotime']
            # sanity check to make sure we have expected number of events:
            while ((indx < (len(sessevs_raw_rec)-1)) and
                   (sessevs_raw_rec[indx]['code'] != start_phase_code)):
                if sessevs_raw_rec[indx]['code'] in stim_dict:
                    raise ValueError('Stim in phase 1: ' +
                                     str(sessevs_raw_rec[indx]))
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
                            ev_dict['first_resp'] = key_dict[sessevs_raw_rec[
                                indx]['code']]['response']
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
            #if sessevs_raw_rec[indx]['param'] != 2:
            #    raise ValueError('Phase 2 expected: ' +
            #                     str(sessevs_raw_rec[indx]))
        if sessevs_raw_rec[indx]['param'] == 2:
            # stim_counter = 0
            phase = 2
            ev_dict['phase2_time'] = sessevs_raw_rec[indx]['neurotime']
            indx += 1
            while ((indx < (len(sessevs_raw_rec)-1)) and
                   (sessevs_raw_rec[indx]['code'] != start_phase_code)):
                # if sessevs_raw_rec[indx]['code'] in stim_dict:
                #     raise ValueError('Stim in phase 2: ' +
                #                      str(sessevs_raw_rec[indx]))
                if sessevs_raw_rec[indx]['code'] in stim_dict:
                    if stim_counter > 0:
                        raise ValueError('stim counter > 0: ' +
                                         str(stim_counter))
                    ev_dict.update(stim_dict[sessevs_raw_rec[indx]['code']])
                    ev_dict['stim_time'] = sessevs_raw_rec[indx]['neurotime']
                    stim_counter += 1
                elif sessevs_raw_rec[indx]['code'] in key_dict:
                    key_count[phase][key_dict[sessevs_raw_rec[indx]['code']][
                        'response']] += 1
                    if ((key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] == 0) or
                        (key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] >
                         sessevs_raw_rec[indx]['neurotime'])):
                        if ev_dict['first_resp'] == '':
                            ev_dict['first_resp'] = key_dict[
                                sessevs_raw_rec[indx]['code']]['response']
                            ev_dict['first_resp_time'] = sessevs_raw_rec[
                                indx]['neurotime']
                        key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] = sessevs_raw_rec[indx][
                                'neurotime']
                    if ((key_last_time[phase][key_dict[sessevs_raw_rec[
                            indx]['code']]['response']] == 0) or
                        (key_last_time[phase][key_dict[sessevs_raw_rec[
                            indx]['code']]['response']] <
                         sessevs_raw_rec[indx]['neurotime'])):
                        key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                            'code']]['response']] = sessevs_raw_rec[
                                indx]['neurotime']
                elif sessevs_raw_rec[indx]['code'] in exp_dict:
                    raise ValueError(
                        'exp code when key code was expected -- 2: ' +
                        str(sessevs_raw_rec[indx]))
                else:
                    raise ValueError('unkown code: ' +
                                     str(sessevs_raw_rec[indx]))
                indx += 1
            if indx == len(sessevs_raw_rec):
                raise ValueError('Incomplete trial -- 2')
            assert sessevs_raw_rec[indx]['code'] == start_phase_code, \
                'sessevs_raw_rec[indx][\'code\'] == start_phase_code'
            # if sessevs_raw_rec[indx]['param'] != 3:
            #     raise ValueError('Phase 3 expected: ' +
            #                     str(sessevs_raw_rec[indx]))
        if sessevs_raw_rec[indx]['param'] == 3:
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
                        key_last_time[phase][key_dict[sessevs_raw_rec[
                            indx]['code']]['response']] = sessevs_raw_rec[
                                indx]['neurotime']
                    indx += 1
                elif sessevs_raw_rec[indx]['code'] in exp_dict:
                    if sessevs_raw_rec[indx]['code'] != start_phase_code:
                        indx += 1
                else:
                    raise ValueError('unkown code: ' +
                                     str(sessevs_raw_rec[indx]))
                # indx += 1
            if indx == len(sessevs_raw_rec):
                raise ValueError('Incomplete trial -- 3')
        assert sessevs_raw_rec[indx]['code'] == start_phase_code, \
            'sessevs_raw_rec[indx][\'code\'] == start_phase_code'
        if sessevs_raw_rec[indx]['param'] != 4:
            raise ValueError('Phase 4 expected: ' +
                             str(sessevs_raw_rec[indx]))
        phase = 4
        ev_dict['phase4_time'] = sessevs_raw_rec[indx]['neurotime']
        indx += 1
        while ((indx < (len(sessevs_raw_rec)-1)) and
               (sessevs_raw_rec[indx]['code'] != start_phase_code)):
            if sessevs_raw_rec[indx]['code'] in stim_dict:
                raise ValueError('Stim in phase 4: ' +
                                 str(sessevs_raw_rec[indx]))
            elif sessevs_raw_rec[indx]['code'] in key_dict:
                key_count[phase][key_dict[sessevs_raw_rec[indx]['code']][
                    'response']] += 1
                ev_dict['resp_late'] = True
                if ((key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                        'code']]['response']] == 0) or
                    (key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                        'code']]['response']] > sessevs_raw_rec[indx][
                            'neurotime'])):
                    if ev_dict['first_resp'] == '':
                        ev_dict['first_resp'] = key_dict[sessevs_raw_rec[
                            indx]['code']]['response']
                        ev_dict['first_resp_time'] = sessevs_raw_rec[
                            indx]['neurotime']
                    key_first_time[phase][key_dict[sessevs_raw_rec[indx][
                        'code']]['response']] = sessevs_raw_rec[indx][
                            'neurotime']
                if ((key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                        'code']]['response']] == 0) or
                    (key_last_time[phase][key_dict[sessevs_raw_rec[indx][
                        'code']]['response']] <
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
        if not 'stim_time' in ev_dict:
            continue
        resp_count_all = 0
        for phase in key_count:
            for resp in [key_dict[r]['response'] for r in key_dict]:
                if phase <= 3:
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
            (exp, subj, expnum, sess, expstr, ev_dict['stim_time'],
             ev_dict['phase1_time'], ev_dict['phase2_time'],
             ev_dict['phase3_time'], ev_dict['phase4_time'], 'sr_trial',
             ev_dict['stim_file'], ev_dict['stim_word'],
             ev_dict['stim_version'], ev_dict['stim_mf'],
             ev_dict['stim_mf_version'],  # ev_dict['stim_num'],
             ev_dict['stim_freq'], ev_dict['resp_early'],
             ev_dict['resp_late'], ev_dict['resp_mult'],
             ev_dict['first_resp'], ev_dict['first_resp_time'],
             ev_dict['keys_reversed'], ev_dict['analysis_dir'],
             ev_dict['resp_count_old_1'], ev_dict['resp_count_old_2'],
             ev_dict['resp_count_old_3'], ev_dict['resp_count_old_4'],
             ev_dict['resp_count_new_1'], ev_dict['resp_count_new_2'],
             ev_dict['resp_count_new_3'], ev_dict['resp_count_new_4'],
             ev_dict['resp_count_other_1'], ev_dict['resp_count_other_2'],
             ev_dict['resp_count_other_3'],ev_dict['resp_count_other_4'],
             ev_dict['resp_firsttime_old_1'],
             ev_dict['resp_firsttime_old_2'],
             ev_dict['resp_firsttime_old_3'],
             ev_dict['resp_firsttime_old_4'],
             ev_dict['resp_firsttime_new_1'],
             ev_dict['resp_firsttime_new_2'],
             ev_dict['resp_firsttime_new_3'],
             ev_dict['resp_firsttime_new_4'],
             ev_dict['resp_firsttime_other_1'],
             ev_dict['resp_firsttime_other_2'],
             ev_dict['resp_firsttime_other_3'],
             ev_dict['resp_firsttime_other_4'],
             ev_dict['resp_lasttime_old_1'],
             ev_dict['resp_lasttime_old_2'],
             ev_dict['resp_lasttime_old_3'],
             ev_dict['resp_lasttime_old_4'],
             ev_dict['resp_lasttime_new_1'],
             ev_dict['resp_lasttime_new_2'],
             ev_dict['resp_lasttime_new_3'],
             ev_dict['resp_lasttime_new_4'],
             ev_dict['resp_lasttime_other_1'],
             ev_dict['resp_lasttime_other_2'],
             ev_dict['resp_lasttime_other_3'],
             ev_dict['resp_lasttime_other_4']))
    if indx < len(sessevs_raw_rec):
        raise ValueError('Incomplete trial -- end')
    ev_dtypes = [('exp', 'U8'), ('subject', 'U8'), ('expnum', int),
                 ('session', int), ('expstr', 'U16'), ('stim_time', float),
                 ('phase1_time', float), ('phase2_time', float),
                 ('phase3_time', float), ('phase4_time', float),
                 ('evtype', 'U16'),
                 ('stim_file', 'U32'), ('stim_word', 'U32'),
                 ('stim_version', int), ('stim_mf', 'U1'),
                 ('stim_mf_version', int),  # ('stim_num', int),
                 ('stim_freq', int), ('resp_early', bool), ('resp_late', bool),
                 ('resp_mult', bool), ('first_resp', 'U16'),
                 ('first_resp_time', float), ('keys_reversed', bool),
                 ('analysis_dir', 'U128'), ('resp_count_old_1', int),
                 ('resp_count_old_2', int), ('resp_count_old_3', int),
                 ('resp_count_old_4', int),
                 ('resp_count_new_1', int), ('resp_count_new_2', int),
                 ('resp_count_new_3', int), ('resp_count_new_4', int),
                 ('resp_count_other_1', int), ('resp_count_other_2', int),
                 ('resp_count_other_3', int), ('resp_count_other_4', int),
                 ('resp_firsttime_old_1', float),
                 ('resp_firsttime_old_2', float),
                 ('resp_firsttime_old_3', float),
                 ('resp_firsttime_old_4', float),
                 ('resp_firsttime_new_1', float),
                 ('resp_firsttime_new_2', float),
                 ('resp_firsttime_new_3', float),
                 ('resp_firsttime_new_4', float),
                 ('resp_firsttime_other_1', float),
                 ('resp_firsttime_other_2', float),
                 ('resp_firsttime_other_3', float),
                 ('resp_firsttime_other_4', float),
                 ('resp_lasttime_old_1', float),
                 ('resp_lasttime_old_2', float),
                 ('resp_lasttime_old_3', float),
                 ('resp_lasttime_old_4', float),
                 ('resp_lasttime_new_1', float),
                 ('resp_lasttime_new_2', float),
                 ('resp_lasttime_new_3', float),
                 ('resp_lasttime_new_4', float),
                 ('resp_lasttime_other_1', float),
                 ('resp_lasttime_other_2', float),
                 ('resp_lasttime_other_3', float),
                 ('resp_lasttime_other_4', float)]
    sessevs_rec = np.array(sessevs, ev_dtypes)

    presented_words = {}
    first = np.empty(sessevs_rec['stim_word'].shape, bool)
    delay = np.empty(sessevs_rec['stim_word'].shape, float)
    delay.fill(np.nan)
    for i, word in enumerate(sessevs_rec['stim_word']):
        if i < (len(sessevs_rec['stim_word'])-1):
            delay[i] = sessevs_rec['stim_time'][i+1]-sessevs_rec[
                'phase4_time'][i]
        if word in presented_words:
            first[i] = False
            if 'lag' in presented_words[word]:
                if presented_words[word]['first_indx'] != 0:
                    raise ValueError(
                        'first_indx ' + word +'\n'+str(presented_words[word]))
                #### raise ValueError('Word repeated more than once: '+word)
                print('Word repeated more than once: '+word)
                presented_words[word]['third_indx'] = i
            else:
                presented_words[word]['lag'] = i-presented_words[word][
                    'first_indx']
                presented_words[word]['second_version'] = sessevs_rec[
                    'stim_version'][i]
                presented_words[word]['second_mf'] = sessevs_rec['stim_mf'][i]
                presented_words[word]['second_mf_version'] = sessevs_rec[
                    'stim_mf_version'][i]
        else:
            first[i] = True
            presented_words[word] = {}
            presented_words[word]['first_indx'] = i
            presented_words[word]['first_version'] = sessevs_rec[
                'stim_version'][i]
            presented_words[word]['first_mf'] = sessevs_rec[
                'stim_mf'][i]
            presented_words[word]['first_mf_version'] = sessevs_rec[
                'stim_mf_version'][i]
    lag = np.zeros(sessevs_rec['stim_word'].shape, int)
    paired = np.empty(sessevs_rec['stim_word'].shape, 'U16')
    for i, word in enumerate(sessevs_rec['stim_word']):
        if 'lag' in presented_words[word]:
            lag[i] = presented_words[word]['lag']
            if (presented_words[word]['first_mf'] ==
                presented_words[word]['second_mf']):
                if (presented_words[word]['first_mf_version'] ==
                    presented_words[word]['second_mf_version']):
                    paired[i] = 'identical'
                else:
                    paired[i] = 'same_mf'
            else:
                paired[i] = 'different'
        else:
            lag[i] = -1
            paired[i] = 'unpaired'
        if 'third_indx' in presented_words[word]:
            if presented_words[word]['third_indx'] == i:
                paired[i] = 'third'
    sessevs_rec = append_fields(
        sessevs_rec, 'first', first,
        dtypes=bool, usemask=False, asrecarray=True)
    sessevs_rec = append_fields(
        sessevs_rec, 'delay', delay,
        dtypes=float, usemask=False, fill_value=np.nan, asrecarray=True)
    sessevs_rec = append_fields(
        sessevs_rec, 'lag', lag,
        dtypes=int, usemask=False, fill_value=-999, asrecarray=True)
    sessevs_rec = append_fields(
        sessevs_rec, 'paired', paired,
        dtypes='U16', usemask=False, fill_value='', asrecarray=True)
    if evs is None:
        evs = sessevs_rec
    else:
        evs = np.r_[evs, sessevs_rec]



# #####################################
# #####################################
# #####################################
# ###
# ### Validate against sliceInfo.txt
# ###


# # sliceInfo.txt is coded wrong for this exp?
# # bad_key_info = ['s20e8cr']



# sliceinfo_missing = ['s12e10cr', 's14e17cr', 's14e4cr', 's27e10cr', 's28e13cr']



# slice_infos = sorted(glob(
#     rhino_root +
#     '/scratch/josh/BniData/Subjects/s*/analysis/s*cr/sliceInfo.txt'))
# # slice_infos.extend(sorted(glob(
# #     rhino_root +
# #     '/scratch/josh/BniData/Subjects/s*/analysis/s*cr/crmTrials.txt')))

# crm_trials = sorted(glob(
#     rhino_root +
#     '/scratch/josh/BniData/Subjects/s*/analysis/s*cr/crmTrials.txt'))

# # if we have sliceInfo.txt and crmTrials.txt, remove the
# # crmTrials.txt:
# for crmt  in crm_trials:
#     expstr = crmt.split(os.path.sep)[-2]
#     for si in slice_infos:
#         if expstr in si:
#             crm_trials.remove(crmt)

# slice_infos.extend(crm_trials)


# for sf in slice_infos:
#     slice_info = pd.read_csv(sf, sep='\t')
#     expstr = sf.split(os.path.sep)[-2]
#     if expstr in bad_exps:
#         continue
#     print(sf)    
#     evs_sess = evs[evs['expstr'] == expstr]
#     exp_dur = (evs_sess['phase1_time'][-1]-evs_sess['phase1_time'][0])
#     # print('Experiment', expstr, 'lasted', int(np.round(exp_dur/60e6)), 'minutes.')
#     assert (exp_dur/60e6) > 8, '(exp_dur/60e6) > 8'
#     # slice_info['name'] == evs_sess['stim_file']
#     assert np.alltrue(slice_info['name'] == evs_sess['stim_file']), \
#         "np.alltrue(slice_info['name'] == evs_sess['stim_file'])"
#     assert np.allclose(slice_info['stTime'], evs_sess['phase1_time']), \
#         'np.allclose(slice_info[\'stTime\'], evs_sess[\'phase1_time\'])'
#     assert np.max(evs_sess['stim_time']-evs_sess['phase1_time']) <= 2000, \
#         'np.max(evs_sess[\'stim_time\']-evs_sess[\'phase1_time\']) <= 2000'
#     assert np.allclose(slice_info['endTime'], evs_sess['phase2_time']), \
#         'np.allclose(slice_info[\'endTime\'], evs_sess[\'phase2_time\'])'
#     if expstr in sliceinfo_missing:
#         # if expstr not in bad_times:
#         assert np.allclose(slice_info['firstResp']-slice_info['stTime'],
#                            evs_sess['first_resp_time']-evs_sess['stim_time'],
#                            rtol=0.005, atol=0.0005), \
#                            """slice_info['firstResp']-slice_info['stTime'],
#                            evs_sess['first_resp_time']-evs_sess['stim_time'],
#                            rtol=0.005, atol=0.0005)"""
#     else:
#         tmp1 = slice_info['firstResp']*1000.0
#         tmp2 = evs_sess['first_resp_time']-evs_sess['stim_time']
#         # tmp2 = evs_sess['first_resp_time']-evs_sess['phase1_time']
#         tmp2[~np.isfinite(tmp2)] = -1000
#         assert np.allclose(tmp1[(~slice_info['keyEarly']) &
#                                 (slice_info['firstResp'] > 0)],
#                            tmp2[(~slice_info['keyEarly']) &
#                                 (slice_info['firstResp'] > 0)],
#                            rtol=0.005, atol=0.0005), \
#                            """np.allclose(tmp1[(~slice_info['keyEarly']) &
#                            (slice_info['firstResp'] > 0)],
#                            tmp2[(~slice_info['keyEarly']) &
#                            (slice_info['firstResp'] > 0)],
#                            rtol=0.005, atol=0.0005)"""
#     assert np.alltrue(slice_info['keyEarly'] == evs_sess['resp_early']), \
#         "np.alltrue(slice_info['keyEarly'] == evs_sess['resp_early'])"
#     assert np.alltrue(slice_info['multiPress'] == evs_sess['resp_mult']), \
#         "np.alltrue(slice_info['multiPress'] == evs_sess['resp_mult'])"
#     # if expstr not in bad_key_info:
#     assert np.alltrue(evs_sess['first_resp'][slice_info['oldKey']] == 'old'), \
#         "np.alltrue(evs_sess['first_resp'][slice_info['oldKey']] == 'old')"
#     # else:
#     # assert np.alltrue(evs_sess['first_resp'][slice_info['oldKey']] == 'new'), \
#     #     "np.alltrue(evs_sess['first_resp'][slice_info['oldKey']] == 'new')"
#     assert np.alltrue(evs_sess['first_resp'][
#         slice_info['otherKey']] == 'other'), \
#         """np.alltrue(evs_sess['first_resp'][
#         slice_info['otherKey']] == 'other')"""
#     # if expstr not in bad_times:
#     assert np.abs(np.round(evs_sess['delay']/1000)[
#         :-1] - slice_info['delay'][:-1]).max() <= 2, \
#         """np.abs(np.round(evs_sess['delay']/1000)[
#         :-1] - slice_info['delay'][:-1]).max() <= 2"""
#     assert np.alltrue(evs_sess['paired'][
#         slice_info['isPaired']] != 'unpaired'), \
#         """np.alltrue(evs_sess['paired'][
#         slice_info['isPaired']] != 'unpaired')"""
#     assert np.alltrue(evs_sess['paired'][
#         ~slice_info['isPaired']] == 'unpaired'), \
#         """np.alltrue(evs_sess['paired'][
#         ~slice_info['isPaired']] == 'unpaired')"""
#     assert np.alltrue(evs_sess['paired'][
#         slice_info['pairedWithDup']] == 'same'), \
#         """np.alltrue(evs_sess['paired'][
#         slice_info['pairedWithDup']] == 'same')"""
#     assert np.alltrue(evs_sess['paired'][
#         (~slice_info['pairedWithDup']) &
#         (slice_info['isPaired'])] == 'different'), \
#         """np.alltrue(evs_sess['paired'][
#         (~slice_info['pairedWithDup']) &
#         (slice_info['isPaired'])] == 'different')"""
#     assert np.alltrue(evs_sess['first'] == slice_info['isFirst']), \
#         "np.alltrue(evs_sess['first'] == slice_info['isFirst'])"
#     assert np.alltrue(evs_sess['lag'] == slice_info['lag']), \
#         "np.alltrue(evs_sess['lag'] == slice_info['lag'])"



np.savez_compressed('./data/bni_evs_sr.npz', evs=evs)
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


