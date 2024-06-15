import numpy as np
import re
import os
from scipy import signal
import h5py
from glob import glob


class Ncs:
    """
    Read and pre-process NCS data.
    
    Parameters
    ----------
    analysis_dir: str
        Path to directory containing NCS files.
    header_size: int, optional
        Size of the NCS file headers in bytes.
    header_size: int, optional
        Size of the NCS file blocks in samples.
    """
    def __init__(self, analysis_dir, header_size=16*1024, block_size=512):
        self.analysis_dir = analysis_dir
        self.header_size = header_size
        self.block_size = block_size

        self.ncs_channels = sorted(glob(os.path.join(
            self.analysis_dir, 'CSC?.Ncs')))
        self.ncs_channels.extend(sorted(glob(
            os.path.join(self.analysis_dir, 'CSC??.Ncs'))))
        self.ncs_channels.extend(sorted(glob(
            os.path.join(self.analysis_dir, 'CSC???.Ncs'))))
        assert len(glob(os.path.join(
            self.analysis_dir, 'CSC*.Ncs'))) == len(self.ncs_channels), \
            """len(glob(os.path.join(
            self.analysis_dir, 'CSC*.Ncs'))) == len(self.ncs_channels)"""

        self.headers = [self._stat_ncs(chan) for chan in self.ncs_channels]
        # these keys vary across channels:
        variable_keys = ['ADChannel', 'Channel']
        # If 'ADBitVolts' differs across channels, scaling for
        # converting to microvolts needs to be adjusted in the
        # resample_scale method below. It seems to be constant across
        # channels.
        for key in self.headers[0]:
            if key in variable_keys:
                continue
            assert np.alltrue([h[key] == self.headers[0][key]
                               for h in self.headers]), 'key match error: '+key
        # threshold used for determining the time range across channels:
        self.timediff_thresh = 1000/self.headers[0]['SamplingFrequency']*1e6
        # theshold for differences in aligned timestamps:
        self.timestampdiff_thresh = 1/(
            self.headers[0]['SamplingFrequency']/1e6)/1.9
        # threshold for differences in aligned timestamps across blocks:
        self.timestampdiff_block = (1/(
            self.headers[0]['SamplingFrequency']/1e6))*self.block_size
        self.timestampdiff_block_tolerance = 5000  # 5 ms
        self.timestampdiff_block_allowed = (self.timestampdiff_block +
                                            self.timestampdiff_block_tolerance)
        # minimal duration for a good chunck of data:
        self.min_good_dur = 30*1e6  # 30 s in micro-seconds
        # data type for NCS files:
        self.ncs_dtype = [('timestamp', 'uint64'), ('channel', 'uint32'),
                          ('sample_rate', 'uint32'), ('nb_valid', 'uint32'),
                          ('samples', 'int16', (self.block_size,))]        
        
    def read_ncs_data(self, verbose=False):
        """
        Read in the ncs data.
        
        Parameters
        ----------
        verbose: boolean, optional
            Verbose output
        """
        if verbose:
            print('Calculating time range')
        self._get_timerange()
        if verbose:
            print('Read data')
        self._get_dat(verbose)
        self.data_times = []
        for di in self.data_indices:
            self.data_times.append(
                [(self.times[indx[0]],
                  self.times[indx[1]]) for indx in di])
    def resample_scale(self, freq, scale=True, verbose=False):
        """
        Resample data and scale to microvolts.
        
        Parameters
        ----------
        freq: float
            New sampling frequency
        scale: boolean, optional
            Whether to scale to microvolts
        verbose: boolean, optional
            Verbose output
        """
        # Resampling
        # divide by 1e6 to convert to sec.:
        num = int(np.round((self.minmax_time - self.maxmin_time) / 1e6 * freq))
        # # doing it all at once results in memory errors:
        # if verbose:
        #     print('Resampling!')
        # self.data = signal.resample(self.dat, num = num, axis=0)
        dat_resampled = np.empty((len(self.data), num))
        dat_resampled.fill(np.nan)
        for i in range(len(self.data)):
            if verbose:
                print('Resampling', i+1, '/', len(self.data))
            dat_resampled[i] = signal.resample(
                self.data[i], num=num)
        assert np.sum(~np.isfinite(dat_resampled)) == 0,\
            'np.sum(~np.isfinite(dat_resampled)) == 0'
        self.times = (
            np.arange(0, num) * (
                self.times[1] - self.times[0]) *
            self.data.shape[1] / num + self.times[0])
        self.data = dat_resampled
        for i in range(len(self.headers)):
            self.headers[i]['SamplingFrequency'] = freq
        # update indices for "good" data:
        self.data_indices = []
        for dt in self.data_times:
            self.data_indices.append(
                [(np.min(np.where(self.times >= tme[0])[0]),
                  np.max(np.where(self.times <= tme[1])[0]))
                 for tme in dt])
        if scale:
            # convert to microvolts
            assert np.alltrue([h['ADBitVolts'] == self.headers[
                0]['ADBitVolts'] for h in self.headers]), \
                'ADBitVolts check'
            self.data *= (self.headers[0]['ADBitVolts'] * 1e6)

    def to_hdf(self, datfolder='./data/ncs/', verbose=True):
        """
        Save data to HDF file.
        
        Parameters
        ----------
        datfolder: str, optional
            Path to folder to store the HDF file in
        verbose: boolean, optional
            Verbose output
        """
        assert np.alltrue([h['ExpStr'] == self.headers[0]['ExpStr']
                           for h in self.headers]), 'expstr test'
        filename = datfolder+self.headers[0]['ExpStr']+'_ncs.hdf'
        if verbose:
            print('Creating', filename)
        try:
            f = h5py.File(filename, 'w-')
        except IOError as e:
            print(e)
            return
        ncs_data = f.create_dataset('ncs/data', data=self.data,
                                    compression="gzip", compression_opts=9)
        ncs_times = f.create_dataset(
            'ncs/times', data=self.times, compression="gzip",
            compression_opts=9)
        ncs_channels = f.create_dataset(
            'ncs/channels',
            data=np.array([h['Channel'].encode() for h in self.headers]),
            compression="gzip", compression_opts=9)
        ncs_adchannels = f.create_dataset(
            'ncs/ad_channels',
            data=np.array([int(h['ADChannel']) for h in self.headers]),
            compression="gzip", compression_opts=9)
        ncs_times.make_scale('times')
        ncs_channels.make_scale('channels')
        ncs_adchannels.make_scale('ad_channels')
        ncs_data.dims[0].attach_scale(ncs_channels)
        ncs_data.dims[0].attach_scale(ncs_adchannels)
        ncs_data.dims[1].attach_scale(ncs_times)
        for key in self.headers[0]:
            if (key == 'Channel') or (key == 'ADChannel'):
                continue
            for h in self.headers[1:]:
                assert h[key] == self.headers[0][key],\
                    'h[key] == self.headers[0][key]'
            if isinstance(self.headers[0][key], str):
                ncs_data.attrs[key] = self.headers[0][key].encode()
            else:
                ncs_data.attrs[key] = self.headers[0][key]
        # to encode as HDF data_indices and data_times need to be
        # castable into rectangular arrays (rather than lists of
        # potentially different lengths) so we're filling potentially
        # empty spaces with (0, 0):
        maxlen = 0
        for di in self.data_indices:
            if len(di) > maxlen:
                maxlen = len(di)
        for i in range(len(self.data_indices)):
            for j in range(maxlen-len(self.data_indices[i])):
                self.data_indices[i].append((0, 0))
                self.data_times[i].append((0, 0))
        ncs_data.attrs['data_indices'] = self.data_indices
        ncs_data.attrs['data_times'] = self.data_times
        f.close()

    def _stat_ncs(self, channel_file):
        """
        Generate dictionary with the file's header parameters.
        
        Parameters
        ----------
        channel_file: str
            Path to an ncs neuralynx file
        """
        header_keys = [('NLX_Base_Class_Type', None),
                       ('AmpHiCut', float),
                       ('ADChannel', None),
                       ('ADGain', float),
                       ('AmpGain', float),
                       ('SubSamplingInterleave', int),
                       ('ADMaxValue', int),
                       ('ADBitVolts', float),
                       ('SamplingFrequency', float),
                       ('AmpLowCut', float),
                       ('HardwareSubSystemName', None),
                       ('HardwareSubSystemType', None)]
        
        # load header from file
        with open(channel_file, 'rb') as f:
            txt_header = f.read(self.header_size)
        txt_header = txt_header.strip(b'\x00').decode('latin-1')
        
        # find values and make dict
        info = {}
        for k, type_ in header_keys:
            pattern = '-(?P<name>' + k + ')\t(?P<value>[\S ]*)'
            matches = re.findall(pattern, txt_header)
            for match in matches:
                name = match[0]
                val = match[1].rstrip(' ')
                if type_ is not None:
                    val = type_(val)
                info[name] = val
        info['Subject'] = channel_file.split(os.path.sep)[-4]
        info['ExpStr'] = channel_file.split(os.path.sep)[-2]
        info['Channel'] = channel_file.split(os.path.sep)[-1].split('.')[0]
        return info
        
    def _get_timerange(self):
        """Calculate time range that's common across channels"""
        self.maxmin_time = None  # -np.inf
        self.minmax_time = None  # np.inf
        for channel_file in self.ncs_channels:
            data = np.memmap(channel_file, dtype=self.ncs_dtype, mode='r',
                             offset=self.header_size)
            if self.maxmin_time is None:
                self.maxmin_time = data['timestamp'][0]
                # the last timestamp corresponds to the first sample in
                # the last block. We have an additional
                # (data['nb_valid'][-1]-1) samples in that block and add
                # the corresponding time increments to the self.minmax_time:
                self.minmax_time = data['timestamp'][-1]+(
                    data['nb_valid'][-1]-1)*(
                        1/self.headers[0]['SamplingFrequency']*1e6)
                continue
            if self.maxmin_time < data['timestamp'][0]:
                timediff = data['timestamp'][0]-self.maxmin_time
                # only adjust if time difference is below threshold
                # (otherwise difference is likely to be due to corrupt
                # samples):
                if timediff < self.timediff_thresh:
                    self.maxmin_time = data['timestamp'][0]
            # if first time-stamp is substantially below previous
            # ones, previous ones could be corrupt, so we're resetting
            # self.maxmin_time to the earlier start time:
            if data['timestamp'][0] < (self.maxmin_time-self.timediff_thresh):
                self.maxmin_time = data['timestamp'][0]
            if self.minmax_time > data['timestamp'][-1]+(
                    data['nb_valid'][-1]-1)*(
                        1/self.headers[0]['SamplingFrequency']*1e6):
                timediff = self.minmax_time - data['timestamp'][-1]+(
                    data['nb_valid'][-1]-1)*(
                        1/self.headers[0]['SamplingFrequency']*1e6)
                # only adjust if time difference is below threshold
                # (otherwise difference is likely to be due to corrupt
                # samples):
                if timediff < self.timediff_thresh:
                    self.minmax_time = data['timestamp'][-1]+(
                        data['nb_valid'][-1]-1)*(
                            1/self.headers[0]['SamplingFrequency']*1e6)
            # if last time-stamp is substantially larger than previous
            # ones, previous ones could be corrupt, so we're resetting
            # self.minmax_time to the later end time:
            if data['timestamp'][-1] > (self.minmax_time+self.timediff_thresh):
                self.minmax_time = data['timestamp'][-1]+(
                    data['nb_valid'][-1]-1)*(
                        1/self.headers[0]['SamplingFrequency']*1e6)

    def _get_dat(self, verbose):
        """Load data"""
        dur = self.minmax_time-self.maxmin_time
        assert dur/60e6 > 6, 'dur/60e > 6'  # 6 min. minimum
        self.data = np.zeros((len(self.ncs_channels), int(np.round(
            (dur/1e6)*self.headers[0]['SamplingFrequency']))), dtype=np.int16)
        self.times, self.step = np.linspace(
            self.maxmin_time, self.minmax_time, num=self.data.shape[1],
            retstep=True)
        assert np.allclose(self.step, 1/(
            self.headers[0]['SamplingFrequency']/1e6)),\
            'np.allclose(step, 1/(self.headers[0][\'SamplingFrequency\']/1e6))'
        # the start and end indices in self.data for good periods of data:
        self.data_indices = []
        for c, channel_file in enumerate(self.ncs_channels):
            self.data_indices.append([])
            if verbose:
                print(c+1, '/', len(self.ncs_channels), channel_file)
            data = np.memmap(channel_file, dtype=self.ncs_dtype, mode='r',
                             offset=self.header_size)
            block_indx = 0
            while block_indx < len(data):
                start_indx = block_indx
                last_time = data[start_indx]['timestamp']
                # find a streak of blocks with time stamps within the
                # tolerance (i.e. without interruption):
                for b, this_block in enumerate(data[start_indx:]):
                    if ((this_block['timestamp'] - last_time) >
                        self.timestampdiff_block_allowed):
                        break
                    last_time = this_block['timestamp']
                block_indx += b
                dur = data[block_indx]['timestamp'] - data[
                    start_indx]['timestamp']
                # dismiss good streak if too short:
                if dur < self.min_good_dur:
                    block_indx += 1
                    continue
                good_samples = np.sum([d['nb_valid']
                                       for d in data[start_indx:block_indx]])
                # the timestamps for the good data (one for each
                # sample, rather than only one per block as in the NCS
                # file):
                good_ts = np.linspace(
                    data[start_indx]['timestamp'],
                    data[start_indx]['timestamp'] + self.step * good_samples,
                    num=good_samples, endpoint=True)
                assert np.allclose(
                    good_ts[::self.block_size],
                    [d['timestamp']
                     for d in data[start_indx:block_indx]],
                    rtol=3e-05), 'ts comp'
                # align indices for good streak to the full good data
                # array:
                goodts_indices = {'start': 0, 'end': len(good_ts)-1,
                                  'align': int(len(good_ts)/2)}
                goodts_absdiff = np.abs(self.times-good_ts[
                    goodts_indices['align']])
                goodts_minabsdiff = goodts_absdiff.min()
                assert goodts_minabsdiff < self.timestampdiff_thresh,\
                    'goodts_minabsdiff'
                datts_indices = {
                    'align': np.where(
                        goodts_absdiff == goodts_minabsdiff)[0][0]}
                datts_indices['start'] = datts_indices['align']-int(
                    len(good_ts)/2)
                datts_indices['end'] = datts_indices['start']+len(good_ts) - 1
                # account for the possibility that the good streak
                # extends beyond the overal boundaries for good data
                # and truncate accordingly:
                valid_samples_adjust = 0
                valid_samples_offset = 0
                if datts_indices['start'] < 0:
                    goodts_indices['start'] -= datts_indices['start']
                    valid_samples_offset -= datts_indices['start']
                    datts_indices['start'] = 0
                if datts_indices['end'] >= len(self.times):
                    goodts_indices['end'] -= (datts_indices['end'] -
                                              len(self.times)) + 1
                    valid_samples_adjust += (datts_indices['end'] -
                                             len(self.times) + 1)
                    datts_indices['end'] = len(self.times) - 1
                assert goodts_indices['start'] >= 0,\
                    "goodts_indices['start'] >= 0"
                assert goodts_indices['start'] < len(good_ts),\
                    "goodts_indices['start'] < len(good_ts)"
                assert goodts_indices['end'] >= 0, "goodts_indices['end'] >= 0"
                assert goodts_indices['end'] < len(good_ts),\
                    "goodts_indices['end'] < len(good_ts)"
                assert goodts_indices['end'] > goodts_indices['start'],\
                    "goodts_indices['end'] > goodts_indices['start']"
                assert datts_indices['start'] >= 0, "datts_indices['start']"
                assert datts_indices['start'] < len(self.times),\
                    "datts_indices['start'] < len(self.times)"
                assert datts_indices['end'] > 0, "datts_indices['end']"
                assert datts_indices['end'] < len(self.times),\
                    "datts_indices['end'] < len(self.times)"
                assert datts_indices['end'] > datts_indices['start'],\
                    "datts_indices['end'] > datts_indices['start']"
                # # save indices of good streak:
                # self.data_indices[-1].append((datts_indices['start'],
                #                               datts_indices['end']))
                
                # with alignment we erred on side of lower index (when
                # len(good_ts)/2 was not an int). Check if moving up
                # one index improves things (could limit to cases with
                # uneven len(good_ts), but might as well always
                # check):
                ts_diff = self.times[
                    datts_indices['start']:(datts_indices['end']+1)]-good_ts[
                        goodts_indices['start']:(goodts_indices['end']+1)]
                if (datts_indices['end'] < (len(self.times) - 1) and
                    goodts_indices['end'] < (len(good_ts) - 1)):
                    ts_diff_test = (
                        self.times[(datts_indices['start']+1):
                                   (datts_indices['end']+1)] -
                        good_ts[goodts_indices['start']:
                                (goodts_indices['end']+1)])
                    if ts_diff_test.max() < ts_diff.max():
                        assert ts_diff_test.mean() < ts_diff.mean(),\
                            'ts_diff_test.max() < ts_diff.max()'
                        datts_indices['start'] += 1
                        datts_indices['end'] += 1
                        ts_diff = ts_diff_test
                assert ts_diff.mean() < self.timestampdiff_thresh,\
                    'ts_diff.mean() < self.timestampdiff_thresh'
                assert ts_diff.max() < (self.timestampdiff_block/2),\
                    'ts_diff.max() < (self.timestampdiff_block/2) '+str(
                        ts_diff.max())+' '+str(self.timestampdiff_block)
                valid_samples = np.array([d['nb_valid'] for d in data[
                    start_indx:block_indx]])
                assert np.alltrue(valid_samples[:-1] == self.block_size),\
                    'np.alltrue(valid_samples[:-1] == self.block_size)'
                assert valid_samples[-1] <= self.block_size,\
                    'valid_samples[-1] <= self.block_size'
                self.data[c][datts_indices['start']:
                             (datts_indices['end']+1)] = np.ravel(
                                [d['samples'] for d in data[
                                    start_indx:block_indx]])[
                                        int(valid_samples_offset):
                                        int(np.sum(valid_samples) -
                                            valid_samples_adjust)]
                # save indices of good streak:
                self.data_indices[-1].append((datts_indices['start'],
                                              datts_indices['end']))
                block_indx += 1
