#!/usr/bin/env python


# This script is used to download data from the Allen Brain Observatory. It is based on the Allen SDK and
# first downloads the data, then selects a subset of the data, and finally saves the data into Syncopy data files
# and, additionally, in a MATLAB file for usage with Fieldtrip.

import os
import numpy as np
import pandas as pd
from scipy.io import savemat
from tqdm import tqdm, trange
from copy import deepcopy
import syncopy as spy


from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache


# The data in the JSON file determines where downloaded data will be stored. See Allen SDK.
manifest_path = os.path.join(os.path.curdir, 'manifest.json')
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)
sessions = cache.get_session_table()

# Data Parameters
stimType = "flashes"  # 'flashes' 'drifting_gratings'  'gabors'
bipolar = 1  # 1:bipolar derivation   0:unipolar
sample_freq = 1250
sample_freq_rawdata = 30_000

sesnum = 7
ses_id = sessions.index.values[sesnum]
session = cache.get_session_data(ses_id)
probes = cache.get_probes()


# Cut trials
print(pd.unique(session.stimulus_presentations.stimulus_name))
presentation_table = session.stimulus_presentations[session.stimulus_presentations.stimulus_name == stimType]
presentation_times = presentation_table.start_time.values
presentation_ids = presentation_table.index.values
win = np.arange(-.25, .25, 1 / sample_freq)
stim_lbl = session.get_stimulus_table([stimType])['color']





prob_range = session.probes.index.values

"""
Read Allen data and save the spikes in Syncopy file format.
"""
def read_spy_spikes():

    print("Creating Syncopy spike data file...")

    offset = -0.25 * sample_freq  # in samples

    # first assemble units of interest
    units_oi = []
    for prob_ind in tqdm(prob_range):
        units_oi += list(session.units[(session.units.probe_id == prob_ind) &
                                       (session.units.firing_rate > 10) &
                                       (session.units.nn_hit_rate > 0.95)].index)

    # now loop over trials defined by window around presentation_times

    spy_spikes = []
    trldef_rows = []
    active_units = []
    for stim_time in tqdm(presentation_times):

        start = stim_time - 0.25
        stop = stim_time + 0.25

        # 2d array with columns: samples, channel, unit.
        # Channel is not given/interesting for spike data and has to be 0.
        channel_id = 0

        spy_trl = []

        for unit_idx, unit_id in enumerate(units_oi):

            unit_spike_times = session.spike_times[unit_id]
            # index time window
            bool_index = (unit_spike_times > start) & (unit_spike_times < stop)
            # relevant spikes in samples
            unit_spike_samples = (
                unit_spike_times[bool_index] + start) * sample_freq

            # Skip empty arrays, totally (for all trials) inactive
            # units have to be skipped entirely.
            if unit_spike_samples.size == 0:
                continue

            active_units.append(unit_id)

            spy_trl.append(np.array([(sp_sample, channel_id, unit_idx) for
                                     sp_sample in unit_spike_samples], dtype=np.int64))

        spy_trl = np.vstack(spy_trl)
        # sort the samples
        spy_trl = spy_trl[spy_trl[:, 0].argsort()]
        spy_spikes.append(spy_trl)

        trldef_rows.append([spy_trl[0, 0], spy_trl[-1, 0], offset])

    trldef = np.vstack(trldef_rows)
    spy_spikes = spy.SpikeData(data=spy_spikes, samplerate=sample_freq)

    spy_spikes.trialdefinition = trldef
    # unit labels have to be strings
    spy_spikes.unit = list(map(str, set(active_units)))

    spy_spikes.info['fsample_raw'] = sample_freq_rawdata   # Add custom metadata: raw sampling rate.

    # Finally save to disc.
    spy_spikes.save(os.path.join(os.path.curdir, 'allen_spike'), overwrite=True)
    # return Spikes



"""
Read Allen data and save the LFP data in Syncopy file format.
"""
def read_spy_lfp():

    print("Creating Syncopy LFP data file...")

    loc = []
    ds = []
    twin = []
    lfps = []
    lfp_channel_id = []

    
    for prob_ind in tqdm(prob_range):
        lfps.append(session.get_lfp(prob_ind))

    for tr in trange(len(presentation_times)):
        lfp_buf = []
  
        # Append LFP of all channels
        for _, lfp in enumerate(lfps):
            ch_range = lfp.channel.values
            for z in range(len(ch_range)):
                channel_ids = ch_range[z]
                try:
                    if bipolar == 1:
                        lfp_tmp = (lfp.sel(time=(presentation_times[tr]+win), method='nearest', channel=ch_range[z]))-(
                            lfp.sel(time=(presentation_times[tr]+win), method='nearest', channel=ch_range[z+1]))
                        loc_tmp = (session.channels.ecephys_structure_acronym[channel_ids])
                    elif bipolar == 0:
                        lfp_tmp = (lfp.sel(time=(presentation_times[tr]+win), method='nearest', channel=channel_ids))
                        loc_tmp = (session.channels.ecephys_structure_acronym[channel_ids])

                    if tr == 0:
                        loc.append(loc_tmp)
                        if bipolar == 1:
                            lfp_channel_id.append([ch_range[z], ch_range[z+1]])
                        elif bipolar == 0:
                            lfp_channel_id.append(channel_ids)
                    lfp_buf.append(lfp_tmp)
                    lfp_buf[-1] = np.nan_to_num(lfp_buf[-1].values.reshape(1, -1))
                except:
                    pass

        ds.append(np.concatenate(lfp_buf, axis=0).T)
        twin.append(deepcopy(win))

    # Put data into Syncopy objects.
    lfp = spy.AnalogData(data=ds, samplerate=sample_freq)
    lfp.info['Label'] = loc
    lfp.info['stim_lbl'] = list(stim_lbl)
    lfp.info['Twin'] = [arr.tolist() for arr in twin]



"""
Read Allen data and save it for loading into Fieldtrip, in Matlab MAT file format.
"""
def read_ft():

    print("Creating Fieldtrip data file...")

    spikes = []
    loc = []
    locspk = []
    ds = []
    ds_spk = []
    twin = []
    lfps = []
    spk_channel_id = []
    lfp_channel_id = []

    fieldtrip = {"Label": [],  # cell(1*chNum)
                 "Labelspk": [],
                 "chIDlfp": [],
                 "chIDspk": [],
                 "trialinfo": {},  # double(trialNum*info)
                 "trial": [],          # cell(1*trialNum) >> (chNum*sample)
                 "trial_spike": [],
                 "time": [],           # cell(1*trialNum) >> (1*sample)
                 "fsample": sample_freq,        # sampling rate  1250
                 "fsample_raw": sample_freq_rawdata,
                 # 1*1 structure (other infor) >> not necessary
                 "cfg": []
                 }

    # Append spikes of all channels
    for prob_ind in tqdm(prob_range):
        lfps.append(session.get_lfp(prob_ind))
        units_of_interest = (session.units[(session.units.probe_id == prob_ind) &
                                           (session.units.firing_rate > 10) &
                                           (session.units.nn_hit_rate > 0.95)])
        for i in range(len(units_of_interest)):
            unit_id = units_of_interest.index.values[i]
            channel_index = units_of_interest .loc[unit_id].probe_channel_number
            spk_channel_id.append(session.channels[(session.channels.probe_channel_number == channel_index) &
                                                   (session.channels.probe_id == prob_ind)].index.values[0])
            spikes.append(session.spike_times[unit_id])

    for tr in trange(len(presentation_times)):
        lfp_buf = []
        spk_buf = []
        for unitNum in range(len(spikes)):
            tmp = np.zeros(len(win))
            sp = (spikes[unitNum][(spikes[unitNum] > presentation_times[tr]+win[0]) &
                                  (spikes[unitNum] < presentation_times[tr]+win[-1])])-(presentation_times[tr]+win[0])
            for i in sp:
                tmp[int(i*fieldtrip['fsample'])] = 1
            spk_buf.append(tmp)
            if tr == 0:
                locspk.append(
                    session.channels.ecephys_structure_acronym[spk_channel_id[unitNum]])

        # Append LFP of all channels
        for _, lfp in enumerate(lfps):
            ch_range = lfp.channel.values
            for z in range(len(ch_range)):
                channel_ids = ch_range[z]
                try:
                    if bipolar == 1:
                        lfp_tmp = (lfp.sel(time=(presentation_times[tr]+win), method='nearest', channel=ch_range[z]))-(
                            lfp.sel(time=(presentation_times[tr]+win), method='nearest', channel=ch_range[z+1]))
                        loc_tmp = (
                            session.channels.ecephys_structure_acronym[channel_ids])
                    elif bipolar == 0:
                        lfp_tmp = (lfp.sel(
                            time=(presentation_times[tr]+win), method='nearest', channel=channel_ids))
                        loc_tmp = (
                            session.channels.ecephys_structure_acronym[channel_ids])
                    if tr == 0:
                        loc.append(loc_tmp)
                        if bipolar == 1:
                            lfp_channel_id.append([ch_range[z], ch_range[z+1]])
                        elif bipolar == 0:
                            lfp_channel_id.append(channel_ids)
                    lfp_buf.append(lfp_tmp)
                    lfp_buf[-1] = np.nan_to_num(
                        lfp_buf[-1].values.reshape(1, -1))
                except:
                    pass

        ds.append(np.concatenate(lfp_buf, axis=0).T)
        print(ds[-1].shape)
        twin.append(deepcopy(win))
        ds_spk.append(list(map(list, zip(*spk_buf))))

    fieldtrip["trial"] = ds
    fieldtrip["spike"] = ds_spk
    fieldtrip["Label"] = loc
    fieldtrip["Labelspk"] = locspk
    fieldtrip["time"] = twin
    fieldtrip["trialinfo"] = {'stim_lbl': stim_lbl,
                              'chIDlfp': lfp_channel_id,
                              'chIDspk': spk_channel_id
                              }
    fieldtrip2 = deepcopy(fieldtrip)
    fieldtrip2["trialinfo"] = []
    fieldtrip2["stim_lbl"] = list(stim_lbl)

    savemat(os.path.join(os.path.curdir, 'allen_FT.mat'), fieldtrip2)


if __name__ == "__main__":
    read_spy_spikes()
    read_spy_lfp()
    read_ft()
