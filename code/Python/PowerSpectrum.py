#!/usr/bin/env python

import pickle
import numpy as np
import syncopy as spy
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

# Parameters
saveFig = 1             # 0:not savefigs  1: save figs
freq_range = [15, 90]   # frequency range
fs = 1250               # sampling rate
areas = {
    "V1L": ['VISl'],
    "V1RL": ['VISrl'],
    "Thal-MG": ['MGv'],
    "Thal-LP": ['LP']
}

# Load data
lfpdata = spy.load('allen_lfp.spy')

# Power estimation
L = len(lfpdata.time[0])
base_time = np.arange(0, np.ceil(L/2)-1, 1)
stim_time = np.arange(np.floor(L/2), L-1, 1)
TW_base = [lfpdata.time[0][int(base_time[0])],
           lfpdata.time[0][int(base_time[-1])]]
TW_stim = [lfpdata.time[0][int(stim_time[0])],
           lfpdata.time[0][int(stim_time[-1])]]
data_base = lfpdata.selectdata(latency=TW_base)
data_stim = lfpdata.selectdata(latency=TW_stim)

cfg = []
cfg = spy.get_defaults(spy.freqanalysis)
cfg.method = 'mtmfft'
cfg.output = 'pow'
cfg.pad = 1
cfg.taper = 'hann'
cfg.nTaper = 1
cfg.keeptrials = True
cfg.foilim = freq_range

# Baseline correction
psd_stim = spy.freqanalysis(cfg, data_stim)
psd_base = spy.freqanalysis(cfg, data_base)
psd = []
psd2 = []
psd = psd_stim
powbs = (np.mean(psd_base.data, axis=0))
psd.data = psd_stim.data/powbs
selectch = np.arange(0, len(psd.channel), 1)
pow = psd.data

# Outlier detection
outl = np.nanmean(np.nanmean(pow, axis=2), axis=0) > (np.nanmean(
    np.reshape(pow, (1, -1)), axis=1)+3*np.std(np.reshape(pow, (1, -1)), axis=1))
indices = np.where(outl)[0]
selch = np.delete(selectch, indices)
spy.selectdata(psd, channel=list(lfpdata.channel[selch]))

cfg = spy.StructDict()
psd2 = spy.selectdata(cfg, psd)

# Grand power
channelLBL = lfpdata.info['Label']

GrandPSD = spy.mean(psd2, dim='trials')
iter = 1
Pow = []

# for area_iter in areas.keys():
area_iter = list(areas.keys())[0]
area = []
y = []
y_tmp = []
indexes = []
for i in areas[f"{area_iter}"]:
    indexes.append([index for index in range(
        len(channelLBL)) if channelLBL[index] == i])
chindx = list(map(int, np.concatenate(indexes)))

y_tmp = spy.selectdata(GrandPSD, channel=chindx)
y = spy.mean(y_tmp, dim='channel')
fig, ax = y.singlepanelplot(logscale=False)
ax.set_title(f"ch= {area_iter}")
ax.set_xlabel('Frequency(Hz)')
ax.set_ylabel('Syncopy power')
ax.set_xlim([1, 90])
Pow.append(np.squeeze(y.data))
iter = iter+1

if saveFig == 1:
    plt.savefig(f'power{area_iter}.png')
    plt.savefig(f'power{area_iter}.pdf')

pickle.dump(Pow, open(f"SavedData/PyPow.plk", "wb"))
