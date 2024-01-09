#!/usr/bin/env python

import pickle
import numpy as np
import syncopy as spy
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

plt.rcParams['pdf.fonttype'] = 42

# parameters
cohppcGC = 1  # 1: coh  2:ppc 3:GC
saveFig = 1             # 0:not savefigs  1: save figs
freq_range = [1, 100]   # frequency range
fs = 1250               # sampling rate
areas = {
    "V1L": ['VISl'],
    "V1RL": ['VISrl'],
    "Thal-MG": ['MGv'],
    "Thal-LP": ['LP']
}

# Load data
lfpdata = spy.load('allen_lfp.spy')
channelLBL = lfpdata.info['Label']

#  connectivity estimation ..................................
L = len(lfpdata.time[0])
stim_time = np.arange(np.floor(L/2), L-1, 1)
TW_stim = [lfpdata.time[0][int(stim_time[0])],
           lfpdata.time[0][int(stim_time[-1])]]
data_stim_all = lfpdata.selectdata(latency=TW_stim)  # , trials=All


cfg = []
cfg = spy.StructDict()
if cohppcGC == 1:
    cfg.method = 'coh'
elif cohppcGC == 2:
    cfg.method = 'ppc'
elif cohppcGC == 3:
    cfg.method = 'granger'
    np.seterr(divide='ignore', invalid='ignore')

cfg.pad = 1  # 'maxperlen'
# cfg.taper = 'None' #'hann'  # 'dpss' 'hann'
# cfg.nTaper = 1
# cfg.tapsmofrq = 3

area_num = list(areas.keys())
plt.close('all')
plt.figure()
iter = 0

for sender in range(0, 1):  # range(len(area_num)):  # -1
    for receiver in range(1, 2):  # range(len(area_num)):  # sender+1,
        iter = iter+1
        if receiver > sender:
            for sender_ch in areas[area_num[sender]]:
                for receiver_ch in areas[area_num[receiver]]:
                    Con12 = []
                    Con21 = []
                    ch1 = []
                    ch2 = []
                    ppc_val = []
                    ppc_val_FB = []
                    for index, item in enumerate(channelLBL):
                        if item == sender_ch:
                            ch1.append(index)
                        if item == receiver_ch:
                            ch2.append(index)
                    # ch1 = np.nonzero(channelLBL == sender_ch)[0]
                    # ch2 = np.nonzero(channelLBL == receiver_ch)[0]
                    for c1 in ch1:
                        for c2 in ch2:
                            print(c1)
                            chsV1V2 = np.concatenate((c1, c2), axis=None)
                            data_stim_subGroup = spy.selectdata(
                                data_stim_all, channel=chsV1V2)
                            ppc = spy.connectivityanalysis(
                                data_stim_subGroup, cfg)
                            fr = np.where((ppc.freq >= 15) &
                                          (ppc.freq <= 90))[0]
                            fr = fr+1
                            c11 = []
                            c22 = []
                            for index, item in enumerate(chsV1V2):
                                if channelLBL[item] == sender_ch:
                                    c11.append(index)
                                if channelLBL[item] == receiver_ch:
                                    c22.append(index)
                            # ppc_val = spy.selectdata(ppc, frequency=[min(fr), max(fr)], channel_i=c11, channel_j=c22)
                            tmp1 = ppc.data[0]
                            tmp2 = tmp1[fr, c11, c22]
                            Con12.append(tmp2)
                            # # Con12.append(np.squeeze(ppc_val.data))
                            if cohppcGC == 3:
                                ppc_val_FB = spy.selectdata(
                                    ppc, frequency=[min(fr), max(fr)], channel_i=c22, channel_j=c11)
                                # # Con21.append(np.squeeze(ppc_val_FB.data))
                                tmp1 = []
                                tmp2 = []
                                tmp1 = ppc.data[0]
                                tmp2 = tmp1[fr, c22, c11]
                                Con21.append(tmp2)

                    print(sender, receiver)
                    plt.subplot(4, 4, iter)
                    plt.plot(np.squeeze(np.nanmean(Con12, axis=0)), '-',
                             lw=2, c='red', alpha=0.5)
                    if cohppcGC == 3:
                        plt.plot(np.squeeze(np.nanmean(Con21, axis=0)),
                                 '--', lw=2, c='red', alpha=0.5)
                        plt.title(
                            f"__{area_num[sender]}->{area_num[receiver]}")
                    else:
                        plt.title(f"{area_num[sender]},-,{area_num[receiver]}")

                    plt.xticks([], [])
                    if iter == 1:
                        plt.xlabel('Frequency(Hz)')
                        if cohppcGC == 1:
                            plt.ylabel('Coherence')
                        elif cohppcGC == 2:
                            plt.ylabel('PPC')
                        elif cohppcGC == 3:
                            plt.ylabel('GC')
                    if saveFig == 1:
                        plt.savefig(f'coh_ppc_gc_{cohppcGC}.png')
                        plt.savefig(f'coh_ppc_gc_{cohppcGC}.pdf')
Con = {'fr': fr,
       'Con12': Con12,
       'Con21': Con21,
       }

pickle.dump(Con, open(f"SavedData/PyCohPpcGC{cohppcGC}.plk", "wb"))
