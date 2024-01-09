#!/usr/bin/env python

import syncopy as spy
import matplotlib.pyplot as plt

# Settings
plt.rcParams['pdf.fonttype'] = 42
saveFig = 1             # 0:not savefigs  1: save figs

# Load data
spikeData = spy.load('allen_spike.spy')
data_spk = spikeData.trials
samplerate = 1250
offset = -0.25 * samplerate  # Offset in samples for aligning samples to stimulus onset.

# Raster plot of the data.
unitID = '951009344'
unitname = 'unit30'
plt.close('all')
plt.figure()
spikeData.singlepanelplot(unit=unitID)

if saveFig == 1:
    plt.savefig("Raster.png")
    plt.savefig("Raster.pdf")


# Run the PSTH.
cfg = spy.StructDict()
cfg.latency = "minperiod"  # to avoid NaNs
cfg.output = "rate"
cfg.binsize = 0.01  # in seconds
spike_rate = spy.spike_psth(spikeData, cfg)

# Remove channel0 from channel names.
chan_names = [old_name.split('_')[1] for old_name in spike_rate.channel]
spike_rate.channel = chan_names

# Compute the mean.
rate_trl_mean = spy.mean(spike_rate, dim='trials')

# Plot.
fig, ax = rate_trl_mean.singlepanelplot(channel=unitname)
ax.set_ylabel('spike rate (1/s)')
fig.tight_layout()

if saveFig == 1:
    plt.savefig("PSTH.png")
    plt.savefig("PSTH.pdf")
