import sys
import numpy as np
import obspy
from obspy import read
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import z_detect
from obspy.signal.trigger import carl_sta_trig
from obspy.signal.trigger import recursive_sta_lta
import matplotlib.pyplot as plt



#Read the seisomic data
st = read("EC.RETU..SHZ.D.2013.mseed")

st = st.select(component="Z")

#Trace
tr = st[0]

df = tr.stats.sampling_rate


print("Trace is: " + str(st[0]))
print(st[0].stats)
print(st[0].data)

#Day Plot of this thing
st.plot(type="dayplot", interval=60, right_vertical_labels=False, vertical_scaling_range=5e3, one_tick_per_line=True, color=['k', 'r', 'b', 'g'], show_y_UTC_label=False, events={'min_magnitude': 7})

#Making a copy of the original data
tr_filt = tr.copy()

#bandpass filter
tr_filt.filter('bandpass', freqmin=0.2, freqmax=0.21, corners=2, zerophase=False)

# Plotting the raw and filtered data
t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
plt.plot(t, tr_filt.data, 'k')
plt.ylabel('Bandpassed Data (2HZ --> 6HZ)')
plt.xlabel('Time')
plt.suptitle(tr.stats.starttime)
plt.show()


#Trigger Picker

cft = classic_sta_lta(tr_filt.data, int(5 * df), int(10 * df))

plot_trigger(tr_filt, cft, 1.5, 1)
