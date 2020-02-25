import sys
import numpy as np

import obspy
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import pk_baer

import pprint

import os

from obspy import read
import matplotlib.pyplot as plt

st = obspy.Stream()
cwd = os.getcwd()

#Zips together data from file
files = os.listdir(cwd)
#print(files)
#print("length is " + str(len(files)))

for filename in os.listdir(cwd):
    if filename.endswith('.htm'):
        st += read(filename)


#Getting stream using select

st2 = st.select(station="MBGE", component = 'Z')
#print(st2.__str__(extended=True))

#Selecting event time and slicing
dt = obspy.UTCDateTime("1997-02-11T00:05:00")
ft = obspy.UTCDateTime("1997-02-11T23:45:00")
st2 = st2.slice(starttime = dt, endtime = ft)

#Trace of sliced data
tr2 = st2[0]
df = tr2.stats.sampling_rate

#Spectrogram Plotter

#st2.plot()
#st2.spectrogram(log=False, title='MBGE STATION' + str(st[0].stats.starttime))

#Making a copy of the original data
tr_filt = tr2.copy()

#bandpass filter
tr_filt.filter('bandpass', freqmin=2, freqmax=6, corners=2, zerophase=False)

print(tr_filt)

"""
# Plotting the raw and filtered data
t = np.arange(0, tr2.stats.npts / tr2.stats.sampling_rate, tr2.stats.delta)
plt.plot(t, tr_filt.data, 'k')
plt.ylabel('Bandpassed Data (2HZ --> 6HZ)')
plt.xlabel('Time')
plt.suptitle(tr2.stats.starttime)
plt.show()
"""
#Getting the trigger
trig = classic_sta_lta(tr_filt.data, int(5 * df), int(10 * df))

#plot_trigger(tr_filt, trig, 1.25, 0.75)

#Getting a list of trigger start times
time_list_on = (obspy.signal.trigger.trigger_onset(trig, 1.25, 0.75, max_len=9e+99)/df)[:,0]
time_list_off = (obspy.signal.trigger.trigger_onset(trig, 1.25, 0.75, max_len=9e+99)/df)[:,1]

print(time_list_on)
print(time_list_off)

#Looping through Sliced Stream Objects

for i in range(len(time_list_on)):
    st = obspy.Stream()
    st[i] = st.slice(starttime = time_list_on[i], endtime = time_list_off[i])
    st[i].plot()
