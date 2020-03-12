import sys
import numpy as np
import math

import obspy
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import pk_baer
from obspy import UTCDateTime

from scipy.fftpack import next_fast_len

import pprint

import os

from obspy import read

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

#
# converting seconds into UTC, hopefully
#
def convert(seconds):
    seconds = seconds % (24 * 3600)
    hour = seconds/3600
    seconds %= 3600
    minutes = seconds/60
    seconds %= 60

    return "%d:%02d:%02d" % (hour, minutes, seconds)


st = obspy.Stream()
cwd = os.getcwd()

#
#Zips together data from file
#
files = os.listdir(cwd)
#print(files)
#print("length is " + str(len(files)))

for filename in os.listdir(cwd):
    if filename.endswith('.htm'):
        st += read(filename)
print(st.__str__(extended=True))
#
#Getting stream using select
#
st2 = st.select(station = "MBGE", component = 'Z')
#print(st2.__str__(extended=True))

#
#Selecting start time and slicing, so we only look at one days worth of data
#
dt = obspy.UTCDateTime("1997-02-12T00:00:00")

print("dt is: " + str(dt))

ft  = obspy.UTCDateTime("1997-02-13T00:00:00")

print("ft is: " + str(ft))

dt_in_seconds = dt.second

st2 = st2.slice(starttime = dt, endtime = ft)
print(st2)
st2.merge(fill_value = 'interpolate')
st2.split()

#
#Dayplot
#
#st2.plot(type="dayplot", interval=60, title = "Day 1 Plot for MBGA Station", right_vertical_labels=False, vertical_scaling_range=5e3, one_tick_per_line=True, color=['k', 'r', 'b', 'g'], show_y_UTC_label=False, events={'min_magnitude': 6.5})


total_blip_list = [] #Total number of events at each point
blip_list = [] #Number of events in a given window of time
i_list = [] #number of times we loop


second_count = 0 #counter for seconds
number_of_blips = 0

print(st2)

tr_filt = st2[0]
df = tr_filt.stats.sampling_rate

#Spectrogram Plotter
st2.plot()
#st2.spectrogram(log=False, title='MBGE STATION' + str(st[0].stats.starttime))

#bandpass filter
tr_filt.filter('bandpass', freqmin=2, freqmax=6, corners=2, zerophase=True)

"""
# Plotting the raw and filtered data
t = np.arange(0, tr2.stats.npts / tr2.stats.sampling_rate, tr2.stats.delta)
plt.plot(t, tr_filt.data, 'k')
plt.ylabel('Bandpassed Data (2HZ --> 6HZ)')
plt.xlabel('Time')
plt.suptitle(tr2.stats.starttime)
plt.show()
"""

#
#Getting the trigger
#
trig = classic_sta_lta(tr_filt.data, int(5 * df), int(10 * df))
plot_trigger(tr_filt, trig, 1.65, 0.65)

n_picks = len(obspy.signal.trigger.trigger_onset(trig, 1.65, 0.65, max_len=9e+99))

max_amp_list = np.zeros(n_picks)

rms_list = np.zeros(n_picks)

final_time_list = np.zeros(n_picks)

var_list = np.zeros(n_picks)

print(n_picks)
print(trig)



for pick in range(n_picks):

    start = (obspy.signal.trigger.trigger_onset(trig, 1.65, 0.65, max_len=9e+99)[pick][0])/df

    start = dt + start

    event = st2.slice(starttime = start - 2, endtime = start + 10)
    event = event.split()
    event = event.detrend('demean')

    #final_time_list[pick] = mdates.date2num(start.datetime)

    max_amp_list[pick] = np.absolute(event[0]).max()

    #rms_list[pick] = np.sqrt(np.mean(np.square(event[0])))

    #var_list[pick] = np.var(event[0])

#RMS Values Plotting
"""plt.plot_date(final_time_list, rms_list, c ='crimson', xdate = True)
plt.title("Day 2 10:30 RMS/Time")
plt.xlabel("Time")
plt.ylabel("RMS")
plt.show()

#Variance Plotting
plt.plot_date(final_time_list, var_list, c ='crimson', xdate = True)
plt.title("Day 2 10:30 Var/Time")
plt.xlabel("Time")
plt.ylabel("Var")
plt.show()

#RMS vs Max Amp
plt.scatter(rms_list, max_amp_list, c = 'crimson')
plt.title("Day 2 10:30 RMS/AMP")
plt.xlabel("RMS")
plt.ylabel("Maximum amplitude")
plt.show()
"""
#np.save("Day 2 10 Time List", final_time_list)
#np.save("Day 2 10 RMS List", rms_list)
np.save("Day 2 10 MAXAMP List", max_amp_list)
#np.save("Day 2 10 Variance List", var_list)

"""
#total_blip_list plot
plt.scatter(i_list, blip_list, label="Number of events at given time", c = 'r')
plt.plot(i_list, total_blip_list, label="Total Events", c = 'g')
plt.xlabel("Time")
plt.ylabel("Event Number")
plt.show()
"""
