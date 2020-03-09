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

#
#Getting stream using select
#
st2 = st.select(station = "MBGE", component = 'Z')
#print(st2.__str__(extended=True))

#
#Selecting start time and slicing, so we only look at one days worth of data
#
dt = obspy.UTCDateTime("1997-02-11T00:00:00")
print("dt is: " + str(dt))
ft = dt + 86400
print("ft is: " + str(ft))
dt_in_seconds = dt.second
st2 = st2.slice(starttime = dt, endtime = ft)

#
#list and counters
#

max_amp_list = []
rms_list = []
var_list = []
final_time_list = []

total_blip_list = [] #Total number of events at each point
blip_list = [] #Number of events in a given window of time
i_list = [] #number of times we loop

second_count = 0 #counter for seconds
number_of_blips = 0

#
#looping ove
#
for i in range(len(st2)):

    tr2 = st2[i]
    tr_filt = tr2.copy()
    df = tr_filt.stats.sampling_rate

    i_list.append(i*1200)

    #Spectrogram Plotter
    #st2.plot()
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

    #Getting the trigger
    trig = classic_sta_lta(tr_filt.data, int(5 * df), int(10 * df))
    #plot_trigger(tr_filt, trig, 1.6, 0.65)

    #Getting a list of trigger start and end times
    time_list_on = []
    time_list_off = []
    time_list = (obspy.signal.trigger.trigger_onset(trig, 1.6, 0.65, max_len=9e+99))



    for i in range(len(time_list)):
        time_list_on.append(time_list[i][0])
        time_list_off.append(time_list[i][1])

    blip_list.append(len(time_list_on))
    number_of_blips += len(time_list_on)
    total_blip_list.append(number_of_blips)

    for i in range(len(time_list_on)):

        dt2 = (dt + time_list_on[i]/df)
        ft2 = (dt + time_list_off[i]/df)

        second_count += (dt2.second)*0.67

        final_time_list.append(second_count)
        st3 = st2.slice(starttime = dt2, endtime = ft2)
        tr3 = st3[0]

        waveform_data = (tr3.data) #Check that this is the amplitude

        max_amp_list.append(abs(tr3.max()))
        rms_list.append(np.sqrt(np.mean(waveform_data**2)))
        var_list.append(np.var(waveform_data))

    dt = dt2
    ft = ft2 #resetting the times

"""#RMS Values Plotting
plt.scatter(final_time_list, rms_list, c ='crimson')
plt.xlabel("Time in Seconds from _____")
plt.ylabel("RMS")
plt.show()

#Variance Plotting
plt.scatter(final_time_list, var_list, c ='crimson')
plt.xlabel("Time in Seconds from _____")
plt.ylabel("Var")
plt.show()

#RMS vs Max Amp
plt.scatter(rms_list, max_amp_list, c = 'crimson')
plt.xlabel("RMS")
plt.ylabel("Maximum amplitude")
plt.show()
"""
#total_blip_list plot
plt.scatter(i_list, blip_list, c = 'r', label = "Number of events at given time")
plt.plot(i_list, total_blip_list, c = 'g', label = "Total Events")
plt.xlabel("Time")
plt.ylabel("Event Number")
plt.show()
