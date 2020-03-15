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
dt = obspy.UTCDateTime("1997-02-13T00:00:00")
print("dt is: " + str(dt))
ft = dt + 86400
print("ft is: " + str(ft))
dt_in_seconds = dt.second
st2 = st2.slice(starttime = dt, endtime = ft)
st2.split()
#
#list and counters
#

final_time_list = []

total_blip_list = [] #Total number of events at each point
blip_list = [] #Number of events in a given window of time
i_list = [] #number of times we loop

second_count = 0 #counter for seconds
number_of_blips = 0

#
#looping over
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

    #
    #Getting the trigger
    #
    trig = classic_sta_lta(tr_filt.data, int(5 * df), int(10 * df))
    #plot_trigger(tr_filt, trig, 1.6, 0.65)

    #
    #Getting a list of trigger start and end times
    #
    time_list_on = []
    time_list_off = []
    time_list = (obspy.signal.trigger.trigger_onset(trig, 1.62, 0.65, max_len=9e+99))

    for i in range(len(time_list)):
        time_list_on.append(time_list[i][0])


    blip_list.append(len(time_list_on))
    number_of_blips += len(time_list_on)
    total_blip_list.append(number_of_blips) #Counting number of blips


#total_blip_list plot
fig,ax = plt.subplots()
plt.title("Timeseries")
ax.bar(i_list, blip_list, width = 1200, align = 'edge', color = 'purple')
ax.set_xlabel("time",fontsize=14)
ax.set_ylabel("Number of Events per 20min",color="purple",fontsize=14)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
ax2.plot(i_list, total_blip_list, label="Total Events", c = 'k')
ax2.set_ylabel("Total Number of Events",color="k",fontsize=14)
plt.show()
