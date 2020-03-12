import sys
import numpy as np
import math

import os

import obspy
from obspy.core import read
from numpy import (absolute, int16)
from scipy.io import wavfile
from obspy import UTCDateTime
from obspy import read_inventory
from obspy.io.xseed import Parser



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
station = "MBGE"
st2 = st.select(station = "MBGE", component = 'Z')
#print(st2.__str__(extended=True))


#
#Selecting start time and slicing, so we only look at one days worth of data
#
dt = obspy.UTCDateTime("1997-02-13T00:00:00")
ft  = obspy.UTCDateTime("1997-02-14T00:00:00")


st2 = st2.slice(starttime = dt, endtime = ft)
print(st2)
st2.merge(fill_value = 'interpolate')

str3=st2.copy()
str3.detrend('demean')



fs = st2[0].stats.sampling_rate

xspd  = 180.                               # speed up by factor of 360 (1hr = 10 seconds)
fswav = fs*xspd
wscal = str3[0].data/max(absolute(str3[0].data))    # scale to -1 <= p <= 1; p is vector length M
wscal = 0.99*wscal                         # avoid clipping
convert_16_bit = float(2**15)              # convert to int16 format for wavwrite
wscal = int16( wscal * convert_16_bit ) # convert to int16 format for wavwrite
wavfile.write("Day 3 Sound.wav", int(fswav), wscal)

print('Audio file written')

del wscal
