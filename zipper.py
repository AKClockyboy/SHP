import sys
import numpy as np
import obspy
import os
from obspy import read
import matplotlib.pyplot as plt

st = obspy.Stream()
cwd = os.getcwd()

files = os.listdir(cwd)
print(files)
print("length is " + str(len(files)))

for filename in os.listdir(cwd):
    if filename.endswith('.htm'):
        st += read(filename)

#Trace
tr = st[0]

df = tr.stats.sampling_rate

print("Trace is: " + str(st[0]))
print(st[0].stats)
print(st[0].data)

tr.plot()

#Day Plot of this thing
#st.plot(type="dayplot", interval=60, right_vertical_labels=False, vertical_scaling_range=5e3, one_tick_per_line=True, color=['k', 'r', 'b', 'g'], show_y_UTC_label=False, events={'min_magnitude': 7})
