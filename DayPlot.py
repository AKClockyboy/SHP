import sys
import numpy as np
import obspy
from obspy import read
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import plot_trigger
import matplotlib.pyplot as plt


st = read("EC.RETU..SHZ.D.2013.mseed")

#Trace
tr = st[0]

tr_filt = tr.copy()

tr_filt.filter('bandpass', freqmin=2, freqmax=6, corners=2, zerophase=True)

#STA/LTA algorithms
cft = classic_sta_lta(.data, int(2 * df), int(10 * df))

plot_trigger(tr, cft, 1.5, 0.5)
