import numpy as np
from scipy import signal

def calc_norm_spec(data,fs):
    try:
        x = data
        f, t, Sxx = signal.spectrogram(x, fs = fs, nperseg = 256, noverlap= 128)
        datamin = Sxx.min(axis = 0)
        datamax = Sxx.max(axis = 0)
        datarange = datamax - datamin
        dataNorm = (Sxx-datamin) / datarange
        medianPower = np.nanmedian(dataNorm, axis = 1)
    except Exception as e:
        print(e)
        medianPower = np.nan
        f= np.nan
    return(medianPower, dataNorm, f)

def mean_freq_ratio(median_freq,f, low_f,high_f,init_f,final_f):
    try:
        ratio=np.mean(median_freq[(f>low_f) & (f<high_f)])
    except:
        ratio=np.nan
    return ratio

def split_by_0(arr,cutoff,above_below, tolerance):
    equals = arr
    l = []
    t=0
    if above_below=='below':
        for i in range(len(arr)):
            if arr[i]<cutoff:
                l.append(i)
                t=0
            if arr[i]>=cutoff and len(l)>0:
                l.append(i)
                t+=1
                if t>=tolerance:
                    yield l[:-tolerance]
                    l = []
                    t=0
        yield l
    if above_below=='above':
        for i in range(len(arr)):
             if arr[i]>=cutoff:
                 l.append(i)
                 t=0
             if arr[i]<cutoff and len(l)>0:
                 l.append(i)
                 t+=1
                 if t>=tolerance:
                     yield l[:-tolerance]
                     l = []
                     t=0
        yield l
