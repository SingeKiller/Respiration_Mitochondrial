#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

datas = [np.load("data14h16.npy"),
         np.load("data15h12.npy"),
         np.load("data15h37.npy"),
         np.load("data16h00.npy"),
         np.load("data17h03.npy")]
labels = ["data14h16.npy",
          "data15h12.npy",
          "data15h37.npy",
          "data16h00.npy",
          "data17h03.npy"]

fs = 12

plt.ion()
plt.close('all')
plt.figure()
plt.xlabel("Time (min)",fontsize=fs)
plt.ylabel("Normalized oxygen concentration",fontsize=fs)

for data,label in zip(datas,labels):

  plt.plot(data[0],data[1],label=label)
plt.gca().tick_params(labelsize=0.8*fs)
plt.gca().tick_params(labelsize=0.8*fs)
plt.tight_layout()
plt.legend()

plt.savefig("resp.svg")

