import numpy as np
from scipy.ndimage.filters import convolve, gaussian_filter
import matplotlib.pyplot as plt
import math
import os
import re


ArcToCm = math.pi / 180.0 / 3600.0 * 1.49597870e13
frequency = [
    1000000000,
    2000000000,
    3750000000,
    9400000000,
    17000000000,
    35000000000,
    55000000000,
    80000000000,
    110000000000,
    140000000000,
    190000000000,
    230000000000,
    250000000000
]

prefix = 'c7DataHR_pol'
title = 'C7 Atmosphere High-Res'
numVox = 32
VoxToCm = 181317735.10211

files = [x for x in os.listdir() if x.endswith('.csv') and x.startswith(prefix)]
indices = [int(re.match(prefix+'_(.*?).csv', x).group(1)) for x in files]
files = [x for _,x in sorted(zip(indices, files))]

for i in range(len(files)):
    mat = np.genfromtxt(files[i], delimiter=',').T
    beamSize = 0
    if frequency[i] < 35e9:
        beamSize = 10 
    elif frequency[i] < 100e9:
        beamSize = 5
    else:
        beamSize = 2

    fovArcSec = numVox * VoxToCm / ArcToCm
    pxArcSec = fovArcSec / mat.shape[0]
    kernelSize = beamSize / pxArcSec
    conv = gaussian_filter(mat, kernelSize, mode='constant', cval=0.0)
    plt.figure(figsize=(10,8))
    plt.imshow(conv, origin='bottom left', extent=[0, 1, 0, 1], cmap='plasma')
    plt.title("Total Brightess Temperature at %.2f GHz (%.1f' beam, %s)" % (frequency[i] / 1e9, beamSize, title))
    locs, labels = plt.xticks()
    labels = ['%.0f' % x for x in np.linspace(0, fovArcSec, num=len(locs))]
    plt.xticks(locs, labels)
    plt.yticks(locs, labels)
    # plt.xlabel("'")
    # plt.ylabel("'")
    plt.colorbar()
    plt.tight_layout()
    filename = 'TotTbConv_'
    if prefix.endswith('pol'):
        filename = 'PolConv_'
    plt.savefig(filename+prefix+'_'+str(i)+'.png')
    # plt.show()

