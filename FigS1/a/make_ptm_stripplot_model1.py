#! /usr/local/bin/python3.10

import sys
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd

font = {'size'   : 12}

matplotlib.rc('font', **font)

if __name__ == '__main__':

    f = open('af2_ptm_scores.txt').read().splitlines()
    f2 = open('good_RFAH_autoinhibited_plDDTS_curated.txt').read().splitlines()

    nums = []
    numsm = []
    nums2 = []

    l70 = 0
    g80 = 0

    for i in f:
        info = i.split()
        if 'model_1' in info[0]:
            numsm.append(float(info[1][:-1]))
        nums.append(float(info[1][:-1]))
        if nums[-1] < 70:
            l70 += 1
        elif nums[-1] >= 80:
            g80 += 1

    for i in f2:
        info = i.split()
        nums2.append(float(info[1][:-1]))
        
    print(min(nums))
    

    dict_info = {'Whole MSA':nums,'Whole MSA, model 1': numsm, 'AF-cluster':nums2}

    labels = []
    for i in range(len(nums)):
        labels.append('Whole MSA')
    for i in range(len(numsm)):
        labels.append('Whole MSA\nmodel 1')
    for i in range(len(nums2)):
        labels.append('AF-cluster')


    all_nums = nums+numsm+nums2

    dat = {'Method':labels,'pLDDT':all_nums}

    fig= plt.figure(figsize=(2.2,4))
    ax = fig.add_subplot(111)

    sns.stripplot(data=dat,x='Method',y='pLDDT',zorder=0)
    plt.plot([-0.15,0.15],[68.6,68.6],'r-',zorder=10)
    plt.setp(ax.get_xticklabels(), ha="right", rotation=45)


    avg = np.average(nums)
    med = np.median(nums)

    avg2=np.average(numsm)

    avg3 = np.average(nums2)

    
    print(avg)
    print(avg2)
    print(avg3)

    plt.plot([-0.15,0.15],[avg,avg],'k-',zorder=10)
    plt.plot([0.85,1.15],[avg2,avg2],'k-',zorder=10)
    plt.plot([1.85,2.15],[avg3,avg3],'k-',zorder=10)
    plt.ylabel('plDDT')

    plt.tight_layout()

    plt.savefig('ptm_plot_all_and_model1.png',dpi=1200)
