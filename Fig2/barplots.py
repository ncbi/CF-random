from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np


def make_fig(proteins,att,ylabel='',name='test'):

    x = np.arange(len(proteins))
    width = 0.25
    spacing = 0.05
    multiplier = 0
    c = ['white','#555555']

    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size
    mpl.rcParams['figure.figsize'] = [3,4]
    if name == 'num_preds':
        mpl.rcParams['figure.figsize'] = [3,3.3]

    plt.rcParams['font.family'] = 'Helvetica'

    fig, ax = plt.subplots(layout='constrained')

    ax.set_prop_cycle(color=c)

    for attribute, measurement in att.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute,edgecolor='k')
        ax.bar_label(rects, padding=2,fontsize=8)
        multiplier += 1

    ax.set_ylabel(ylabel,fontsize=16)
    ax.set_xticks(x+width,proteins)

    plt.savefig(name+'.png',dpi=600,transparent=True)

    

if __name__ == '__main__':

    proteins = ('KaiB','Mad2','RfaH')

    MCCs = {
        'CF-random': (1.0,0.87,0.62),
        'AF-cluster': (0.6,0.7,0.4)}

    Accuracies = {
        'CF-random': (68,24,64),
        'AF-cluster': (29,2,4)}

    Efficiencies = {
        'CF-random': (2,1,2),
        'AF-cluster': (329,95,250)}

    Structures = {
        'CF-random': (330,95,250),
        'AF-cluster': (329,95,250)}

    make_fig(proteins,MCCs,'Matthews Correlation Coefficient','MCC')
    make_fig(proteins,Accuracies,'%Success','accuracy')
    make_fig(proteins,Efficiencies,'Number of AF2/CF Runs','efficiency')
    make_fig(proteins,Structures,'#Predictions','num_preds')
    
    

    

    
