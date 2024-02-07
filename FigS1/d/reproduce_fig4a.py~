#! /usr/local/bin/python3

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import seaborn as sns
sns.set_style('ticks')
sns.set_context('paper')

if __name__ == '__main__':

    df = pd.read_json('Rfah_3recycles.json')

    def get_type(pdb):
        end = pdb.split('_')[-1]
        if end=='REF.pdb':
            return 'REF'
        elif end.startswith('U'):
            return end.split('-')[0]
        else:
            return 'Tree'

    df_confident = df[df['mean_pLDDT'] > 80]

    print(df_confident['mean_pLDDT'],df_confident['pdb'],df_confident['PC 1'],df_confident['PC 2'])
    
    df['Type'] = df.apply(lambda row: get_type(row['pdb']), axis=1)

    df = df.loc[df.Type=='Tree']

    plt.figure(figsize=(4,4))
    df = df.sort_values('mean_pLDDT')
    plt.scatter(df['PC 1'], df['PC 2'], c=df['mean_pLDDT'], cmap='rainbow_r',vmin=50, vmax=90)

    plt.xlabel('PC 1 (on contacts)')
    plt.ylabel('PC 2 (on contacts)')

    plt.savefig('rfah_landscape.pdf',bbox_inches='tight')
