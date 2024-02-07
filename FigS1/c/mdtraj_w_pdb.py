#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:45:38 2024

@author: schaferjw

Create a dual-fold contact map for a fold-switching protein.
Two structures are needed to run this, the structures need to have indices that are similar (i.e. if one chain goes from 25-200 and the other goes from 1000-1175 this will fail)
This script will find an optimal alignment of the contact maps between the two structures.

Load in MSATransformer data in csv format to overlay the predictions onto the contacts derived from structures

"""
import re
from collections import Counter
import argparse
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

class PDB_STRUCTURE():
    """
    Get crystal structure contacts using mdtraj
    """
    def __init__(self,structure):
        print('\n')
        print('Loading crystallographic information.')
        # This section loads in the two pdbs and saves information needed for identifying intra-chain contacts
        self.pdb = md.load_pdb(f'{structure}')
        self.topo = self.pdb.topology
        table,bonds = self.topo.to_dataframe()
        self.indices = np.array(self.topo.select('protein and chainid 0')) #index of atom number
        self.pdb = self.pdb.atom_slice(self.indices)
        self.last_res = re.search('[0-9]+', str(self.topo.atom(self.indices[-1]))) #identify last residue
        self.first_res = re.search('[0-9]+', str(self.topo.atom(self.indices[1]))) #identify first residue

        #_______________________________________________________________________________________________________________________________

        self.N = len(table['resSeq'].unique())

        self.structure = structure

    def Get_contacts(self,distance):
        """
        read every contact that is within distance (set to 8A by default)
        this code finds all contacts by atom identifies residue numbers and removes redunant contacts from list
        this code also removes any contacts within 3 residues of origin
        """
        pairs,bfactor = [],[]
        for atom in zip(self.indices,self.pdb.xyz[0,self.indices,2]):
            contacts = md.compute_neighbors(self.pdb,cutoff=distance,query_indices=[atom[0]]) #contacts within 5A of index
            origin = re.search('[0-9]+',str(self.topo.atom(atom[0])))
            b = atom[1]
            for i in zip(contacts[0],self.pdb.xyz[0,contacts,2][0]):
                name = re.search('[0-9]+', str(self.topo.atom(i[0]))) #match the id of the sequence
                if int(name.group()) < int(origin.group())-2 or int(name.group()) > int(origin.group())+2:
                    if int(name.group()) > int(origin.group()):
                        if [int(origin.group()),int(name.group())] not in pairs:
                            pairs.append([int(origin.group()),int(name.group())])
                            bfactor.append(min(b,i[1]))

        return pairs
        #_______________________________________________________________________________________________________________________________

    def Multimer_contact(self,distance):
        """
        This section loads in the two pdbs and saves information needed for identifying inter-chain contacts
        """
        multi = md.load_pdb(f'{self.structure}')
        multi_topo = multi.topology
        multi_indices = np.array(multi_topo.select('protein and chainid 0')) #index of atom number
        last_res = re.search('[0-9]+', str(multi_topo.atom(multi_indices[-1]))) #identify last atom index number
        pairs = []
        chain_pairs = []

        #_______________________________________________________________________________________________________________________________

        #read every contact that is within distance (set to 10A by default)
        #this code finds all contacts by atom identifies residue numbers and removes redunant contacts from list
        #this code also removes any contacts within 3 residues of origin
        for atom in zip(multi_indices):
            contacts = md.compute_neighbors(multi,cutoff=distance,query_indices=[atom[0]]) #contacts within 10A of index
            origin = re.search('[0-9]+',str(multi_topo.atom(atom[0])))
            for i in zip(contacts[0]):
                name = re.search('[0-9]+', str(multi_topo.atom(i[0]))) #match the id of the sequence
                if int(name.group()) < int(origin.group())-4 or int(name.group()) > int(origin.group())+4 and int(name.group()) <int(last_res.group())+1:
                    if int(name.group()) > int(origin.group()):
                        if [int(origin.group()),int(name.group())] not in pairs:
                            pairs.append([int(origin.group()),int(name.group())])
                if i > multi_indices[-1]: #identify any residue number larger than origin chains largest index
                    if int(name.group()) > int(origin.group()):
                        if (int(origin.group()),int(name.group())) not in chain_pairs and (int(origin.group()),int(name.group())) not in pairs:
                            if len(name.group()) >3:
                                chain_pairs.append((int(origin.group()),int(name.group()[1:])))
                            else:
                                chain_pairs.append((int(origin.group()),int(name.group())))

        #_______________________________________________________________________________________________________________________________
        #create a separate  list of contacts within 10A that are between chain A and another chain (this will be empty in pdb's that only have one chain)
        chain_pairs = list(set(chain_pairs)) #using set is an efficient way of dealing with duplicates, note it does not retain the order of the original list
        chain_pairs = [self.reverse(i) if self.triangle(i[1],i[0]) > 0  else i for i in chain_pairs]
        chain_pairs = [[i[0],i[1]] for i in chain_pairs if i[1]-3 > i[0]]

        return chain_pairs
        #_______________________________________________________________________________________________________________________________


    def triangle(self,x,y):return y-x
    def reverse(self,x):return x[::-1]


class xOPT():
    """
   	search for optimal alignment of two pdb structures contacts
    """
    def __init__(self,contacts_pdb1,contacts_pdb2,pdb1N,pdb2N,structure_1,structure_2):
        print('Searching for best alignment of PDBs...')
        opt = [0,0]
        contacts_pdb1 = [[contact[0],contact[1]] for contact in contacts_pdb1]
        contacts_pdb2 = [[contact[0],contact[1]] for contact in contacts_pdb2]

        if pdb1N < pdb2N:
            for increment in range(-int(pdb1N/2),int(pdb1N/2)):

                a,b = [],[]
                pairs = [[i[0]+increment,i[1]+increment] for i in contacts_pdb1]
                for i in pairs:
                    if i in contacts_pdb2:
                        a.append('common')    #contact exists in both folds
                    else:
                        a.append(f'{structure_1}') #contact only exists in fold 1

                for i in contacts_pdb2:
                    if i in pairs:
                        b.append('common')          #contact exists in both folds
                    else:
                        b.append('{structure_2}') #contact only exists in fold 2
                concat = a + b
                count = Counter(concat) #creates a dictionary of each designation with the associated count
                if opt[0] < count['common']:
                    opt = (count['common'], increment) #optimal alignment is stored here as (number_of_both, number_to_add_to_index)
        else:
            for increment in range(-int(pdb2N/2),int(pdb2N/2)):
                a,b = [],[]
                pairs = [[i[0]+increment,i[1]+increment] for i in contacts_pdb2]
                for i in pairs:
                    if i in contacts_pdb1:
                        a.append('common')         #contact exists in both folds
                    else:
                        a.append(f'{structure_2}') #contact only exists in fold 1

                for i in contacts_pdb1:
                    if i in pairs:
                        b.append('common')     #contact exists in both folds
                    else:
                        b.append(f'{structure_1}') #contact only exists in fold 2
                concat = a + b
                count = Counter(concat) #creates a dictionary of each designation with the associated count
                if opt[0] < count['common']:
                    opt = (count['common'], increment) #optimal alignment is stored here as (number_of_both, number_to_add_to_index)

        print(f'Best alignment between pdbs found with adjustment of {opt[1]} yielding {opt[0]} redundant contacts.\n')

        self.opt = opt

    def OPT(self,manual,pdb,contacts_pdb1,contacts_pdb2,pdb1N,pdb2N,structure_1,structure_2):
        a,b = [],[]
        #add the perturbation to create the optimal alignment to the shorter pdb
        if manual == 'n':
            if pdb1N < pdb2N:
                contacts_pdb1 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in contacts_pdb1]
                contacts_pdb2 = [(i[0],i[1]) for i in contacts_pdb2]
            else:
                contacts_pdb2 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in contacts_pdb2]
                contacts_pdb1 = [(i[0],i[1]) for i in contacts_pdb1]

        else:
            if pdb1N < pdb2N:
                contacts_pdb1 = [(i[0]+int(manual),i[1]+int(manual)) for i in contacts_pdb1]
                contacts_pdb2 = [(i[0],i[1]) for i in contacts_pdb2]
            else:
                contacts_pdb2 = [(i[0]+int(manual),i[1]+int(manual)) for i in contacts_pdb2]
                contacts_pdb1 = [(i[0],i[1]) for i in contacts_pdb1]

        #_______________________________________________________________________________________________________________________________
        #recalculate the number of unique/nonunique contacts
        for i in contacts_pdb1:
            if i in contacts_pdb2:
                a.append('both')  #this means this contact exitst in both folds
            else:a.append(f'{structure_1}')  #this means this contact only exits in fold 1
        for i in contacts_pdb2:
            if i in contacts_pdb1:
                b.append('both')  #this means this contact exitst in both folds
            else:b.append(f'{structure_2}')  #this means this contact only exits in fold 2
        #_______________________________________________________________________________________________________________________________
        if manual == 'n':
            if pdb1N < pdb2N:
                contacts_pdb2 = [(i[1],i[0]) for i in contacts_pdb2]
            else:
                contacts_pdb1 = [(i[1],i[0]) for i in contacts_pdb1]
        else:
            if pdb1N < pdb2N:
                contacts_pdb2 = [(i[1],i[0]) for i in contacts_pdb2]
            else:
                contacts_pdb1 = [(i[1],i[0]) for i in contacts_pdb1]
        concat = a + b

        # save xcrystal information
        self.switch = 'n'
        if contacts_pdb1[0][0] < contacts_pdb1[0][1]:
            contacts_pdb1 = [(contact[1],contact[0]) for contact in contacts_pdb1]
            contacts_pdb2 = [(contact[1],contact[0]) for contact in contacts_pdb2]
            self.switch = 'y'

        plddt = ['none' for i in range(len(concat))]
        xcontact = contacts_pdb1 + contacts_pdb2
        df = pd.DataFrame(xcontact,columns=list('ij'))
        df['Fold'] = concat

        return xcontact, self.opt[1], df, plddt

    def OPT_Multi(self,dist_contacts_pdb1,dist_contacts_pdb2,manual,pdb,pdb1N,pdb2N,structure_1,structure_2):
        a,b = [],[]
       #add the perturbation to create the optimal alignment to the shorter pdb
        if manual == 'n':
            if pdb1N < pdb2N:
                dist_contacts_pdb1 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in dist_contacts_pdb1]
                dist_contacts_pdb2 = [(i[0],i[1]) for i in dist_contacts_pdb2]
            else:
                dist_contacts_pdb2 = [(i[0]+self.opt[1],i[1]+self.opt[1]) for i in dist_contacts_pdb2]
                dist_contacts_pdb1 = [(i[0],i[1]) for i in dist_contacts_pdb1]
        else:
            if pdb1N < pdb2N:
                dist_contacts_pdb1 = [(i[0]+int(manual),i[1]+int(manual)) for i in dist_contacts_pdb1]
                dist_contacts_pdb2 = [(i[0],i[1]) for i in dist_contacts_pdb2]
            else:
                dist_contacts_pdb2 = [(i[0]+int(manual),i[1]+int(manual)) for i in dist_contacts_pdb2]
                dist_contacts_pdb1 = [(i[0],i[1]) for i in dist_contacts_pdb1]

       #_______________________________________________________________________________________________________________________________
       #recalculate the number of unique/nonunique contacts
        for i in dist_contacts_pdb1:
            if i in dist_contacts_pdb2:
                a.append('both')  #this means this contact exitst in both folds
            else:
                a.append(f'{structure_1}')  #this means this contact only exits in fold 1
        for i in dist_contacts_pdb2:
            if i in dist_contacts_pdb1:
                b.append('both')  #this means this contact exitst in both folds
            else:
                b.append(f'{structure_2}')  #this means this contact only exits in fold 2
       #_______________________________________________________________________________________________________________________________
        if manual == 'n':
            if pdb1N < pdb2N:
                dist_contacts_pdb2 = [(i[1],i[0]) for i in dist_contacts_pdb2]
            else:
                dist_contacts_pdb1 = [(i[1],i[0]) for i in dist_contacts_pdb1]
        else:
            if pdb1N < pdb2N:
                dist_contacts_pdb2 = [(i[1],i[0]) for i in dist_contacts_pdb2]
            else:
                dist_contacts_pdb1 = [(i[1],i[0]) for i in dist_contacts_pdb1]
        concat = a + b

        # save xcrystal information
        # if dist_contacts_pdb1[0][0] < dist_contacts_pdb1[0][1]:
        if self.switch == 'y':
            dist_contacts_pdb1 = [(contact[1],contact[0]) for contact in dist_contacts_pdb1]
            dist_contacts_pdb2 = [(contact[1],contact[0]) for contact in dist_contacts_pdb2]

        plddt = ['none' for i in range(len(concat))]

        # save xcrystal information
        xcontact = dist_contacts_pdb1 + dist_contacts_pdb2
        df = pd.DataFrame(xcontact,columns=list('ij'))
        df['Fold'] = concat
        return xcontact, df

def FS_plot(pdb1,pdb2,df_dist,df,df_sorted):
    """
    Create a dual-fold contact map
    """
    colors = {'both':'#767676'}
    counts = df_sorted['sort'].value_counts()
    if 'pdb_1' not in counts:
        counts['pdb_1']  = 0
    if 'pdb_2' not in counts:
        counts['pdb_2']  = 0
    if counts['pdb_1'] > counts['pdb_2']:
        df_sorted = df_sorted.rename({'i':'j', 'j':'i'}, axis=1)
        df_dist = df_dist.rename({'i':'j', 'j':'i'}, axis=1)
        df = df.rename({'i':'j', 'j':'i'}, axis=1)
        colors = {'both':'#767676', f'{pdb1}':'#682860', f'{pdb2}':'#000000'}
        name1 = pdb2
        name2 = pdb1
    if counts['pdb_1'] <= counts['pdb_2']:
        colors = {'both':'#767676', f'{pdb2}':'#682860', f'{pdb1}':'#000000'}
        name1 = pdb1
        name2 = pdb2

    df_plot = df_sorted.copy()
    f, ax = plt.subplots(1,1,figsize=(11,9))

    if df_dist.empty is False:
        ax.scatter(x=df_dist['i'], y=df_dist['j'], s=30, c=df_dist['Fold'].map(colors),linewidth=0)
    ax.scatter(x=df['i'], y=df['j'],c=df['Fold'].map(colors), s=100, linewidth=0, label=df['Fold'])

    df_plot = df_plot.drop(df_plot[df_plot.group == -1].index)
    df_other = df_plot.loc[df_plot["sort"] == 'noise']
    df_temp = df_other.rename({'i':'j', 'j':'i'}, axis=1)
    df_other = pd.concat([df_other,df_temp])
    df_plot = df_plot.drop(df_plot[df_plot.sort == 'noise'].index)
    df_c = df_plot.loc[df_plot["sort"] == 'common']
    df_c = df_c.rename({'i':'j', 'j':'i'}, axis=1)
    df_s = df_plot.loc[df_plot["sort"] == 'pdb_1']
    df_s = df_s.rename({'i':'j', 'j':'i'}, axis=1)
    df_plot = df_plot.drop(df_plot[df_plot.sort == 'pdb_1'].index)
    df_plot = pd.concat([df_plot,df_s,df_c])

    #redefine colormap
    cmap = cm.get_cmap('Blues')
    blues = cm.get_cmap('Blues', 128)
    ncmap = blues(np.linspace(0, 1.0, 128))
    ncmap = ncmap[0:85]
    ncmap = ListedColormap(ncmap)

    high_contacts = df_plot[df_plot['zscore'] > 0.5]

    print(high_contacts.to_string())

    df_plot = df_plot.sort_values('zscore') #make sure largest values are plotted last so that they aren't buried 
    plot = ax.scatter(x=df_plot['i'], y=df_plot['j'], alpha=df_plot['zscore'], s=30, linewidth=0, c=df_plot['zscore'], cmap=ncmap, vmin=-0.0, vmax=1.0)
    plt.colorbar(plot)
    plt.plot([115,115],[0,162],'r--')
    plt.plot([0,162],[115,115],'r--')
    df_other = df_other.drop_duplicates()
    ax.scatter(x=df_other['i'], y=df_other['j'],c='#008080', s=20,edgecolor='None',marker='D', alpha=0.3)

    #_______________________________________________________________________________________________________________________________
    #Final adjustments to image
    ax.set_xlabel(r'RfaH $\beta$-sheet', fontsize=30)
    ax.set_ylabel(r'RfaH $\alpha$-helix', fontsize=30)
    plt.xticks(fontsize=11, weight = 'bold')
    plt.yticks(fontsize=11, weight = 'bold')

    #_______________________________________________________________________________________________________________________________
    plt.savefig('df_cmap.png',dpi=200)

def Import_msatransformer(csv):
    """
    MSATransformer predictions loaded in as csv file of NxN matrix
    returned as a pandas data frame to be plotted
    """
    df_sorted = {'i':[],'j':[],'sort':[],'zscore':[], 'group':[]}
    mtx = np.loadtxt(f"{csv}", delimiter=",", dtype=float)
    for row_idx, i in enumerate(mtx):
        for col_idx, j in enumerate(i):
            df_sorted['i'].append(row_idx)
            df_sorted['j'].append(col_idx)
            df_sorted['zscore'].append(j)
            df_sorted['sort'].append('')
            df_sorted['group'].append(0)
    df_sorted = pd.DataFrame().from_dict(df_sorted)
    return df_sorted

def main():

    """
    compare pdb structures to contact predictions from MSATransformer
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--Manual_PDB", type=str, default='n', help='Turn off auto-align algorithm and manually align PDBs')
    parser.add_argument("--pdb1",  type=str, help='PDB structure for the first conformation (remove hydrogens, ligands, and HOH')
    parser.add_argument("--pdb2",  type=str, help='PDB structure for the second conformation (remove hydrogens, ligands, and HOH')
    parser.add_argument("--MSA_T", type=str, help='MSATransformer prediction data in csv format')

    args = parser.parse_args()

    #___________________________________________________________________________________________________
    #   Load crystallagraphic information
    #___________________________________________________________________________________________________
    #MDtraj uses deprecated np.int instead of the new np.int32 or np.int64
    #comment these lines out to see full deprecation warning
    import warnings
    warnings.filterwarnings('ignore')

    pdb1 = PDB_STRUCTURE(args.pdb1)
    contacts_pdb1 = pdb1.Get_contacts(0.8)       #define list of intrachain contacts for comparison
    dist_contacts_pdb1 = pdb1.Multimer_contact(1.0)   #create list of inerchain contacts (based on above)
    pdb2 = PDB_STRUCTURE(args.pdb2)
    contacts_pdb2 = pdb2.Get_contacts(0.8)       #define list of intrachain contacts for comparison
    dist_contacts_pdb2 = pdb2.Multimer_contact(1.0)   #create list of inerchain contacts (based on above)

    xopt = xOPT(contacts_pdb1,contacts_pdb2,pdb1.N,pdb2.N,args.pdb1,args.pdb2)
    xcontact, opt, df_intra, plddt = xopt.OPT(args.Manual_PDB,args.pdb1[:-4],contacts_pdb1,contacts_pdb2,pdb1.N,pdb2.N,args.pdb1,args.pdb2)
    dist_xcontact, df_dist = xopt.OPT_Multi(dist_contacts_pdb1,dist_contacts_pdb2,args.Manual_PDB,args.pdb1[:-4],pdb1.N,pdb2.N,args.pdb1,args.pdb2)

    #___________________________________________________________________________________________________
    #   Load MSATransformer information
    #___________________________________________________________________________________________________

    df_sorted = Import_msatransformer(args.MSA_T)

    #Create final plot
    FS_plot(args.pdb1, args.pdb2, df_dist, df_intra, df_sorted)

if __name__ == '__main__':
    main()
