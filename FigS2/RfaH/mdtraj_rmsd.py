from glob import glob
from typing import List
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import PDB
import sys, math

plt.rcParams['font.family'] = 'Helvetica'


def Bfactor(pdb: str) -> List:
    """
    take in a pdb file and identify the index of every alpha carbon
    """
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb)

    _bfactor = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    _bfactor.append(residue['CA'].bfactor)
    return _bfactor


def Alpha_Carbon_Indices(pdb: str) -> List:
    """
    take in a pdb file and identify the index of every alpha carbon
    """
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb)

    alpha_carbons = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    resid = residue.resname
                    alpha_carbons.append([resid, residue['CA'].get_serial_number() - 1])
    return alpha_carbons

def Match_Alpha_Carbons(pdb_1: str, pdb_2: str) -> List[int]:
    """
    Take in two pdb structure files and search through them for matching alpha carbons
    This should identify positions correctly even if sequences are not identical
    """
    alpha_c_1 = Alpha_Carbon_Indices(pdb_1)
    alpha_c_2 = Alpha_Carbon_Indices(pdb_2)

    matching_alpha_carbons1 = []
    matching_alpha_carbons2 = []

    for i, (resname_1, ca_index1) in enumerate(alpha_c_1):
        for j, (resname_2, ca_index2) in enumerate(alpha_c_2):
            if resname_2 == resname_1 and ca_index1 not in [_[1] for _ in matching_alpha_carbons1] and ca_index2 not in [_[1] for _ in matching_alpha_carbons2]:
                #prevent erroneous match at NTD
                if i > 0 and j > 0:
                    if alpha_c_1[i-1][0] != alpha_c_2[j-1][0]: #check previous matches
                        continue
                # prevent erroneous backtracking
                if len(matching_alpha_carbons1) > 2 and len(matching_alpha_carbons2) > 2:
                    if ca_index2 < matching_alpha_carbons2[-1][-1]:
                        continue
                #prevent erroneous match at CTD
                if i < len(alpha_c_1) - 1 and j < len(alpha_c_2) - 1:
                    if alpha_c_1[i+1][0] != alpha_c_2[j+1][0]: #check next matches
                        continue

                matching_alpha_carbons1.append([resname_1, ca_index1])
                matching_alpha_carbons2.append([resname_2, ca_index2])
                break
    #skip first residue to avoid erroneous glycine match
    return matching_alpha_carbons1[1:], matching_alpha_carbons2[1:]


def Calculate_RMSD(structure_1: str, structure_2: str, structure_1_index: List[int], structure_2_index: List[int]) -> int:
    """
    calculate the RMSD between two structures using MDtraj library
    this script will fail if mdtraj is not loaded in your python environment
    recommend python 3.10
    """

    #with warnings.catch_warnings(action="ignore"):
    #    turn_off_warnings()
    #load structure information in mdtraj
    pdb = md.load(structure_1)
    pdb_ca = pdb.atom_slice(structure_1_index) #select only CA atoms

    #load structure information in mdtraj
    reference = md.load(structure_2)
    reference_ca = reference.atom_slice(structure_2_index) #select only CA atoms

    # Calculate RMSD of CA atoms
    pdb_ca.superpose(reference_ca)
    return md.rmsd(pdb_ca, reference_ca, frame=0)

def MCC(tp,fp,tn,fn):

    return(((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))


if __name__ == '__main__':

    reference_6c6s = '6c6s.pdb'
    reference_2oug = '2oug.pdb'

    ###Comment top pdb_structures for CF-random figure, bottom for AF-cluster 
    pdb_structures = glob('00_RfaH_preds_3r/RFAH_[0-9]??*pdb') #AF-cluster
    #pdb_structures = glob('RfaH_250/*/*pdb') #CF-random


    results_dict = {'2oug rmsd': [], '2oug rmsd CTD': [], '6C6S rmsd': [], '6C6S rmsd CTD':[], 'bfactor': [], 'structure': []}
    rmsd_via_2oug = []
    bfactor = []

    alpha = 0
    beta = 0
    other = 0
    fp = 0
    fn = 0
    tn = 0

    if len(sys.argv) < 2:
        sys.exit('Usage: python mdtraj_rmsd.py FigureName')
    
    for structure in pdb_structures:
        structure_index,reference_index = Match_Alpha_Carbons(structure, reference_2oug)
        results_dict['structure'].append(structure)
        results_dict['2oug rmsd'].append(Calculate_RMSD(structure, reference_2oug, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0])
        results_dict['2oug rmsd CTD'].append(Calculate_RMSD(structure, reference_2oug, [_[1] for _ in structure_index][103:],[_[1] for _ in reference_index][103:])[0])
        results_dict['bfactor'].append(np.average(Bfactor(structure)[103:]))

        structure_index,reference_index = Match_Alpha_Carbons(structure, reference_6c6s)
        results_dict['6C6S rmsd'].append(Calculate_RMSD(structure, reference_6c6s, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0])
        results_dict['6C6S rmsd CTD'].append(Calculate_RMSD(structure, reference_6c6s, [_[1] for _ in structure_index][103:],[_[1] for _ in reference_index][103:])[0])

        if results_dict['2oug rmsd CTD'][-1]*10 <= 3 and results_dict['bfactor'][-1] >= 55:
            alpha += 1
        elif results_dict['6C6S rmsd CTD'][-1]*10 <= 3 and results_dict['bfactor'][-1] >= 55:
            beta += 1
        else:
            other += 1

        if results_dict['2oug rmsd CTD'][-1]*10 <= 3 and results_dict['bfactor'][-1] < 55:
            fn += 1
        elif results_dict['6C6S rmsd CTD'][-1]*10 <= 3 and results_dict['bfactor'][-1] < 55:
            fn += 1

        if results_dict['2oug rmsd CTD'][-1]*10 >3 and results_dict['6C6S rmsd CTD'][-1]*10 >3 and results_dict['bfactor'][-1] >= 55:
            fp += 1
        elif results_dict['2oug rmsd CTD'][-1]*10 >3 and results_dict['6C6S rmsd CTD'][-1]*10 >3 and results_dict['bfactor'][-1] < 55:
            tn += 1
            
    df = pd.DataFrame.from_dict(results_dict)
    df['2oug rmsd']     = df['2oug rmsd']*10 # nm -> A
    df['2oug rmsd CTD'] = df['2oug rmsd CTD']*10 # nm -> A
    df['6C6S rmsd']     = df['6C6S rmsd']*10 # nm -> A
    df['6C6S rmsd CTD'] = df['6C6S rmsd CTD']*10 # nm -> A

    df = df.sort_values('bfactor')

    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #    print(df)
    #for idx,_ in enumerate(df['structure']):
    #    print(idx,_)

    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size
    plt.scatter(data=df, x='6C6S rmsd CTD',y='2oug rmsd CTD', c=df['bfactor'], cmap='rainbow_r',vmin=40,vmax=90)
    plt.ylabel(r'Autoinhibited RMSD ($\AA$)',fontsize=16)
    plt.xlabel(r'Active RMSD ($\AA$)',fontsize=16)
    #pltxlim((0,int(max(df['6C6S rmsd CTD'].max(),df['2oug rmsd CTD'].max()) + 1.5)))
    #plt.ylim((0,int(max(df['6C6S rmsd CTD'].max(),df['2oug rmsd CTD'].max()) + 1.5)))
    plt.xlim([0,18])
    plt.ylim([0,13])
    
    plt.colorbar()
    plt.savefig(sys.argv[1]+'.png',dpi=600)

    s = alpha+beta+other

    print('alpha: %1.4f, beta: %1.4f, other: %1.4f, all: %i' %(alpha/s,beta/s,other/s,s))

    d = fp+tn+fn+alpha+beta
    tp = alpha+beta

    print('TP: %1.2f, FP: %1.2f, TN: %1.2f, FN: %1.2f, all: %i, tp: %i' %(tp/d,fp/d,tn/d,fn/d,d,tp))

    print('MCC: %1.2f' %(MCC(tp,fp,tn,fn)))
