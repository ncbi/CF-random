from glob import glob
from typing import List
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import PDB
import sys, math
import matplotlib as mpl

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
    
    #Return Matthews Correlation Coefficient
    return(((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))


if __name__ == '__main__':

    reference_1s2h = '1s2h.pdb'
    reference_3gmh = '3GMH_L.pdb'

    ###Comment top pdb_structures for af-cluster predictions; bottom for cf-random
    pdb_structures = glob('Mad2_8_16_95_5model/*pdb') #cf random
    #pdb_structures = glob('01_Mad2_1S2H_preds/1S2H_[0-9]??.pdb') #af_cluster

    if len(sys.argv) < 2:
        sys.exit('Usage: python mdtraj_rmsd.py FigureName')

    results_dict = {'3GMH rmsd': [], '3GMH rmsd CTD': [], '1S2H rmsd': [], '1S2H rmsd CTD':[], 'bfactor': [], 'structure': []}
    rmsd_via_3gmh = []
    bfactor = []

    openC = 0
    closedC = 0
    other = 0
    fp = 0
    fn = 0
    tn = 0
    
    
    for structure in pdb_structures:
        structure_index,reference_index = Match_Alpha_Carbons(structure, reference_3gmh)
        results_dict['structure'].append(structure)
        results_dict['3GMH rmsd'].append(Calculate_RMSD(structure, reference_3gmh, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0])
        results_dict['3GMH rmsd CTD'].append(Calculate_RMSD(structure, reference_3gmh, [_[1] for _ in structure_index][0:],[_[1] for _ in reference_index][0:])[0])
        results_dict['bfactor'].append(np.average(Bfactor(structure)[0:]))

#    rmsd_via_1s2h = []
#    for structure in pdb_structures:
        structure_index,reference_index = Match_Alpha_Carbons(structure, reference_1s2h)
        results_dict['1S2H rmsd'].append(Calculate_RMSD(structure, reference_1s2h, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0])
        results_dict['1S2H rmsd CTD'].append(Calculate_RMSD(structure, reference_1s2h, [_[1] for _ in structure_index][0:],[_[1] for _ in reference_index][0:])[0])

        if results_dict['3GMH rmsd'][-1]*10 <= 5.5 and results_dict['bfactor'][-1] >= 70:
            openC += 1
        elif results_dict['1S2H rmsd'][-1]*10 <= 5.5 and results_dict['bfactor'][-1] >= 70:
            closedC += 1
        else:
            other += 1

        if results_dict['3GMH rmsd'][-1]*10 <= 5.5 and results_dict['bfactor'][-1] < 70:
            fn += 1
        elif results_dict['1S2H rmsd'][-1]*10 <= 5.5 and results_dict['bfactor'][-1] < 70:
            fn += 1

        if results_dict['3GMH rmsd'][-1]*10 > 5.5 and results_dict['1S2H rmsd'][-1]*10 > 5.5 and results_dict['bfactor'][-1] >= 70:
            fp += 1
        elif results_dict['3GMH rmsd'][-1]*10 > 5.5 and results_dict['1S2H rmsd'][-1]*10 > 5.5 and results_dict['bfactor'][-1] < 70:
            tn += 1

    df = pd.DataFrame.from_dict(results_dict)
    df['3GMH rmsd']     = df['3GMH rmsd']*10 # nm -> A
    df['3GMH rmsd CTD'] = df['3GMH rmsd CTD']*10 # nm -> A
    df['1S2H rmsd']     = df['1S2H rmsd']*10 # nm -> A
    df['1S2H rmsd CTD'] = df['1S2H rmsd CTD']*10 # nm -> A

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)
    for idx,_ in enumerate(df['structure']):
        print(idx,_)

    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size

    plt.scatter(data=df, x='1S2H rmsd CTD',y='3GMH rmsd CTD', c=df['bfactor'], cmap='rainbow_r', vmin=50,vmax=90)
    plt.xlabel(r'Closed RMSD ($\AA$)',fontsize=16)
    plt.ylabel(r'Open RMSD ($\AA$)',fontsize=16)
    #plt.xlim((0,int(max(df['1S2H rmsd CTD'].max(),df['3GMH rmsd CTD'].max()) + 1.5)))
    #plt.ylim((0,int(max(df['1S2H rmsd CTD'].max(),df['3GMH rmsd CTD'].max()) + 1.5)))
    plt.xlim([0,23])
    plt.ylim([0,21])
    plt.colorbar()
    plt.savefig(sys.argv[1],transparent=True,dpi=600)

    s = openC+closedC+other

    print('open: %1.2f, closed: %1.2f, other: %1.2f, all: %i' %(openC/s,closedC/s,other/s,s))

    d = fp+tn+fn+openC+closedC
    tp = openC+closedC

    print('TP: %1.2f, FP: %1.2f, TN: %1.2f, FN: %1.2f, all: %i' %(tp/d,fp/d,tn/d,fn/d,d))

    print('MCC: %1.2f' %(MCC(tp,fp,tn,fn)))
