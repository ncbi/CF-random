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

    reference_2qke = '2qke.pdb'
    reference_5jwq = '5jyt.pdb'

    if len(sys.argv) < 2:
        sys.exit('Please enter a name for your output figure')

    ####Uncomment bottom pdb_structures for af-cluster figure###
    pdb_structures = glob('2QKE_329/*/*pdb') #cf-random
    #pdb_structures = glob('00_KaiB_preds_3r/2QKEE_[0-9]??*pdb') #af-cluster


    GS = 0
    FS = 0
    other = 0
    fp = 0
    fn = 0
    tn = 0


    results_dict = {'5jwq rmsd': [], '2qke rmsd': [], 'bfactor': [], 'structure': []}
    rmsd_via_5jwq = []
    bfactor = []
    for structure in pdb_structures:
        structure_index,reference_index = Match_Alpha_Carbons(structure, reference_5jwq)

        results_dict['structure'].append(structure)
        results_dict['5jwq rmsd'].append(Calculate_RMSD(structure, reference_5jwq, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0])
        results_dict['bfactor'].append(np.average(Bfactor(structure)[0:]))

        structure_index,reference_index = Match_Alpha_Carbons(structure, reference_2qke)
        results_dict['2qke rmsd'].append(Calculate_RMSD(structure, reference_2qke, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0])

        if results_dict['2qke rmsd'][-1]*10 <3:
            print(structure,np.average(Bfactor(structure)[0:]),results_dict['2qke rmsd'][-1]*10)

        if results_dict['2qke rmsd'][-1]*10 < 3 and results_dict['bfactor'][-1] >= 70:
            GS += 1
        elif results_dict['5jwq rmsd'][-1]*10 < 3 and results_dict['bfactor'][-1] >= 70:
            FS += 1
        else:
            other += 1

        if results_dict['2qke rmsd'][-1]*10 < 3 and results_dict['bfactor'][-1] < 70:
            fn += 1
        elif results_dict['5jwq rmsd'][-1]*10 < 3 and results_dict['bfactor'][-1] < 70:
            fn += 1

        if results_dict['5jwq rmsd'][-1]*10 > 3 and results_dict['2qke rmsd'][-1]*10 > 3 and results_dict['bfactor'][-1] >= 70:
            fp += 1
        elif results_dict['5jwq rmsd'][-1]*10 > 3 and results_dict['2qke rmsd'][-1]*10 > 3 and results_dict['bfactor'][-1] < 70:
            tn += 1

    df = pd.DataFrame.from_dict(results_dict)
    df['5jwq rmsd']     = df['5jwq rmsd']*10 # nm -> A
    df['2qke rmsd']     = df['2qke rmsd']*10 # nm -> A


    df = df.sort_values('bfactor')

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)
    for idx,_ in enumerate(df['structure']):
        print(idx,_)

    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size
    plt.scatter(data=df, x='2qke rmsd',y='5jwq rmsd', c=df['bfactor'], cmap='rainbow_r',vmin=50,vmax=90)
    plt.xlabel(r'RMSD, Ground state ($\AA$)',fontsize=16)
    plt.ylabel(r'RMSD, FS state ($\AA$)',fontsize=16)
    plt.xlim([0,12])
    plt.ylim([0,12])
    plt.colorbar()
    plt.savefig(sys.argv[1]+'.png',dpi=600)

    s = GS+FS+other

    print('GS: %1.2f, FS: %1.2f, other: %1.2f, all: %i, GS: %i' %(GS/s,FS/s,other/s,s,GS))

    d = fp+tn+fn+GS+FS
    tp = GS+FS

    print('TP: %1.2f, FP: %1.2f, TN: %1.2f, FN: %1.2f, all: %i' %(tp/d,fp/d,tn/d,fn/d,d))

    print('MCC: %1.2f' %(MCC(tp,fp,tn,fn)))
