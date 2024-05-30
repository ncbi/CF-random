import itertools
import json
import glob
from collections import defaultdict

import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, Selection
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

class CMAP():
    """
    Take in a PDB structure that was loaded using BioPython and produce a contact map
    """

    def __init__(self,structure):
        target_file = structure; pdb = structure.split('/')[-1].replace('.pdb', '')
        structure = PDBParser().get_structure(pdb, target_file)
        self.structure = structure


    def Create_Map(self, cutoff):
        """
        Create the contact map from the PDB structure
        gaps should be filled in and homo-oligomer chain lengths should match
        even if crystal density was not present
        """
        chains = {}
        for chain in self.structure[0].get_chains():
            residues = Selection.unfold_entities(chain, 'R')
            residues = sorted(residues, key=lambda r:r.get_id()[1])
            chains[chain] = [residue['CA'] for residue in residues if 'CA' in residue.child_dict.keys()]

        coor_sep = np.array([[atom.coord for atom in chain[1]] for chain in chains.items()],dtype=object)
        coor_idx = [[i for i,a in enumerate(chain)] for chain in coor_sep]
        tot_idx,offset = [],0
        for idx, coor in enumerate(coor_idx):
            if idx == 0:
                tot_idx.append(list(coor_idx[idx]))
            else:
                offset = offset + coor_idx[idx - 1][-1] +1
                tot_idx.append([pos+offset for pos in coor_idx[idx]])
        tot_idx_flat = [pos for chain in tot_idx for pos in chain]

        coor = np.array([pos for chain in coor_sep for pos in chain])
        dist_matrix = distance_matrix(coor[tot_idx_flat,:],coor[tot_idx_flat,:],p=2)

        #contact map
        contact_map = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        contact_map[dist_matrix < cutoff] = 1
        contact_map[dist_matrix == 0] = 0

        res_list = [res for key in chains for res in chains[key]]
        combo = list(itertools.permutations(res_list,2))
        res_idx = list(itertools.permutations(list(range(len(res_list))),2))

        fill = []
        for key,idx in zip(combo,res_idx):
            key_a,key_b = key[0],key[1]
            if key_a.parent is not None and key_b.parent is not None and idx[1] > idx[0]:
                a = np.array([key_a.parent.child_dict[key].coord for key in key_a.parent.child_dict])
                b = np.array([key_b.parent.child_dict[key].coord for key in key_b.parent.child_dict])
                dist_temp = distance_matrix(a,b,p=2)
                if np.sum(dist_temp <= cutoff).astype(bool):
                    fill.append(idx)
        for i,j in fill:
            contact_map[i][j] = 1

        #last remove hits within +-3 of the diagonal
        mask = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        mask = np.abs(np.arange(len(tot_idx_flat)) - np.arange(len(tot_idx_flat))[:,np.newaxis]) <= 3
        contact_map[mask] = 0

        #upper triangle is a contact map that represents every contact between any heavy atom of two residues
        #lower triangle is a contact map that represents contacts between CAs of two residues
        tu = np.triu_indices(contact_map.shape[0])
        contact_map[tu[::-1]] = contact_map[tu]

        #create a mask that makes oligomers easier to visualize
        #as is this will only handle up to 5 chains, makes interchain squares darker
        mask = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        perm = list(range(len(tot_idx)))
        perm = [p for p in itertools.permutations(perm, r=2) if p[0] < p[1]]

        value = 0.2
        for i in perm:
            mask_idx = np.array(list(itertools.product(tot_idx[i[0]],tot_idx[i[1]])))

            mask[mask_idx[:,0],mask_idx[:,1]] = 1 - int(i[1]-i[0])*value
            mask[mask_idx[:,1],mask_idx[:,0]] = 1 - int(i[1]-i[0])*value

        return contact_map,mask

def pad_with(vector,pad_width,iaxis,kwargs):
    """
    custom paramter to create a padded edge of a given value around an np.array using np.pad()
    """
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

def overlap_definition(A, B, mtx_return=False):
    #pad edges for sliding window consistency
    padded_A = np.pad(A, ((1,1),(1,1)), mode='constant', constant_values=0)
    #Extract all possible windows
    windows = np.lib.stride_tricks.sliding_window_view(padded_A, (3,3))
    #positions with contacts
    mask_1 = np.where(B == 1)
    mask_2 = np.where(B == 2)
    mask_3 = np.where(B == 3)
    #find number of contacts in B that match a (3,3) window in A that has a nonzero elements
    matches_1 = np.any(windows[mask_1] != 0, axis=(-1,-2))
    matches_2 = np.any(windows[mask_2] != 0, axis=(-1,-2))
    matches_3 = np.any(windows[mask_3] != 0, axis=(-1,-2))

    # return B array with +s where a +-1 match in A exists and -s where they don't exist
    # 0s where no contact exists, this function retains the type of contact
    if mtx_return == True:
        B_true = np.copy(B)
        B_true[mask_1] = matches_1 * 1
        B_true[mask_2] = matches_2 * 2
        B_true[mask_3] = matches_3 * 3

        B_false = np.copy(B)
        B_false[mask_1] = ~matches_1 * 1
        B_false[mask_2] = ~matches_2 * 2
        B_false[mask_3] = ~matches_3 * 3
        B_false = B_false*-1
        return B_true + B_false
    return sum([*matches_1, *matches_2, *matches_3])

def COMPARE(_cmap1,_cmap2,pdb1_name,pdb2_name):
    """
    Compare two contact maps and return a contact map alignment based on number of symmetric contacts
    """
    A = _cmap1
    A_name = pdb1_name
    B = _cmap2
    B_name = pdb2_name

    if B.shape[0] < A.shape[0]:
        B = np.pad(B, ((0,A.shape[0]-B.shape[0]),(0,A.shape[0]-B.shape[0])), 'constant')
       
    #pad first matrix to allow for sliding window
    n = int(B.shape[0]*0.1)
    B = np.pad(B, n, pad_with, padder=0)
    mask = np.triu(np.ones(B.shape))
    B = B * mask

    #initialize the maximum number of overlapping 1s
    max_overlap = 0

    #initialize the starting idices of the optimal alignment
    start_i,start_j = 0,0

    beg = range(B.shape[0]-A.shape[0]+1)
    end = range(B.shape[0]-A.shape[0]+1)

    for i,j in zip(beg,end):
        subset = B.copy()
        subset = subset[i:i+A.shape[0], j:j+A.shape[1]]
        # only consider exact matches
        # overlap = len(np.argwhere(np.multiply(A, subset) == 1))

        # consider matches within +- 1
        overlap = overlap_definition(A, subset)
        #update the maximum number of overlapping 1s and the best starting indices if needed
        if overlap > max_overlap:
            max_overlap = overlap
            start_i,start_j = i,j

    lower_triangle_idx = np.tril_indices(A.shape[0],0)

    #offset indices to create the aligned dualfold cmap
    lower_triangle_idx_align = np.transpose(lower_triangle_idx)
    lower_triangle_idx_align = lower_triangle_idx_align + start_i
    lower_triangle_idx_align = np.transpose(lower_triangle_idx_align)
    lower_triangle_idx_align = tuple(np.array(i) for i in lower_triangle_idx_align)

    #create the duafold cmap
    B[lower_triangle_idx_align] = A[lower_triangle_idx]

    return B, B_name, n, A_name, start_i, max_overlap

def QUAD_PLOT(cmap_right, cmap_up, msa_tr, states, state, name):
    """
    Plot comparison of coevolutionary information compared to structures that are predefined to be
    the best representation of each conformation of the fold-switching pair

    Calculate overlap of coevolutionary signal with right, up, and common and return these values
    """

    right_rgba = np.full((*cmap_right.shape, 4), (72/255, 110/255, 158/255, 0), dtype='f')
    right_rgba[np.where(cmap_right == 1)] = (72/255, 110/255, 158/255, 1)

    up_rgba = np.full((*cmap_up.shape, 4), (216/255, 75/255, 89/255 , 0), dtype='f')
    up_rgba[np.where(cmap_up == 1)] = (216/255, 75/255, 89/255, 1)

    tri_upper_idx = np.triu_indices(msa_tr.shape[0])
    tri_lower_idx = np.tril_indices(msa_tr.shape[0])

    right_rgba[np.where((right_rgba[:,:,-1] == 1) & (up_rgba[:,:,-1] == 1))] = (0,0,0,1)
    up_rgba[np.where((right_rgba[:,:,-1] == 1) & (up_rgba[:,:,-1] == 1))] = (0,0,0,1)


    right_rgba[tri_lower_idx[0], tri_lower_idx[1], -1] = right_rgba[tri_lower_idx[0], tri_lower_idx[1], -1] * msa_tr[tri_lower_idx]
    up_rgba[tri_upper_idx[0], tri_upper_idx[1], -1] = up_rgba[tri_upper_idx[0], tri_upper_idx[1], -1] * msa_tr[tri_upper_idx]

    _ = right_rgba[tri_lower_idx[0], tri_lower_idx[1], :]
    right_z_sum = np.sum(_[np.where(((_[:,0] != 0) & (_[:,-1] != 0)))][:, -1])
    _ = up_rgba[tri_upper_idx[0], tri_upper_idx[1], :]
    up_z_sum = np.sum(_[np.where(((_[:,0] != 0) & (_[:,-1] != 0)))][:, -1])
    _ = up_rgba[tri_upper_idx[0], tri_upper_idx[1], :]
    common_z_sum = np.sum(_[np.where(((_[:,0] == 0) & (_[:,-1] != 0)))][:, -1])


    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,8))

    #top left
    ax[0][0].imshow(up_rgba, cmap='gist_yarg', interpolation='none')
    ax[0][0].invert_yaxis()

    #top right
    ax[1][0].imshow(msa_tr, cmap='gist_yarg', interpolation='none')
    ax[1][0].invert_yaxis()

    #bottom left
    labels = states
    sizes = [right_z_sum, up_z_sum]
    labels = list(map(lambda x: "{:}: {:.1f}".format(x[0], x[1]), zip(labels, sizes)))
    ax[0][1].pie(sizes, labels=labels, 
                 colors = [(72/255, 110/255, 158/255, 1),(216/255, 75/255, 89/255, 1)],
                 labeldistance = 0.5)

    #bottom right
    ax[1][1].imshow(right_rgba, cmap='gist_yarg', interpolation='none')
    ax[1][1].invert_yaxis()

    #save quad plot
    fig.suptitle(f'Contat map of predicted structure is closer to {state}')
    plt.savefig(name)
    plt.clf()
    plt.close()

    #save the coevolution in its own plot without the other three panels
    rgba_zscore = np.empty_like(right_rgba)
    rgba_zscore[tri_lower_idx] = right_rgba[tri_lower_idx]
    rgba_zscore[tri_upper_idx] = up_rgba[tri_upper_idx]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
    ax.imshow(rgba_zscore, interpolation='none')
    ax.invert_yaxis()
    plt.savefig(name.replace('.png', '_zscore.png'))
    plt.clf()
    plt.close()
    return common_z_sum, right_z_sum, up_z_sum

def main():
    np.random.seed(0)
    #------------------------------------------------------------------------------------------
    #kaib structures
    # kaib_ground = 'kaib/kaib__1_2/KaiB_unrelaxed_rank_001_alphafold2_model_3_seed_000.r3.pdb'
    # kaib_fs = 'kaib/kaib__16_32/KaiB_unrelaxed_rank_001_alphafold2_model_5_seed_003.r3.pdb'
    # right_structure = kaib_fs; up_structure = kaib_ground; 
    # states = ['fold switch', 'ground']
    # msa_tr_list = glob.glob('kaib/*/*extra_msa_*.npy')
    #__________________________________________________________________________________________

    #------------------------------------------------------------------------------------------
    #mad2 structures
    # mad2_closed = 'mad2/mad2__Full/Mad2_unrelaxed_rank_025_alphafold2_model_4_seed_002.r3.pdb'
    # mad2_open = 'mad2/mad2__2_4/Mad2_unrelaxed_rank_019_alphafold2_model_5_seed_004.r3.pdb'
    # right_structure = mad2_closed; up_structure = mad2_open;
    # states = ['closed', 'open']
    # msa_tr_list = glob.glob('mad2/*/*extra_msa_*.npy')
    #__________________________________________________________________________________________

    #------------------------------------------------------------------------------------------
    #rfah structures
    rfah_beta = 'rfah/rfah__16_32/RfaH_unrelaxed_rank_005_alphafold2_model_1_seed_001.r3.pdb'
    rfah_alpha = 'rfah/rfah__Full/RfaH_unrelaxed_rank_004_alphafold2_model_2_seed_000.r3.pdb'
    right_structure = rfah_beta; up_structure = rfah_alpha;
    states = ['active', 'autoinhibited']
    msa_tr_list = glob.glob('rfah/*/*extra_msa_*.npy')
    #__________________________________________________________________________________________

    #------------------------------------------------------------------------------------------
    #contact map of each conformation of a fold-switching pair
    Cmap_up = CMAP(up_structure) #cmap object now contains structural information!!!
    cmap_up,mask_up = Cmap_up.Create_Map(cutoff=8) # cutoff in angstom
    Cmap_right = CMAP(right_structure) #cmap object now contains structural information!!!
    cmap_right,mask_right = Cmap_right.Create_Map(cutoff=8) # cutoff in angstom
    #__________________________________________________________________________________________

    full_zscore_dict = {}
    #loop through all files containing coevolutionary information from MSATransformer
    for msa_file in msa_tr_list[:]:
        full_zscore_dict[msa_file] = {'common':[], 'right':[], 'up': []}
        np_full = np.load(msa_file)

        #------------------------------------------------------------------------------------------
        #determine which conformation the current runs final pdb is most similar to
        name = msa_file.replace('/rank','/*rank')[:-20] + '.r3.pdb'
        name = glob.glob(name)[0]
        Cmap_current = CMAP(name) #cmap object now contains structural information!!!
        cmap_current,mask_current = Cmap_current.Create_Map(cutoff=8) # cutoff in angstom
        #__________________________________________________________________________________________

        if len(np.where((cmap_current == 1) & (cmap_right == 1))[0]) > len(np.where((cmap_current == 1) & (cmap_up == 1))[0]):
            state = states[0]
        else:
            state = states[1]

        common_z_sum, right_z_sum, up_z_sum = QUAD_PLOT(cmap_right=cmap_right, cmap_up=cmap_up, msa_tr=np_full ,states=states, state=state, name=msa_file.replace('npy', 'png'))

        full_zscore_dict[msa_file]['common'].append(common_z_sum)
        full_zscore_dict[msa_file]['right'].append(right_z_sum)   
        full_zscore_dict[msa_file]['up'].append(up_z_sum)

    #update dictionary for plot
    df_dict = {'name':[],'position':[], 'z_type':[], 'z_score':[]}
    for key, values in full_zscore_dict.items():
        for name, zscores in values.items():
            for z in zscores:
                df_dict['name'].append(key)
                df_dict['position'].append(key.split('/')[1])
                df_dict['z_type'].append(name)
                df_dict['z_score'].append(z)

    df = pd.DataFrame.from_dict(df_dict)
    translate_x_axis = { '1_2':1, '2_4':2, '4_8':3, '8_16':4, '16_32':5, 'full': 7}

    translate_type = {'common': '#000000', 'right': '#486f9e', 'up': '#d84b59'}
    df['x'] = df['position'].map(lambda x: translate_x_axis[x[6:].lower()])

    f, ax = plt.subplots(1,1,figsize=(9,9))
    axes = []
    #each 'type' needs its own handle for matplotlib to give unique legend elements
    for t in df['position'].unique():
        axes.append(ax.scatter(x=df.loc[df['position'].eq(t), 'x'].map(lambda x: x + np.random.uniform(-0.11,0.11)),
                               y= df.loc[df['position'].eq(t), 'z_score'],
                               c=df.loc[df['position'].eq(t), 'z_type'].map(translate_type), 
                               s=20, linewidth=0, linestyle="None"))


    box_dict_up = { 1: [], 2: [], 3: [], 4: [], 5: [],   7: []}
    box_dict_right = { 1: [], 2: [], 3: [], 4: [], 5: [],   7: []}
    box_dict_common = { 1: [], 2: [], 3: [], 4: [], 5: [],   7: []}
    
    for t in df['position'].unique():
        _ = df[df['position'].eq(t)]
    
        df_common_data = _[_['z_type'].eq('common')]
        box_dict_common[df_common_data['x'].unique()[0]].append(df_common_data['z_score'].to_numpy())
    
        df_right_data = _[_['z_type'].eq('right')]
        box_dict_right[df_right_data['x'].unique()[0]].append(df_right_data['z_score'].to_numpy())
    
        df_up_data = _[_['z_type'].eq('up')]
        box_dict_up[df_up_data['x'].unique()[0]].append(df_up_data['z_score'].to_numpy())
    
    np_common = np.zeros((len(list(box_dict_common.keys())) + 1, len(list(box_dict_common.values())[0][0])))
    np_common[:,:] = np.nan
    np_right = np.zeros((len(list(box_dict_right.keys())) + 1, len(list(box_dict_right.values())[0][0])))
    np_right[:,:] = np.nan
    np_up = np.zeros((len(list(box_dict_up.keys())) + 1, len(list(box_dict_up.values())[0][0])))
    np_up[:,:] = np.nan
    
    for key, value in box_dict_common.items():
        np_common[key -1] = value[0]
    axes.append(ax.boxplot(np.transpose(np_common),
                            patch_artist = True,
                            boxprops=dict(facecolor=(0, 0, 0, 0.5), color=(0, 0, 0, 0.5)),
                            capprops=dict(color=(0, 0, 0, 1.0), linewidth = 2),
                            whiskerprops=dict(color=(0, 0, 0, 1.0), linestyle='--' , linewidth=2),
                            flierprops=dict(color=(0, 0, 0, 0.0), markeredgecolor=(0, 0, 0, 0.0)),
                            medianprops=dict(color=(0, 0, 0, 1.0))))
    
    for key, value in box_dict_right.items():
        np_right[key -1] = value[0]
    axes.append(ax.boxplot(np.transpose(np_right),
                            patch_artist = True,
                            boxprops=dict(facecolor=(72/255, 110/255, 158/255, 0.5), color=(72/255, 110/255, 158/255, 0.5)),
                            capprops=dict(color=(0, 0, 0, 1.0), linewidth = 2),
                            whiskerprops=dict(color=(0, 0, 0, 1.0), linestyle='--' , linewidth=2),
                            flierprops=dict(color=(0, 0, 0, 0.0), markeredgecolor=(0, 0, 0, 0.0)),
                            medianprops=dict(color=(0, 0, 0, 1.0))))
    
    for key, value in box_dict_up.items():
        np_up[key -1] = value[0]
    axes.append(ax.boxplot(np.transpose(np_up),
                            patch_artist = True,
                            boxprops=dict(facecolor=(216/255, 75/255, 89/255, 0.5), color=(216/255, 75/255, 89/255, 0.5)),
                            capprops=dict(color=(0, 0, 0, 1.0), linewidth = 2),
                            whiskerprops=dict(color=(0, 0, 0, 1.0), linestyle='--' , linewidth=2),
                            flierprops=dict(color=(0, 0, 0, 0.0), markeredgecolor=(0, 0, 0, 0.0)),
                            medianprops=dict(color=(0, 0, 0, 1.0))))
    
    ax.set_ylabel('Z-score Magnitude')
    ax.set_xlabel('CF-Random Search Depth')
    plt.xticks(sorted(df['x'].unique()), list(translate_x_axis.keys())[:])
    plt.savefig('box_scatter.png')
    plt.clf()
    plt.close()
if __name__ == "__main__":
    main()
