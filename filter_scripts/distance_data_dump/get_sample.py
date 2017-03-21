#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from biopandas import pdb
import pandas as pd
from physical_constants import get_vdW_radius, vdW_bounds
import random
import sys

proteins = list(pd.read_csv(sys.argv[1])['PDB ID'].unique())
proteins = random.sample(proteins, 100)
print(proteins)

def _make_distance_comparator(origin, tolerance=0.2):
    keys = ["x_coord", "y_coord", "z_coord", "chain_id", "residue_name", "atom_name", "residue_number"]
    coords = origin[keys].values[0]
    x_o, y_o, z_o = coords[0], coords[1], coords[2]
    chain_o, residue_o, residue_num_o = coords[3], coords[4], coords[6]
    vdW_keyatom = get_vdW_radius(origin['atom_name'].values[0], residue_o)
    def distance_comparator(point):
        coords = point[keys]
        x, y, z = coords[0], coords[1], coords[2]
        chain, residue, name, residue_num = coords[3], coords[4], coords[5], coords[6]
        if residue_num == residue_num_o:
            return np.nan, np.nan
        net_x, net_y, net_z = abs(x - x_o), abs(y - y_o), abs(z - z_o)
        if net_x < vdW_bounds['lower'] and net_y < vdW_bounds['lower'] and net_z < vdW_bounds['lower']:
            distance = np.sqrt(net_x **2 + net_y ** 2 + net_z ** 2)
            expected_weak_iteraction_dist = get_vdW_radius(point['atom_name'], residue) + vdW_keyatom
            if distance < expected_weak_iteraction_dist + tolerance and distance > expected_weak_iteraction_dist - tolerance:
                return distance, get_label(name, residue)
        return np.nan, np.nan

    return distance_comparator

# Residue interaction categories:
residue_categories = pd.read_pickle('./interaction_labels/interaction_dictionary.pkl')

def get_label(atom_name, residue):
    try:
        res = residue_categories[residue]
    except:
        print("residue not found:", residue)
        return -1
    try:
        res = residue_categories[residue]
        return res[res['Residue Atom'] == atom_name].Code.values[0]
    except:
        print("could not find", atom_name, "in", residue)
        return -1

key_atoms = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C4X', 'N5', 'C5X', 'C6', 'C7', 'C7M', 'C8', 'C9', 'C9A', 'N10', 'C10']
flavins = ['FMN', 'FAD']

pdb_columns = ['record_name', 'atom_number', 'blank_1', 'atom_name', 'alt_loc',
  'residue_name', 'blank_2', 'chain_id', 'residue_number', 'insertion',
  'blank_3', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor',
  'blank_4', 'segment_id', 'element_symbol', 'charge', 'line_idx']

dataset = pd.DataFrame(columns = ['distance']+ ['key_' + x for x in pdb_columns] + ['target_' + x for x in pdb_columns])

for protein in proteins:
    pro = None
    try:
        pro = pdb.PandasPDB().fetch_pdb(protein).df
    except:
        try:
            pro = pdb.PandasPDB().fetch_pdb(protein).df
        except:
            try:
                pro = pdb.PandasPDB().fetch_pdb(protein).df
            except:
                print("UNABLE TO DOWNLOAD: ", protein)
                continue
    pro = pd.concat([pro['ATOM'], pro["HETATM"]])
    atmnums = [[], []]
    for key in key_atoms:
        key_rows = pro[pro['atom_name'] == key]
        for i in range(len(key_rows)):
            if key_rows.iloc[i]['residue_name'] in flavins:
                atmnums[0].append(key_rows.iloc[i]['atom_number'])
                atmnums[1].append(key)

    for num in atmnums[0]:
        atom_data = pro[pro.atom_number == num]
        func = _make_distance_comparator(atom_data)
        distance_from_atom_df = pro.apply(func, axis=1)
        pro['distance'] = [x[0] for x in distance_from_atom_df]
        pro['interaction_label'] = [x[1] for x in distance_from_atom_df]
        valid_atom_indices = pro['distance'].notnull()
        valid_key_atoms = pro[valid_atom_indices].sort_values(by=['distance'], ascending=True)
        atom_name = atom_data.atom_name.values[0]
        atom_residue = atom_data.residue_name.values[0]
        chain = atom_data.chain_id.values[0]

        # set temp dataframe and add to dataset
        temp_df = pd.DataFrame(columns=dataset.columns)
        temp_df['distance'] = valid_key_atoms['distance']
        temp_df['PDB_ID'] = [protein for _ in range(len(temp_df))]
        temp_df['key_atom_number'] = [num for _ in range(len(temp_df))]
        temp_df['key_atom_name'] = [atom_name for _ in range(len(temp_df))]
        temp_df['key_atom_residue'] = [atom_residue for _ in range(len(temp_df))]
        temp_df['key_atom_chain_id'] = [chain for _ in range(len(temp_df))]
        temp_df['target_atom_residue'] = valid_key_atoms['residue_name']
        temp_df['target_atom_number'] = valid_key_atoms['atom_number']
        temp_df['target_atom_name'] = valid_key_atoms['atom_name']
        temp_df['target_atom_chain_id'] = valid_key_atoms['chain_id']
        temp_df['interaction_label'] = valid_key_atoms['interaction_label']
        dataset = pd.concat([dataset, temp_df])



try:
    dataset.to_csv(sys.argv[1] + " _raw_dataset." + str(random.randint(0, 1000000)) + ".csv")
except:
    print(dataset)
