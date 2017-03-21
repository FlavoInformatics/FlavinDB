import numpy as np
from biopandas import pdb
import pandas as pd
from physical_constants import get_vdW_radius, vdW_bounds

def _make_distance_comparator(origin, tolerance=0.2):
    keys = ["x_coord", "y_coord", "z_coord", "chain_id", "residue_name", "atom_name"]
    vdW_keyatom = get_vdW_radius(origin['atom_name'].values[0])
    coords = origin[keys].values[0]
    x_o, y_o, z_o = coords[0], coords[1], coords[2]
    chain_o, residue_o = coords[3], coords[4]
    def distance_comparator(point):
        coords = point[keys]
        x, y, z = coords[0], coords[1], coords[2]
        chain, residue, name = coords[3], coords[4], coords[5]
        net_x, net_y, net_z = abs(x - x_o), abs(y - y_o), abs(z - z_o)
        if net_x < vdW_bounds['lower'] and net_y < vdW_bounds['lower'] and net_z < vdW_bounds['lower']:
            distance = np.sqrt(net_x **2 + net_y ** 2 + net_z ** 2)
            expected_weak_iteraction_dist = get_vdW_radius(point['atom_name']) + vdW_keyatom
            if distance < expected_weak_iteraction_dist + tolerance and distance > expected_weak_iteraction_dist - tolerance:
                # if they're in the same flavin, don't return them
                if chain_o == chain and residue_o == residue:
                    return np.nan, np.nan
                # if they aren't return the distance
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
        print("could find", atom_name, "in", residue)
        return -1

atmnums = [[], []]
key_atoms = ['N1', 'C2', 'O2', 'N3', 'C4']
flavins = ['FMN', 'FAD']

dor = pdb.PandasPDB().fetch_pdb('2dor').df
dor = pd.concat([dor['ATOM'], dor["HETATM"]])
for key in key_atoms:
    key_rows = dor[dor['atom_name'] == key]
    for i in range(len(key_rows)):
        if key_rows.iloc[i]['residue_name'] in flavins:
            atmnums[0].append(key_rows.iloc[i]['atom_number'])
            atmnums[1].append(key)

neighbours = []
dataset = pd.DataFrame(columns=['PDB_ID', 'key_atom_name', 'key_atom_number', 'key_atom_residue', 'key_atom_chain_id', 'target_atom_residue', 'target_atom_number', 'target_atom_chain_id', 'distance', 'interaction_label'])

for num in atmnums[0]:
    atom_data = dor[dor.atom_number == num]
    func = _make_distance_comparator(atom_data)
    distance_from_atom_df = dor.apply(func, axis=1)
    dor['distance'] = [x[0] for x in distance_from_atom_df]
    dor['interaction_label'] = [x[1] for x in distance_from_atom_df]
    valid_atom_indices = dor['distance'].notnull()
    valid_key_atoms = dor[valid_atom_indices].sort_values(by=['distance'], ascending=True)
    atom_name = atom_data.atom_name.values[0]
    atom_residue = atom_data.residue_name.values[0]
    chain = atom_data.chain_id.values[0]

    # set temp dataframe and add to dataset
    temp_df = pd.DataFrame(columns=dataset.columns)
    temp_df['distance'] = valid_key_atoms['distance']
    temp_df['PDB_ID'] = ['2dor' for _ in range(len(temp_df))]
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

    neighbours.append(valid_key_atoms)

df = pd.DataFrame()
df['atom_number'] = atmnums[0]
df['neighbors'] = neighbours
df.to_pickle('2dor.pkl')

dataset.to_csv('2dor_complete.csv')
