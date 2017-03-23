#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from biopandas import pdb
import pandas as pd
from physical_constants import get_vdW_radius, vdW_bounds
import random
import sys

"""
    Gets a random sample of 100 protines
    arguments:
        ./get_sample.py <PDB_IDS.csv> <output_data.csv> <sample_size>

    for example:
        ./get_sample.py my_FADs.csv my_FADs_data.csv 100

    :PDB_IDS: A CSV file with the heading "PDB ID", will ignore other columns
        in file.
    :data.csv: Name of the file to write data.
    :sample_size: *Optional* If passed in, script will analyze
        min(sample_size, number of unique PDB_IDs). If not passed in script will
        analyze all PDB_IDs in the file.
"""
if len(sys.argv) < 2:
    raise ValueError("Usage: ./get_sample.py <PDB_IDS.csv: Required> <output_data.csv: Required> <sample_size: Optional>")

PDB_IDS = sys.argv[1]
DATAFILENAME = sys.argv[2]
proteins = []
try:
    proteins = list(pd.read_csv(PDB_IDS)['PDB ID'].unique())
except:
    raise ValueError("Unable to read " + PDB_IDS + " please check that this \
            exists and is in correct format and try again.")

SAMPLE_SIZE = len(proteins)
if len(sys.argv) == 4:
    try:
        SAMPLE_SIZE = int(sys.argv[3])
        print("Using sample size of: " + sys.argv[3])
    except:
        raise ValueError("sample_size must be passed in as a valid integer")

proteins = random.sample(proteins, SAMPLE_SIZE)

def _make_distance_comparator(origin, tolerance=0.2):
    """ Returns a distance functor based upon a given atom

        :origin: one row of a pandas.DataFrame with that contains the columns:
            { "x_coord", "y_coord", "z_coord", "chain_id", "residue_name", "atom_name",
            "residue_number" } in standard PDB format
        :tolerance: an integer error term used when calculating distnace in
            angstroms. If the actual distance of the atoms is within
            <tolerance> angstroms it is said that the two atoms under
            examination are interacting.
    """
    keys = ["x_coord", "y_coord", "z_coord", "chain_id", "residue_name",
            "atom_name", "residue_number"]
    coords = origin[keys].values[0]
    x_o, y_o, z_o = coords[0], coords[1], coords[2]
    residue_o, residue_num_o = coords[4], coords[6]
    vdW_keyatom = get_vdW_radius(origin['atom_name'].values[0], residue_o)
    def distance_comparator(point):
        """
        A distance comparator that returns the distance of the between the input "point"
        atom and the atom based in to origin.

        distance is calculated based on an atoms x, y, z coordinates

        Errors/Exceptions: if the atom's name is not listed in
            "physical_constants.py" then on failing to look it up,
            get_vdW_radius will log the atom's name to standard out and
            distance_comparator will return np.nan, np.nan

        Returns: (distance, label)
        :distance: a double in angstroms describing the 3-D distance between
            the atoms under examination
        :label: an integer value that approximates the electrostatic potentials
            between the atoms under consideration, as retrieved in
            interaction_labels.txt
        If these two atoms do not experience electrostatic interactions or an
            error occured while looking up the atom's van Der Waals' forces then
            the return value is a tuple (np.nan, np.nan)
        """
        coords = point[keys]
        x, y, z = coords[0], coords[1], coords[2]
        residue, name, residue_num = coords[4], coords[5], coords[6]
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

# load residue interaction categories:
residue_categories = pd.read_pickle('./interaction_labels/interaction_dictionary.pkl')

def get_label(atom_name, residue):
    """
        Attempts to look up label for an atom in a given residue. If not found
        logs error to stdout and returns np.nan
    """
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


###############################################################################
###############################################################################
## BEGIN SCRIPT DRIVER ~~
###############################################################################
###############################################################################


# anecdotal names of key atoms in the isoalloxazine to lookup
key_atoms = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C4X', 'N5', 'C5X', 'C6', 'C7',
            'C7M', 'C8', 'C9', 'C9A', 'N10', 'C10']
# PDB names for flavins
flavins = ['FMN', 'FAD']
# taken from biopandas.PDB() object
pdb_columns = ['record_name', 'atom_number', 'blank_1', 'atom_name', 'alt_loc',
  'residue_name', 'blank_2', 'chain_id', 'residue_number', 'insertion',
  'blank_3', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor',
  'blank_4', 'segment_id', 'element_symbol', 'charge', 'line_idx']

# DataFrame that aggregates information accross the columns: contains
#  information on both the target atom and key atom
dataset = pd.DataFrame(columns = ['distance']+ ['key_' + x for x in pdb_columns] +
        ['target_' + x for x in pdb_columns])

for protein in proteins:
    pro = None # running out of names for things at this point
    # attempt to download the protein multiple times from the PDB as this can
    #  fail on occassion
    for _ in range(3):
        try:
            pro = pdb.PandasPDB().fetch_pdb(protein).df
            break
        # not finishing this try will cause compile errors on some implementations
        except:
            continue

    if pro == None:
        # totally failed, log the erorr and move on
        print("UNABLE TO DOWNLOAD: ", protein)
        continue

    # adjust the protein so that caculations
    #  are easier since we make no distinction between atom and heteroatom
    pro = pd.concat([pro['ATOM'], pro["HETATM"]])
    atmnums = [[], []]
    for key in key_atoms:
        key_rows = pro[pro['atom_name'] == key]
        for i in range(len(key_rows)):
            if key_rows.iloc[i]['residue_name'] in flavins:
                atmnums[0].append(key_rows.iloc[i]['atom_number'])
                atmnums[1].append(key)

    # iterate over all of the atoms add it to the DataFrame
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

# Finished computations; log data into provided file
if len(sys.argv[2]):
    try:
        dataset.to_csv(sys.argv[2])
    except:
        # printing could be really painful/useless, but it's still better than losing
        # hours worth of computation
        print(dataset)
else:
    try:
        # try to log the file using a the input file's name, a random nonce and ".csv"
        dataset.to_csv(PDB_IDS + " _raw_dataset." + str(random.randint(0, 1000000)) + ".csv")
    except:
        # otherwise log it based upon its name
        print(dataset)
