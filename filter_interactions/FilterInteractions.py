import pandas as pd
import numpy as np
from biopandas.pdb import PandasPDB

# module specific packages
# FIXME: generally a bad idea to do relative imports, but doing this for now so PyCharm stops complaining
from physical_constants import vdW_radii, vdW_bounds, get_vdW_radius

def _make_distance_comparator(origin, tolerance=0.2):
    """
        Creates a function to that determines the distance between two heatms in the dataframe
        :param origin pandas.core.series.Series: the origin point, specified as a x, y, z, tuple
        :param tolerance float64: the amount of tollerance (in Angstroms) acceptable in calculating
            whether or not the molecules are within an acceptable distance
    """
    keys = ["x_coord", "y_coord", "z_coord"]
    vdW_keyatom = get_vdW_radius(origin['atom_name']) # FIXME: misclassifies water (ie the O in CO2 is different that O in H2O)
    x_o, y_o, z_o = origin[keys]

    def distance_comparator(point):
        """
            Returns the Euclidean Distance between two points in 3D space
            :param point pandas.core.Series.Series: the point we want to get the distance for
        """
        x, y, z = point[keys]
        net_x, net_y, net_z = abs(x - x_o), abs(y - y_o), abs(z - z_o)

        if net_x < vdW_bounds['lower'] and net_y < vdW_bounds['lower'] and net_z < vdW_bounds['lower']:

            distance = np.sqrt(net_x **2 + net_y ** 2 + net_z ** 2)
            expected_weak_iteraction_dist = get_vdW_radius(point['atom_name']) + vdW_keyatom
            if distance < expected_weak_iteraction_dist + tolerance and distance > expected_weak_iteraction_dist - tolerance:
                return distance

        return np.nan

    return distance_comparator


class Filter():
    """
        This class is designed to take in BioPandas DataFrame and filter
            key hetatms based on iteractions.

        * As of this version Filter() will only ouput atoms based on distance. *

        Based on a ipython notebook written by Akshay Chiwhane (@achiwhane)
        src: https://github.com/antoniomika/FlavinDB/blob/filter-interactions/FilterInteractions.ipynb
    """

    def __init__(self, pandas_pdb=None, filename="", key_atoms=None):
        """
            :paramter pandas_pdb optionally pass in a set of dataframes read in from a PDB file in the BioPandas format
            :filename pass a file on disk or in the pdb to be filtered
            :key_atoms a list of atom numbers to be checked for in the atom
        """
        if pandas_pdb is not None:
            self.protein_dfs = pandas_pdb
        elif filename != "":
            try:
                self.protein_dfs = PandasPDB().read_pdb(filename)
            except FileNotFoundError:
                try:
                    # if it cant be found on disk, try to read it in across the network
                    self.property_dfs = PandasPDB().fetch_pdb(filename)
                except AttributeError:
                    raise ValueError("Filter: PDB file cannot be found on disk or fetch over the network")

        if key_atoms and isinstance(key_atoms, list):
            self.key_atoms = key_atoms # must be a list of atom_numbers
            self.key_atoms.sort()
        else:
            raise ValueError("key_atoms must be passed in as a list() of atom numbers")
        self.neighbors = []


    def filter_distance(self):
        """
            filter_distance() requires that the protein is loaded into the object
            and will return a list of HETATMS that are within distance of the key
            HETRATMS. If key_atoms is empty, then this will assume everything is a
            key atom and return a massive list of neighbors.

            Will cache the results in neighbors and will only return the results
            in neighbors if this function has already been run for this algorithm.
        """
        if len(self.neighbors):
            return self.neighbors

        if not self.protein_dfs:
            raise ValueError("No protein loaded into Filter object.")

        all_atoms = pd.concat([self.protein_dfs.df['ATOM'], self.protein_dfs.df["HETATM"]])

        # list of atom numbers that we want to find neighbors for
        # if no key atoms given, brute force every combination of atoms
        interesting_atoms = self.key_atoms if self.key_atoms else all_atoms.atom_number

        for atom_number in interesting_atoms:
            atom_data = all_atoms[all_atoms.atom_number == atom_number]

            distance_func = _make_distance_comparator(atom_data)

            # reduce the dataframe of x,y,z to distance from the selected atom
            distance_from_atom_df = all_atoms.apply(distance_func, axis=1)
            all_atoms["distance"] = distance_from_atom_df

            # sort relevant atoms based on distance
            valid_atom_indices = all_atoms["distance"].notnull()
            valid_key_atms = all_atoms[valid_atom_indices].sort_values(by=["distance"], ascending=True)

            # we don't want massive dataframes sitting around in memory so just store the indices of the key heteroatoms
            self.neighbors.append((atom_number, valid_key_atms.index.tolist()))

        return self.neighbors
