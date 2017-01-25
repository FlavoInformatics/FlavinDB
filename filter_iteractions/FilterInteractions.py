import pandas as pd
import numpy as np
from biopandas.pdb import PandasPDB


def __make_distance_comparator(origin):
    """
        Creates a function to that determines the distance between two heatms in the dataframe
        :param origin pandas.core.series.Series: the origin point, specified as a x, y, z, tuple
    """
    keys = ["x_coord", "y_coord", "z_coord"]
    x_o, y_o, z_o = origin[keys]

    def distance_comparator(point):
        """
            Returns the Euclidean Distance between two points in 3D space
            :param point pandas.core.Series.Series: the point we want to get the distance for
        """
        x, y, z = point[keys]
        return np.sqrt((x - x_o) **2 + (y - y_o) ** 2 + (z - z_o) ** 2)

    return distance_comparator

class Filter():
    """
        This class is designed to take in BioPandas DataFrame and filter
            key hetatms based on iteractions.

        * As of this version Filter() will only ouput atoms based on distance. *

        Based on a ipython notebook written by Akshay Chiwhane (@achiwhane)
        src: https://github.com/antoniomika/FlavinDB/blob/filter-interactions/FilterInteractions.ipynb
    """

    def __init__(self, pandas_pdb=None, filename="", distance_bands=None, key_atoms=[], verbose = False):
        """
            :pandas_pdb: optionally pass in a set of dataframes read in from a PDB file in the BioPandas format
            :filename:

        """
        if pandas_pdb != None:
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

        self.verbose = verbose
        distance_bands = None # load this in from a config file or other options later, but leaving as 'None' for now

    def filter_distance(self):
        """
            filter_distance() requires that the protein is loaded into the object and will return a list of
            HETATMS that are within distance of the key HETRATMS. If key_atoms is empty, then
        """

        if self.protein_dfs == None:
            raise ValueError("No protein loaded into Filter object.")

        hetatms = hetatms = self.protein_dfs.df["HETATM"]
        coords = hetatms[["x_coord", "y_coord", "z_coord"]]
        neighbors = []

        for row_idx in range(len(coords)):
            atom = coords.iloc[row_idx]

            if key_atoms and not (hetatms['atom_number'].iloc[row_idx] in key_atoms):
                continue

            distance_func = __make_distance_comparator(atom)

            # reduce the dataframe of x,y,z to distance from the selected atom
            distance_from_atom_df = hetatms.apply(distance_func, axis=1)

            # sort by distance
            sorted_distance = np.sort(distance_from_atom_df)
            sorted_distance = sorted_distance[1:]

            # get the original indices from the dataframe after sorting to index into hetratm dataframe
            sorted_indices = np.argsort(distance_from_atom_df)
            sorted_indices = sorted_indices[1:]

            # creating a copy so the runtime shuts up and we don't have any issues
            # with accidentally modifying the original hetatm dataframe
            key_hetatm_data = hetatms.iloc[sorted_indices].copy()
            key_hetatm_data['distance'] = sorted_distance

            # drop the col telling us we're dealing with HETATM
            key_hetatm_data.drop("record_name", axis=1, inplace=True)
            # move the atom_name col to the front for writing it to file
            cols = ["atom_name"] + [col for col in key_hetatm_data if col != "atom_name"]
            key_hetatm_data = key_hetatm_data[cols]
            # add this to our results list
            neighbors.append((hetatms["atom_number"].iloc[row_idx], key_hetatm_data))

    def plot(self, x=protein_dfs.df["HETATM"]['x_coord'],
                    y=protein_dfs.df["HETATM"]['y_coord'],
                    z=protein_dfs.df["HETATM"]['z_coord'], label=None):
        # plot 3D representation of protein, because why not
        fig = plt.figure()
        ax = fg.add_subplot(111, projection='3d')
        ax.scatter(x, y, z)
        ax.set_xlabel("x axis")
        ax.set_ylabel("y axis")
        ax.set_zlabel("z axis")
        if label == "":
            ax.set_title("Visualization of 2-DOR")
        else:
            ax.set_title(label)
        plt.tight_layout()
        plt.show()
