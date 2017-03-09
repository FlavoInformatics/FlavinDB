# testing framework
import unittest
# supporting libraries
from biopandas.pdb import PandasPDB
import pandas as pd
import numpy as np
# Classes to be tested

from .FilterInteractions import Filter, _make_distance_comparator
from .FilterInteractions import physical_constants


'''
Testing schedule:
    1. unit test functions individually
    2. test functions with a few macro examples <-- Not sure if worth evaluating yet
    3. ???
    4. $$$ <-- TODO
'''

class TestDistanceFilter(unitest.TestCase):
    def test_distance_comparator(self):
        # simple handpicked smoke test
        keys = ["x_coord", "y_coord", "z_coord"]
        vdWs = {'C': 1.91, 'N': 1.82}
        tollerances = [0, 0.1, 0.15, 0.2, 5] # note that this is just below the default tollerance for distance_comparator
        data = pd.DataFrame()

        # nitrogen center double bonded to carbon with vdW on another nitrogen
        # C = N ... N

        data['atom_name'] = ['C', 'N', 'N']
        data["x_coord"] = [-0.2, 0, 2.1]
        data["y_coord"] = [-0.2, 0, 2.1]
        data["z_coord"] = [-0.2, 0, 2.1]

        correct_output = dict({
            0: [np.nan, np.nan, np.nan],
            0.1: [np.nan, np.nan, np.nan],
            0.15: [np.nan, np.nan, (3 * 2.1 ** 2 ) ** 0.5],
            0.2: [np.nan, np.nan, (3 * 2.1 ** 2 ) ** 0.5],
            5: [(3 * 0.2 ** 2) ** 0.5, np.nan, (3 * 2.1 ** 2 ) ** 0.5],
        })

        for toll in tollerances:
            dist_comp = _make_distance_comparator(data[1], tollerance=toll)
            distances = []
            for i in range(len(data)):
                distances.append(dist_comp(data[i]))

            self.assertEqual(distances, correct_output[toll])


if __name__ == '__main__':
    unittest.main()
