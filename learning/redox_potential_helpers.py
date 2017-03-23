#import pandas as pd
#import numpy as np
from biopandas.pdb import PandasPDB

def check_structure_exists(name):
    if not name:
        raise ValueError("Empty name, cannot check if structure is valid")
    pro = None
    for _ in range(3):
        try:
            pro = PandasPDB().fetch_pdb(name)
            if pro: break
        except:
            continue
    if not hasattr(pro, 'df'):
        return False
    if len(pro.df['HETATM']) == 0 or len(pro.df["ATOM"]) == 0:
        return False

    # no errors; this structure is probably fine?
    return True

def trim_list(names, category):
    """ check all of the names in the list :names: to check that they have a PDB structure """
    not_found = []
    found = []
    for name in names:
        try:
            if not check_structure_exists(name):
                not_found.append(name)
            else:
                found.append(name)
        except:
            not_found.append(name)

    return (found, not_found)
