from biopandas.pdb import PandasPDB
import pandas as pd
import numpy as np

# set of PDB IDs loaded from a file
_SCANNED_PDB_IDS = None
# Dictionary of Residue that map to sets of their classes
_RESIDUES = None
# Set to true if written to since last write
_RESIDUES_DIRTY = False


# Not "mission critical" but usable for statistics and testing
_PDB_IDS_TO_NOMENCLATURE_CLASSES = None

def filter_PDB_IDS (protein_list):
    global _SCANNED_PDB_IDS
    """ given a list of PDB IDs return a list of only the PDB IDs that haven't
    been scanned yet.

    :protein_list: a iterable object of PDB IDs as a 4 character string
    returns a list of PDB IDs
    """
    if not _SCANNED_PDB_IDS:
        _SCANNED_PDB_IDS = pd.read_pickle("scanned_pdb_ids.pkl")

    return list(filter(lambda x: x in _SCANNED_PDB_IDS, protein_list))

def add_labels (pdb_id):
    global _RESIDUES, _RESIDUES_DIRTY
    if not _RESIDUES:
        _RESIDUES = pd.read_pickle('residue_classes.pkl')
    def get_PDB(attempts = 3):
        """ try a few times to get the PDB file and return a
        PandasPDB object

        fetching the PDB is usually very successful; however, it can fail
        occassionally and we want perfect data here.

        returns biopandas.pdb.PandasPDB
        """
        for _ in range(attempts):
            try: return PandasPDB().fetch_pdb(pdb_id)
            except: continue

    def check_and_update_labels (residue_group, labels):
        global _RESIDUES_DIRTY
        """ given an iterable collection of labels updates _RESIDUES depending on whether or not this
        particular set of labels is already known for this set of elements

        :residue_group: 3 char string of the residue group
        :labels: an (unsorted) iterable of strings
        returns: an :int: id of label if a label was added, else None
        """
        # cast to numpy as array because there are fast vector compare methods
        #  implemented in C there
        labels = np.array(sorted(labels))
        count = 1
        label_id = 0

        # if residue group hasn't been added yet, add first label_set
        if not residue_group in _RESIDUES:
            _RESIDUES[residue_group].update([labels, label_id, count])
            _RESIDUES_DIRTY = True
            return 0
        else:
            # attempt to match sets
            for label_set in _RESIDUES[residue_group]:
                if np.array_equal(labels, label_set[0]):
                    label_set[2] += 1 # update count
                    return None
            # couldn't match set
            label_id = len(_RESIDUES[residue_group])
            _RESIDUES[residue_group].update([labels, label_id, count])
            _RESIDUES_DIRTY = True
            return label_id
    protein = get_PDB()

    # don't stack protein.df['ATOM'] and protein.df['HETATM'] because residue numbers
    #  aren't unique between them
    for atoms in [protein.df['ATOM'], protein.df['HETATM']]:
        # split into sets of labels
        residue_groups = atoms.groupby('residue_number')
        for group_id in residue_groups.groups.keys():
            group = residue_groups.get_group(group_id)
            group_residue_name = group.iloc[0]['residue_name']
            labels = group['atom_name']
            check_and_update_labels(group_residue_name, labels)

def classify_set(pdb_id_list):
    pdb_id_list = filter_PDB_IDS(pdb_id_list)
    for pdb_id in pdb_id_list:
        add_labels(pdb_id)
        _SCANNED_PDB_IDS.update(pdb_id)

    if _RESIDUES_DIRTY:
        pd.to_pickle('residue_classes.py')


# this is that part where a module is also a script
if __name__ == '__main__':
    import sys
    # initialize list
    filename = sys.argv[1] # name of the file to be read
    column_name = sys.argv[2] # name of column with PDB_IDS
    protein_list = pd.read_csv(filename)[column_name].unique()

    # generate labels
    classify_set(protein_list)

