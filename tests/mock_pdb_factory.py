import pandas as pd
import numpy as np
from biopandas.pdb import PandasPDB
import itertools

'''
    This is a testing framework to generate simple PDB files and objects

    Please note: as of this version - no chemistry has been applied to this
    examples:
        You may find Nitrogen atoms with 6 bonds.
'''
class MockPDB():
    '''
        Design Decisions:
            This module was designed to closely the attributes of a
            biopandas.PandasPDB object.

        Attributes:
            df <dict>: similar to BioPandas, initializing this class will
                return a dictiomary with keys {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}
                that represent the individual parts of a pdb file; where others
                contains all entries that are not parsed as 'ATOM', 'HETATM',
                or 'ANISOU'.

                Upon calling creating an instance of MockPDB, this is what
                will be returned.

        Referrences:
            More information on PDB file format and contents can be found at:
            http://plato.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    '''

    def __init__ (self, random=0, use_real=False, attributes={}):
        '''
            :param random int : if random is set to 0, PDB file stored in sample.
                pdb will be used every time. If random > 0, the number will be
                used as a seed to return a PDB file with random coordinates and
                features.

            :param use_real bool : use_real if set to true will grab a random PDB
                name and download it over the network. Will grab a random protein from
                from a list of proteins in pdb_list.txt (located in this dir)

        '''
        self.df = dict({'ATOM': None, 'HETATM': None, 'ANISOU': None, 'OTHERS': None})
        return
    def _gen_anisou(self):
        ''' Generate information to mimic the ANISOU part of a pdb file.

            Not implemented because not needed as of current version (1/31/2017)


        '''
        self.df = PandasPDB().read_pdb('./test.pdb').df['ANISOU']
        return

    def _gen_others(self):
        ''' Essentialy the 'meta information' at the top of pdb files
            Will add 'OTHER' information to self.df dict


            Returns: None
        '''
        # essentially load the headers from the test.pdb file
        #   which is PDB with some things chagned to make it look obviously wrong
        self.df['OTHERS'] = PandasPDB().read_pdb('./test.pdb').df['OTHERS']
        return

    def _gen_hetatm (self):
        ''' Will build and generate information to mimic HETATM information in a
            normal PDB file

            Format of df['HETATM'] is the same as PandasPDB().df['HETATM']. The columns
            are reporduced below for convenience (examples taken from 2dor.pdb)

                'record_name': 'HETATM'
                'atom_number': integer identifier (see 2dor.pdb for good example)
                'blank_1': empty string for extra information
                'atom_name' multi-character string names of atoms (eg. ['N1', 'C2')
                'alt_loc': 0-len(hetatm) enumeration of the hetatms
                'residue_name': {{ warning: Actual Chemistry }} (eg. 'FMN', 'ORO', 'HOH')
                'blank_2': empty string for extra information
                'chain_id': char denoting the hetatms chain (eg. 'A', 'B')
                'residue_number': {{ warning: Actual Chemistry }} int describing residue group
                'insertion' code insertions for residues
                'blank_3'
                'x_coord': x-coordinate in 3D space
                'y_coord': y-coordinate in 3D space
                'z_coord': z--coordinate in 3D space
                'occupancy': {{ warning: Actual Chemistry }} int
                'b_factor': {{ warning: Actual Chemistry }} int
                'blank_4'
                'segment_id':
                'element_symbol': char of element (eg. 'N', 'C', 'O', 'P')
                'charge': int or np.nan
                'line_idx': int used for enumerating physical structure/etc.

            Returns: None

        '''
        # Caching to save on expense
        if self.df['HETATM'] != None:
            return self.df['HETATM']

        hetatm = pd.DataFrame()
        hetatm['x_coord'] = [0, 1, -1, 2, -2, 4, 2, 1, 2, -2, 3, 1, 2, -1, 3]
        hetatm['y_coord'] = [0, 1, -1, 2, -2, 4, 2, 1, 2, -2, 3, 1, 2, -1, 3]
        hetatm['z_coord'] = [0, 1, -1, 2, -2, 4, 2, 1, 2, -2, 3, 1, 2, -1, 3]
        hetatm['atom_name'] = ['N', 'CA', 'O', 'N', 'CB', 'N', 'CA', 'CB', 'C', 'CB', 'OH', 'CA', 'O', 'N', 'CB']
        hetatm['atom_number'] = [x for x in range(len(hetatm['x_coord']))]
        hetatm['charge'] = [1, np.nan, -1, -1, np.nan, np.nan, -1, np.nan, 1, np.nan, -1, np.nan, np.nan, 1, np.nan]
        self.df['HETATM'] = hetatm
        return hetatm_df

    def _gen_atm (self):
        '''
            See _gen_hetatm for a detailed explanation on these columns
        '''
        atom = pd.DataFrame()
        atom['x_coord'] = [0, 1, -1, 2, -2, 4, 2, 1, 2, -2, 3, 1, 2, -1, 3]
        atom['y_coord'] = [0, 1, -1, 2, -2, 4, 2, 1, 2, -2, 3, 1, 2, -1, 3]
        atom['z_coord'] = [0, 1, -1, 2, -2, 4, 2, 1, 2, -2, 3, 1, 2, -1, 3]
        atom['atom_name'] = ['N', 'CA', 'O', 'N', 'CB', 'N', 'CA', 'CB', 'C', 'CB', 'OH', 'CA', 'O', 'N', 'CB']
        atom['atom_number'] = [x for x in range(len(atom['x_coord']))]
        atom['charge'] = [1, np.nan, -1, -1, np.nan, np.nan, -1, np.nan, 1, np.nan, -1, np.nan, np.nan, 1, np.nan]
        self.df['ATOM'] = atom
        return atom
