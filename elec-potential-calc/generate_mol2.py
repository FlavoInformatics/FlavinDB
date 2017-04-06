import os

from chimera import Molecule
from chimera import openModels
from chimera import selection
from chimera import runCommand
from chimera.fetch import fetchPDB
from chimera.selection import currentAtoms


# TODO: wrap this in an easy-to-use-cli

pdb_code = "1ahv"
 
template = openModels.open(pdb_code, type="PDB")

models = openModels.list(modelTypes=[Molecule])

flavo_ligands = []
for m in models:
    for r in template[0].residues:
        if r.type=="FAD" or r.type=="FMN":
            # set the current selection to this residue
            selection.setCurrent(r)
            runCommand("write format mol2 selected 0 ~/Desktop/{}.{}.mol2".format(pdb_code, r.id))




