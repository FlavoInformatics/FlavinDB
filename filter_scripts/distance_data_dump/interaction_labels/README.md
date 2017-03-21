# Interaction Labels
## categorical labels for interactions in protein atom, chemical codes

---

### Description:
These are CSV files contiaining (Header: "Residue Atom, Code") under the respective amino acid code (eg. ALA) describing the interaction that each atom should be labelled as

---

### Usage:
You can load them one by one if you choose or you can just load the pkl file, which is a dictionary of 3 digit AA to its dictionary of Residue Atom to Code

---

### Example
TBD: example

---

### Updating:
it's essentially the following code (ipython3+)

'''
import pandas as pd
ls = !ls
files = [f for f in ls if len(f) == 3] # hacky, I know
dicts = dict()
for file in files:
    dicts[file] = pd.read_csv(file)

pd.to_pickle(obj=dicts, path='./interaction_dictionary.pkl')
'''

-- WS (@Ramanujamin), March 5 2016
