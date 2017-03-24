# Reduce Labels

This module is designed to scan the PDB, discover different ways residue groups are labelled and 
map them to a common labelling.

The different labellings and their statistics are stored as a CSV for readability and the same info
is maintained as a pandas.DataFrame in pickle v3 format (a serialized python object).

TODO (performance improvements):
    - turn the pandas access calls into their optimized forms (eg. use iloc instead of brackets)
    - evaluate the use SQLite instead of what will likely become large pickle files
    - write common mapping functions (eg. take in a set of labels and their residue and transform it 
into the most common labelling for consistency across the filtering)
