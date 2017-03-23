#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np

if len(sys.argv) < 4:
    raise ValueError('arguments: ./pair_redox_potentials.py sample_redox_potentials_distance.csv sample_redox_potentials.csv sample_dataset.csv')

distances = pd.read_csv(sys.argv[1])
redox = pd.read_csv(sys.argv[2])
output = pd.read_csv(sys.argv[3])

raise Warning("Finish me!") # too tired, leaving this for tomorrow
