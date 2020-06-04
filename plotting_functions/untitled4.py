# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:02:36 2020

@author: mdhurst
"""

import numpy as np
import matplotlib as pyplot

LogL1 = 2980.
LogL2 = 2750.


L1 = np.exp(-LogL1)
L2 = np.exp(-LogL2)

print(L1, L2)