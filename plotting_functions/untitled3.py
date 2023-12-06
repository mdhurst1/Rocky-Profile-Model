# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 22:25:43 2023

@author: mh322u
"""

import numpy as np

x = np.linspace(0, 10, 100)
y1 = np.sin(0)

num_segments = 5
x_values = np.linspace(0, 10, len(x))
y2_segments = [np.cos(2 * np.pi * 0.01* x_values) + i for i in range(num_segments)]
segments = [np.column_stack((x_values, y)) for y in y2_segments]

print(segments)