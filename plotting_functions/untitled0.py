# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 21:30:29 2023

@author: mh322u
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# define workspace
Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
ResultsFolder = Folder + "Results/"

DeltaKappaDF = pd.read_excel(Folder+"RSL_Uplift_DeltaKappa.xlsx",sheet_name="Sheet1",header=0)
print(DeltaKappaDF)

# # Example DataFrame
# data = {'A': [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]}
# df = pd.DataFrame(data)

# # Create a heatmap using Seaborn
# plt.figure(figsize=(8, 6))
# sns.heatmap(df, annot=True, cmap='coolwarm', fmt='d', linewidths=.5, cbar_kws={"shrink": 0.8})

# # Set plot labels and title
# plt.xlabel('Columns')
# plt.ylabel('Rows')
# plt.title('Heatmap of DataFrame')

# # Display the plot
# plt.show()