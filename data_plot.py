# Importing pandas to export data from text file 
# Importing Matplotlib.pyplot to draw graphs
# Importing numpy to to make faster array like structure to store data
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

# Loading the text file
df_txt = pd.read_csv('g.txt',delimiter='\t')
#
# Printing the data file out to check if data is correctly loaded or not
# This can be later removed
print(df_txt)
