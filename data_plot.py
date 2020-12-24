# Importing pandas to export data from text file 
# Importing Matplotlib.pyplot to draw graphs
# Importing numpy to to make faster array like structure to store data
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

# Loading the text file
df_txt = pd.read_csv('output.txt',delimiter='\s+')
#
# Printing the data file out to check if data is correctly loaded or not
# This can be later removed
# print(df_txt)
grid = 100000

f_pot = np.array(df_txt.iloc[0:grid,1])
x = np.array(df_txt.iloc[0:grid,0])

plt.figure(figsize=(8,5), dpi=100)
plt.plot(x,f_pot)
plt.xlabel('x-position')
plt.ylabel('Wave Function')
plt.show()
plt.savefig('psi_2.png', dpi=300)