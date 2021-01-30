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
grid = 1000000

#psi = np.array(df_txt.iloc[0:grid,1])
#x = np.array(df_txt.iloc[0:grid,0])

psi = np.array(df_txt.iloc[42:162,1])
x = np.array(df_txt.iloc[42:162,0])

plt.figure(figsize=(8,5), dpi=100)
plt.plot(x,psi)
plt.xlabel('x-position')
plt.ylabel('Wave Function')
plt.grid()
plt.title("Second Excited State for Harmonic Oscillator")
plt.show()
plt.savefig('second_excited_State.png', dpi=300)