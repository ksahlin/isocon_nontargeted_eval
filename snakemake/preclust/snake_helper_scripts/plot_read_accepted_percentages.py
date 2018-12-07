import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

SMALL_SIZE = 8
MEDIUM_SIZE = 13
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Data
r = [0,1,2,3,4,5,6,7]
mapped = [ 708146, 142156, 246812, 256363, 58679, 365267, 766974, 498817]
aligned = [51168, 12662, 21974, 31670, 29044, 118463, 215491, 303287]
new_cluster = [55353, 30972, 19913, 21716, 12277, 16270, 17535, 88399]

raw_data = {'greenBars': mapped, 'orangeBars': new_cluster, 'blueBars': aligned}
df = pd.DataFrame(raw_data)
 
# From raw value to percentage
totals = [i+j+k for i,j,k in zip(df['greenBars'], df['orangeBars'], df['blueBars'])]
greenBars = [i / j * 100 for i,j in zip(df['greenBars'], totals)]
orangeBars = [i / j * 100 for i,j in zip(df['orangeBars'], totals)]
blueBars = [i / j * 100 for i,j in zip(df['blueBars'], totals)]
 
# plot
barWidth = 0.85
names = ('ALZ','RC0','HUM','ZEB', 'SIM-100k', 'SIM-500k', 'SIM-1000k', 'ONT')
# Create green Bars for accepted with mapping
mapped = plt.bar(r, greenBars, edgecolor='black', color='white', width=barWidth, label= "Minimizer\nmatched")
# Create blue Bars for accepted with Alignment
aligned = plt.bar(r, blueBars, bottom=greenBars,  edgecolor='black', color='grey', width=barWidth, label= "Aligned")
# Create orange Bars for new cluster
new_cl = plt.bar(r, orangeBars, bottom=[i+j for i,j in zip(greenBars, blueBars)], edgecolor='black', color='black', width=barWidth, label="New")




# Custom x axis
plt.xticks(r, names,  rotation=90)
plt.xlabel("Dataset",fontsize=18)
plt.ylabel("%Read assignment",fontsize=18)

# Add a legend
plt.legend(handles=[new_cl, aligned, mapped], loc='upper left', bbox_to_anchor=(1,1), ncol=1)

plt.tight_layout() 
# Show graphic
plt.savefig("/Users/kxs624/Dropbox/qtclust/figures/reads_accepted.pdf")
