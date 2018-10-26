import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Data
r = [0,1,2,3,4,5,6,7]
mapped = [498817, 708146, 142156, 246812, 256363, 766974, 365267, 58679]
aligned = [303287, 51168, 12662, 21974, 31670, 215491, 118463, 29044]
new_cluster = [88399, 55353, 30972, 19913, 21716, 17535, 16270, 12277]

raw_data = {'greenBars': mapped, 'orangeBars': new_cluster, 'blueBars': aligned}
df = pd.DataFrame(raw_data)
 
# From raw value to percentage
totals = [i+j+k for i,j,k in zip(df['greenBars'], df['orangeBars'], df['blueBars'])]
greenBars = [i / j * 100 for i,j in zip(df['greenBars'], totals)]
orangeBars = [i / j * 100 for i,j in zip(df['orangeBars'], totals)]
blueBars = [i / j * 100 for i,j in zip(df['blueBars'], totals)]
 
# plot
barWidth = 0.85
names = ('ONT','ALZ','RC0','HUM','ZEB', 'SIM-100k', 'SIM-500k', 'SIM-1000k')
# Create green Bars for accepted with mapping
mapped = plt.bar(r, greenBars, edgecolor='black', color='white', width=barWidth, label= "Mapped")
# Create blue Bars for accepted with Alignment
aligned = plt.bar(r, blueBars, bottom=greenBars,  edgecolor='black', color='grey', width=barWidth, label= "Aligned")
# Create orange Bars for new cluster
new_cl = plt.bar(r, orangeBars, bottom=[i+j for i,j in zip(greenBars, blueBars)], edgecolor='black', color='black', width=barWidth, label="New cluster")




# Custom x axis
plt.xticks(r, names,  rotation=90)
plt.xlabel("Dataset")
plt.ylabel("%Read assignment")

# Add a legend
plt.legend(handles=[mapped, aligned, new_cl], loc='upper left', bbox_to_anchor=(1,1), ncol=1)

plt.tight_layout() 
# Show graphic
plt.savefig("/Users/kxs624/Dropbox/qtclust/figures/reads_accepted.pdf")
