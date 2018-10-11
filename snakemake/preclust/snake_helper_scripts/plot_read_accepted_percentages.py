import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
 
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
names = ('BHAM_ONT','ALZ_PB','RC0_PB','HUM_PB','ZEB_PB', 'ENS_100k_PB', 'ENS_500k_PB', 'ENS_1M_PB')
# Create green Bars for accepted with mapping
mapped = plt.bar(r, greenBars, color='#b5ffb9', edgecolor='white', width=barWidth, label= "Mapped")
# Create blue Bars for accepted with Alignment
aligned = plt.bar(r, blueBars, bottom=greenBars, color='#a3acff', edgecolor='white', width=barWidth, label= "Aligned")
# Create orange Bars for new cluster
new_cl = plt.bar(r, orangeBars, bottom=[i+j for i,j in zip(greenBars, blueBars)], color='#f9bc86', edgecolor='white', width=barWidth, label="New cluster")
 
# Custom x axis
plt.xticks(r, names,  rotation=90)
plt.xlabel("Datset")
plt.ylabel("%Read assignment")

# Add a legend
plt.legend(handles=[mapped, aligned, new_cl], loc='upper left', bbox_to_anchor=(1,1), ncol=1)

plt.tight_layout() 
# Show graphic
plt.savefig("/Users/kxs624/tmp/QUTE_CLUST_PAPER/reads_accepted.pdf")
