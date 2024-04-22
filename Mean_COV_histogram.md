```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
pd.set_option('display.max_rows', 200)
```


```python
wgsCov_path = 'UAE_wgs_mean_cov_merged.txt'
```


```python
wgsCov = pd.read_csv(wgsCov_path, sep="\t")
wgsCov['mean'] = wgsCov['Content'].str.split(',', expand=True)[1]
wgsCov['Content'] = wgsCov['Content'].str.split(',', expand=True)[0]
wgsCov['mean'] = wgsCov['mean'].astype(float)
wgsCov.head(1)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Content</th>
      <th>Filename</th>
      <th>mean</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Average alignment coverage over wgs</td>
      <td>sample1.wgs_overall_mean_cov.csv</td>
      <td>42.78</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Default configurations
color_pallete = pd.read_csv('color.txt', header=None)
with open('properties.json', 'r') as file:
    properties=json.load(file)
    
############# CHANGE THIS #############
# specifty which column to plot, to choose the number of colors accordingly
colors = color_pallete[0:1][0].values.tolist()

plt.rcParams['figure.figsize']=properties['figsize']['scatter']
legend = ''
properties['legend'] = False

############# PLOT #############


f, ax = plt.subplots()
sns.histplot(wgsCov['mean'] , ax= ax, color= color_pallete.iloc[1])
ax.set(xlabel='Mean Coverage')
ax.set(ylabel='Frequency')
#ax.legend('Mean Coverage')
#plt.savefig('mean_cov_histogram.png', dpi=300)
```
