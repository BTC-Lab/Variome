```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
```

# Read data


```python
cnv_path = 'UAE_cnv_metrics_merged.txt'

cnv = pd.read_csv(cnv_path, sep="\t")
cnv.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sID</th>
      <th>Ploidy</th>
      <th>Ploidy_confidence</th>
      <th>refBases</th>
      <th>avg_align_cov</th>
      <th>num_align_records</th>
      <th>num_filt_records_all</th>
      <th>num_filt_records_all_ratio</th>
      <th>num_filt_records_dup</th>
      <th>num_filt_records_MAPQ</th>
      <th>...</th>
      <th>num_filt_records_umap_ratio</th>
      <th>cov_uniformity</th>
      <th>num_target_int</th>
      <th>num_segments</th>
      <th>num_amplifications</th>
      <th>num_deletions</th>
      <th>num_pass_amplifications</th>
      <th>num_pass_amplifications_ratio</th>
      <th>num_pas_deletions</th>
      <th>num_pas_deletions_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PGP00019334</td>
      <td>XX</td>
      <td>0.9812</td>
      <td>3105715063</td>
      <td>44.04</td>
      <td>997849408</td>
      <td>68851034</td>
      <td>6.90</td>
      <td>0,0</td>
      <td>59989365</td>
      <td>...</td>
      <td>0.89</td>
      <td>0.14</td>
      <td>2430115</td>
      <td>1935</td>
      <td>163</td>
      <td>566</td>
      <td>88</td>
      <td>53.99</td>
      <td>70</td>
      <td>12.37</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PGP0002329</td>
      <td>XY</td>
      <td>0.9989</td>
      <td>3105715063</td>
      <td>58.82</td>
      <td>1334031125</td>
      <td>94321030</td>
      <td>7.07</td>
      <td>0,0</td>
      <td>83829553</td>
      <td>...</td>
      <td>0.79</td>
      <td>0.16</td>
      <td>2430115</td>
      <td>2719</td>
      <td>185</td>
      <td>531</td>
      <td>96</td>
      <td>51.89</td>
      <td>65</td>
      <td>12.24</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PGP0002582</td>
      <td>XY</td>
      <td>0.9951</td>
      <td>3105715063</td>
      <td>55.94</td>
      <td>1268412193</td>
      <td>89556596</td>
      <td>7.06</td>
      <td>0,0</td>
      <td>79895973</td>
      <td>...</td>
      <td>0.76</td>
      <td>0.15</td>
      <td>2430115</td>
      <td>2239</td>
      <td>172</td>
      <td>512</td>
      <td>98</td>
      <td>56.98</td>
      <td>55</td>
      <td>10.74</td>
    </tr>
  </tbody>
</table>
<p>3 rows × 22 columns</p>
</div>



# Plot


```python
# Default configurations
color_pallete = pd.read_csv('color.txt', header=None)
with open('properties.json', 'r') as file:
    properties=json.load(file)

plt.rcParams['figure.figsize']=properties['figsize']['scatter']
############# CHANGE THIS #############
# specifty which column to plot, to choose the number of colors accordingly
color = color_pallete.loc[0:1][0].values.tolist()
legend = ''
properties['legend'] = False
```


```python
f, ax = plt.subplots()
ax.set_title(properties['title'])
if properties['legend']:
    ax.legend(legend, loc=properties['legend_location'])
    
    
ax.scatter(x = cnv['avg_align_cov'],y = cnv['cov_uniformity'], alpha=0.2, color=color_pallete.loc[0])
ax.set_xlabel('Average alignment coverage')
ax.set_ylabel('Coverage uniformity')
ax.legend()
#plt.savefig("cov_uniformity.png")
plt.show()
```

    No artists with labels found to put in legend.  Note that artists whose label start with an underscore are ignored when legend() is called with no argument.
    


    
![png](cov_uniformity_files/cov_uniformity_5_1.png)
    



```python

```
