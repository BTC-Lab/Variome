```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
pd.set_option('display.max_rows', 200)
plt.rcParams['figure.figsize']=[20,10]
```


```python
sv_path = '/egp-data/research-g42-khalifa-pfs/data_Request/Merged_files/APfiles/UAE_sv_metrics_merged.csv'
cnv_path = '/egp-data/research-g42-khalifa-pfs/data_Request/Merged_files/APfiles/UAE_cnv_metrics_merged.txt'
snp_indel_path = '../data/snps_indels.txt'

snp_indel = pd.read_csv(snp_indel_path, sep=",")
sv = pd.read_csv(sv_path, sep=",")
cnv = pd.read_csv(cnv_path, sep="\t")

display(sv.head(1))
display(cnv.head(1))
display(snp_indel.head(1))
```


<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>samples</th>
      <th>Number of deletions (PASS)</th>
      <th>Number of insertions (PASS)</th>
      <th>Number of duplications (PASS)</th>
      <th>Number of breakend pairs (PASS)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PGP00019334</td>
      <td>4881</td>
      <td>5548</td>
      <td>54</td>
      <td>1026</td>
    </tr>
  </tbody>
</table>
</div>



<div>

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
      <td>6.9</td>
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
  </tbody>
</table>
<p>1 rows × 22 columns</p>
</div>



<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>snps</th>
      <th>indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>356964849</td>
      <td>64640220</td>
    </tr>
  </tbody>
</table>
</div>



```python
cnv = cnv.drop(columns = ['Ploidy', 'Ploidy_confidence', 'refBases', 'avg_align_cov', 'num_align_records', 'num_filt_records_all',
       'num_filt_records_all_ratio', 'num_filt_records_dup',
       'num_filt_records_MAPQ', 'num_filt_records_MAPQ_ratio',
       'num_filt_records_umap', 'num_filt_records_umap_ratio',
       'cov_uniformity', 'num_target_int', 'num_segments',
       'num_amplifications', 'num_deletions',
       'num_pass_amplifications_ratio',
       'num_pas_deletions_ratio'])
```


```python
sums_sv = pd.DataFrame(sv.sum(axis=0, numeric_only=True)).reset_index()
sums_sv[0] = sums_sv[0].astype(int)
sums_sv
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>0</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Number of deletions (PASS)</td>
      <td>214116191</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Number of insertions (PASS)</td>
      <td>261542968</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Number of duplications (PASS)</td>
      <td>2250102</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Number of breakend pairs (PASS)</td>
      <td>45969355</td>
    </tr>
  </tbody>
</table>
</div>




```python
sums_cnv = pd.DataFrame(cnv.sum(axis=0, numeric_only=True)).reset_index()
sums_cnv[0] = sums_cnv[0].astype(int)
sums_cnv
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>0</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>num_pass_amplifications</td>
      <td>3990455</td>
    </tr>
    <tr>
      <th>1</th>
      <td>num_pas_deletions</td>
      <td>3309613</td>
    </tr>
  </tbody>
</table>
</div>




```python
color_pallete = pd.read_csv('color.txt', header=None)
with open('properties.json', 'r') as file:
    properties=json.load(file)
    
############# CHANGE THIS #############
# specifty which column to plot, to choose the number of colors accordingly
color = color_pallete.loc[0:len(sums_sv)][0].values.tolist()
legend = sums_sv['index'].values.tolist()
```


```python
#plt.figure()
#f, ax = plt.subplots()
plt.figure(figsize=properties['figsize']['pie'])
plt.title(properties['title'])
if properties['legend']:
    plt.legend(legend, loc=properties['legend_location'])

plt.pie(sums_sv[0], labels=sums_sv['index'].values.tolist(), explode=[0.01,0.01,0.1,0.01], colors=color)

plt.savefig('variant_pies_b.png', dpi=300)
plt.show()
```




```python
#plt.figure()
#f, ax = plt.subplots()
plt.figure(figsize=properties['figsize']['pie'])
plt.title(properties['title'])
if properties['legend']:
    plt.legend(legend, loc=properties['legend_location'])

plt.pie(sums_cnv[0], labels=sums_cnv['index'].values.tolist(), explode=[0.01,0.01], colors=color)

plt.savefig('variant_pies_b.png', dpi=300)
plt.show()
```
   



```python

```
