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
sv_path = 'UAE_sv_metrics_merged.csv'
cnv_path = 'UAE_cnv_metrics_merged.txt'
snp_indel_path = '../data/snps_indels.txt'

snp_indel = pd.read_csv(snp_indel_path, sep=",")
sv = pd.read_csv(sv_path, sep=",")
cnv = pd.read_csv(cnv_path, sep="\t")

display(sv.head(1))
display(cnv.head(1))
display(snp_indel.head(1))
```


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


```python
sums_cnv = pd.DataFrame(cnv.sum(axis=0, numeric_only=True)).reset_index()
sums_cnv[0] = sums_cnv[0].astype(int)
sums_cnv
```


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
