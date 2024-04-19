```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import ConnectionPatch
import json

pd.set_option('display.max_rows', 200)
```


```python
ploidy_path ='UAE_ploidy_estimation_metrics_merged.csv'
ploidy = pd.read_csv(ploidy_path, sep=",")
ploidy.head(1)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>samples</th>
      <th>Autosomal median coverage</th>
      <th>X median coverage</th>
      <th>Y median coverage</th>
      <th>1 median / Autosomal median</th>
      <th>2 median / Autosomal median</th>
      <th>3 median / Autosomal median</th>
      <th>4 median / Autosomal median</th>
      <th>5 median / Autosomal median</th>
      <th>6 median / Autosomal median</th>
      <th>...</th>
      <th>16 median / Autosomal median</th>
      <th>17 median / Autosomal median</th>
      <th>18 median / Autosomal median</th>
      <th>19 median / Autosomal median</th>
      <th>20 median / Autosomal median</th>
      <th>21 median / Autosomal median</th>
      <th>22 median / Autosomal median</th>
      <th>X median / Autosomal median</th>
      <th>Y median / Autosomal median</th>
      <th>Ploidy estimation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PGP00019334</td>
      <td>47.66</td>
      <td>46.96</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.99</td>
      <td>0.98</td>
      <td>1.01</td>
      <td>0.97</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.99</td>
      <td>0.99</td>
      <td>0.0</td>
      <td>XX</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 29 columns</p>
</div>




```python
ploidy = ploidy.drop(columns=['samples', 'Autosomal median coverage', 'X median coverage', 'Y median coverage', 'Ploidy estimation'])

ploidy2 = ploidy.rename(columns = {
    '1 median / Autosomal median': 'Chr1', '2 median / Autosomal median': 'Chr2',     '3 median / Autosomal median': 'Chr3', '4 median / Autosomal median': 'Chr4', 
    '5 median / Autosomal median': 'Chr5', '6 median / Autosomal median': 'Chr6',     '7 median / Autosomal median': 'Chr7', '8 median / Autosomal median': 'Chr8', 
    '9 median / Autosomal median': 'Chr9', '10 median / Autosomal median': 'Chr10',   '11 median / Autosomal median': 'Chr11', '12 median / Autosomal median': 'Chr12', 
    '13 median / Autosomal median': 'Chr13', '14 median / Autosomal median': 'Chr14', '15 median / Autosomal median': 'Chr15', '16 median / Autosomal median': 'Chr16', 
    '17 median / Autosomal median': 'Chr17', '18 median / Autosomal median': 'Chr18', '19 median / Autosomal median': 'Chr19', '20 median / Autosomal median': 'Chr20', 
    '21 median / Autosomal median': 'Chr21', '22 median / Autosomal median': 'Chr22', 'X median / Autosomal median': 'ChrX', 'Y median / Autosomal median': 'ChrY', 
})
ploidy2
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr1</th>
      <th>Chr2</th>
      <th>Chr3</th>
      <th>Chr4</th>
      <th>Chr5</th>
      <th>Chr6</th>
      <th>Chr7</th>
      <th>Chr8</th>
      <th>Chr9</th>
      <th>Chr10</th>
      <th>...</th>
      <th>Chr15</th>
      <th>Chr16</th>
      <th>Chr17</th>
      <th>Chr18</th>
      <th>Chr19</th>
      <th>Chr20</th>
      <th>Chr21</th>
      <th>Chr22</th>
      <th>ChrX</th>
      <th>ChrY</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>0.98</td>
      <td>1.01</td>
      <td>0.97</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>0.99</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.99</td>
      <td>0.99</td>
      <td>0.97</td>
      <td>1.00</td>
      <td>0.96</td>
      <td>0.99</td>
      <td>1.00</td>
      <td>0.98</td>
      <td>0.50</td>
      <td>0.51</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.0</td>
      <td>1.01</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>0.98</td>
      <td>1.01</td>
      <td>0.97</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.98</td>
      <td>0.50</td>
      <td>0.51</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>1.01</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.00</td>
      <td>1.02</td>
      <td>0.99</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1.00</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.02</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.0</td>
      <td>1.01</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.99</td>
      <td>0.97</td>
      <td>0.95</td>
      <td>1.01</td>
      <td>0.92</td>
      <td>0.98</td>
      <td>1.00</td>
      <td>0.94</td>
      <td>0.51</td>
      <td>0.52</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>43606</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.01</td>
      <td>0.50</td>
      <td>0.51</td>
    </tr>
    <tr>
      <th>43607</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.98</td>
      <td>1.00</td>
      <td>0.98</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>0.50</td>
      <td>0.51</td>
    </tr>
    <tr>
      <th>43608</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.50</td>
      <td>0.51</td>
    </tr>
    <tr>
      <th>43609</th>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.01</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.00</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>0.98</td>
      <td>1.00</td>
      <td>0.97</td>
      <td>1.00</td>
      <td>1.00</td>
      <td>0.99</td>
      <td>1.00</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>43610</th>
      <td>0.99</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.01</td>
      <td>1.0</td>
      <td>1.01</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.99</td>
      <td>0.98</td>
      <td>0.96</td>
      <td>1.01</td>
      <td>0.93</td>
      <td>0.99</td>
      <td>1.01</td>
      <td>0.96</td>
      <td>0.51</td>
      <td>0.51</td>
    </tr>
  </tbody>
</table>
<p>43611 rows × 24 columns</p>
</div>




```python
# Default configurations
color_pallete = pd.read_csv('color.txt', header=None)
with open('properties.json', 'r') as file:
    properties=json.load(file)

# specifty which column to plot, to choose the number of colors accordingly
colors = color_pallete.loc[0:len(ploidy2.columns)][0].values.tolist()

f, (autosome, sexChr) = plt.subplots(
    ncols=2, 
    nrows=1, 
    sharey=False, 
    figsize=(20,10), 
    gridspec_kw={"width_ratios": (22, 2)}, 
    layout='constrained'
)

#plt.text(11, -0.2, 'Autosomes', ha='center', fontsize=20)
#plt.text(22.5, -0.2, 'Sex chromosome', ha='center', fontsize=10)
sns.boxenplot(data=ploidy2.iloc[:,:22], ax =autosome, palette='Spectral')#,widths=10, points=100, showmeans=True, showxtrema=True, showmedians=True, bw_method=20)#1, native_scale=True)
sns.violinplot(data=ploidy2.iloc[:,22:24], ax =sexChr, palette={'ChrX':'#FFBFD4', 'ChrY': '#1F77B4'})# widths=10)#, native_scale=True)

lsize=15
autosome.tick_params(labelsize=lsize, axis='y')
autosome.tick_params(labelsize=lsize, axis='x',rotation=45 )
sexChr.tick_params(labelsize=lsize, axis='y')
sexChr.tick_params(labelsize=lsize, axis='x',rotation=45 )

plt.show()
#plt.savefig("mean_chr_coverage_ratio.png")
#plt.savefig("mean_chr_coverage_ratio.pdf")
```


```python

```
