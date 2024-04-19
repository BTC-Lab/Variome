```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
pd.set_option('display.max_rows', 200)
from matplotlib.patches import ConnectionPatch

```


```python
ploidy_path ='/egp-data/research-g42-khalifa-pfs/data_Request/Merged_files/APfiles/UAE_ploidy_estimation_metrics_merged.csv'
ploidy = pd.read_csv(ploidy_path, sep=",")
ploidy.head(1)

ploidy_estimation = pd.DataFrame(ploidy['Ploidy estimation'].value_counts())
ploidy_estimation
```

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Ploidy estimation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>XY</th>
      <td>24535</td>
    </tr>
    <tr>
      <th>XX</th>
      <td>18982</td>
    </tr>
    <tr>
      <th>XXY</th>
      <td>42</td>
    </tr>
    <tr>
      <th>X0</th>
      <td>23</td>
    </tr>
    <tr>
      <th>XYY</th>
      <td>16</td>
    </tr>
    <tr>
      <th>XXX</th>
      <td>12</td>
    </tr>
  </tbody>
</table>
</div>




```python
ploidy_estimation = pd.DataFrame(ploidy['Ploidy estimation'].value_counts())
ploidy_estimation.loc['Others'] = ploidy_estimation['Ploidy estimation']['XXX'] + ploidy_estimation['Ploidy estimation']['XXY'] + ploidy_estimation['Ploidy estimation']['XYY'] +ploidy_estimation['Ploidy estimation']['X0'] 
ploidy_estimation

```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Ploidy estimation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>XY</th>
      <td>24535</td>
    </tr>
    <tr>
      <th>XX</th>
      <td>18982</td>
    </tr>
    <tr>
      <th>XXY</th>
      <td>42</td>
    </tr>
    <tr>
      <th>X0</th>
      <td>23</td>
    </tr>
    <tr>
      <th>XYY</th>
      <td>16</td>
    </tr>
    <tr>
      <th>XXX</th>
      <td>12</td>
    </tr>
    <tr>
      <th>Others</th>
      <td>93</td>
    </tr>
  </tbody>
</table>
</div>




```python
def getData(texts, autotexts):
    per  = []
    for t in autotexts:
        x = str(t)
        rmSp = x.strip()
        rmSp_spl = rmSp.split("'")
        per.append(rmSp_spl[1])
    data  = [z.strip("%") for z  in per]
    
    label  = []
    for t in texts:
        x = str(t)
        rmSp = x.strip()
        rmSp_spl = rmSp.split(',')
        rmSp_spl[2] = rmSp_spl[2].split("'")[1]
        label.append(rmSp_spl[2])
    label
    perLabel =[]
    numLabel = []
    for i in range(len(label)):
        perLabel.append((label[i] +' \n '+ data[i]))
        #numLabel.append
    
    return (per, data, label, perLabel) 
```


```python
x, data, label, perLabel = getData(texts, autotexts)
x1, data1, label1, perLabel1 = getData(texts2, autotexts2)
```


```python
label
```




    ['XX', 'XY', 'Others']




```python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
length = len(data1)
colors = color_pallete[0:length][0].values.tolist()
colors_bigPie = color_pallete[0:3][0].values.tolist()

# Data for the main categories
categories = perLabel
sizes_main =  data # percentages, should sum to 100
 
# Data for the subcategories of SNPs

# common rare singleton
subcategories_snps = perLabel1
sizes_snps = data1  # percentages, should sum to 100
 
# Create the main pie chart
#fig = plt.figure(figsize=properties['figsize']['pie'])
#fig, ax = plt.subplots(figsize=(10,10))
fig = plt.figure(figsize=properties['figsize']['pie'])
ax = fig.add_subplot(121)
ax_snps = fig.add_subplot(122)


# Explode Shoes slice
#explode = [0.02, 0.02, 0.02]
explode = [0.04, 0.04, 0.04]
ax.pie(sizes_main, autopct='', startangle=0, colors=colors_bigPie, labels=categories,
       pctdistance=1.7, labeldistance=0.6, explode=explode,  radius=5,
       textprops={
        'color': 'white',
        'fontsize': 19,
        'fontweight': 'bold',
        'ha': 'left',
        'va': 'top',
        'rotation': 0
    }
)

ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
 
# Create a sub-pie chart for SNPs
explode = [0.02, 0.02, 0.02, 0.02]
#ax_snps = plt.axes([1,0,1,1])#[0,0,0,0], aspect='equal')#, facecolor='lightgrey')
ax_snps.pie(sizes_snps, labels=subcategories_snps, autopct='', startangle=90, explode=explode,
           colors=colors, labeldistance=.7,
                  textprops={
        'color': 'white',
        'fontsize': 16,
        'fontweight': 'bold',
        'ha': 'center',
        'va': 'center',
        'rotation': 0
    },
    wedgeprops={'linewidth': 0, 'edgecolor': 'black'}
)
 
plt.savefig("../sex_distribution_pie.pdf", bbox_inches='tight')
plt.savefig("../sex_distribution_pie.png", bbox_inches='tight')
plt.show()
```


    
![png](sex_distribution_pie_files/sex_distribution_pie_6_0.png)
    

