### ROH distribution:
1. Prepare the data used for the ROH and admixture plots

  a. merge the ROH files for each sample to one file "all.roh.bed" 
    
  b. Calculate the total length of ROH and Number of ROH regions for each sample


```python
df = pd.read_table('/path/all.roh.bed')

## df columns are ['Sample', 'Chromosome', 'Start', 'End', 'Score', '#Homozygous', '#Heterozygous']

## remove chrX
df.drop(df[df['Chromosome'] == 'chrX'].index, axis=0, inplace=True)

samples = df['Sample'].unique()

num_regions_per_sample = []
sum_regions_per_sample = []
allsamples = []

#Calculate the total length of ROH and Number of ROH regions for each sample
for sm in samples:
    allsamples.append(sm)
    data = df[df['Sample'] == sm]
    num_regions_per_sample.append(len(data))
    sum_regions_per_sample.append(np.sum(data['End'] - data['Start']))
    
df = pd.DataFrame(columns=['samples', 'number of regions per samples', 'sum of regions per samples'])
df['samples'] = allsamples
df['number of regions per samples'] = num_regions_per_sample
df['sum of regions per samples']    = sum_regions_per_sample
df.to_csv("num.and.sum.of.regions.csv", index=False)
```

  c. The produced file was then merged with admixture data, and consanguinity data
 


```python
df = pd.read_csv(path+"num.and.sum.of.regions.csv")
df.columns = ['sampleID', 'number_of_roh', 'length_of_roh']

df_admix   = pd.read_csv(path+"metadata_full.csv")
df_admix = pd.merge(df, df_admix, how='inner', on=['sampleID'])
df_admix = df_admix[['sampleID', 'number_of_roh', 'length_of_roh', 'main_component']]

df_consang = pd.read_csv(path+"Consang_families.csv")
df_consang = df_consang.loc[df_consang["consang_poID"].str.contains("Child")]
df_consang = pd.merge(df, df_consang, how='inner', on=['sampleID'])
df_consang = df_consang[['sampleID', 'number_of_roh', 'length_of_roh', 'main_component']]

```

  d. Prepare data for plotting 


```python
colors = ["#2CA02C", "#930000", "#1F77B4", "#FFBFD4", "#82CBFF", "#000000"]
populations = pd.unique(df_admix['main_component'].values)

color_mapping = {}
pop_mapping   = {}
components = {}

for p, c, i in zip(populations, colors, range(1,len(colors)+1)):
    color_mapping[p] = c
    pop_mapping[c] = p
    d = 'Component '+str(i)
    components[c] = d

```

#### example code for plotting


```python
# Admixture

df_admix['color'] = df_admix['main_component'].map(color_mapping)

# Create a scatter plot
fig, ax = plt.subplots(figsize=(15,10))

for main_pop in pd.unique(df_admix['main_component'].values):
    x = df_admix[df_admix['main_component'] == main_pop]['length_of_roh']
    y = df_admix[df_admix['main_component'] == main_pop]['number_of_roh']
    color = pd.unique(df_admix[df_admix['main_component'] == main_pop]['color'].values)[0]
    ax.scatter(x, y, c=color, s=5, label=components.get(color), edgecolors='none')
    ax.set_xticklabels([0, 0, 200, 400, 600, 800, 1000])

# Add legend, labels and title
plt.legend(markerscale=6, fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

plt.ylabel('Number of ROH', fontsize=18)
plt.xlabel('Total length in ROH (Mb)', fontsize=18)

# save and show plot
plt.savefig("../ROH_distribution_b.png")
plt.savefig("../ROH_distribution_b.pdf")
plt.show()
```
