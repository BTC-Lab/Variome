### Steps for extracting annotated variants summary table:
## Step 1:

```python
#!/bin/bash

status="status/status_step1.txt"

# filter files based on repeat maseker regions
for i in {1..22}; do
echo "chr$i start" >> $status
    file_name="maf_VEP_chr${i}.txt" # get file name
    #  varinats in  repeat maseker regions
    bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_In.vcf --targets-file rmsk.bed  fully_annotated/maf_VEP_chr${i}.txt
    
    #  varinats not in  repeat maseker regions
    bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_out.vcf --targets-file ^rmsk.bed  fully_annotated/maf_VEP_chr${i}.txt

    echo "$file_name is done"  >> $status
done 


# the same for chr X,Y,M
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_In.vcf --targets-file rmsk.bed fully_annotated/maf_VEP_chrX.txt # in repeatmasker
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_out.vcf --targets-file ^rmsk.bed  fully_annotated/maf_VEP_chrX.txt # not in repeatmasker

bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_In.vcf --targets-file rmsk.bed  /egp-data/research-g42-khalifa-pfs/maf_g42/fully_annotated/maf_Illumina_VEP_chrY.txt # in repeatmaske
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_out.vcf --targets-file ^rmsk.bed  /egp-data/research-g42-khalifa-pfs/maf_g42/fully_annotated/maf_Illumina_VEP_chrY.txt # in repeatmasker

bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_In.vcf --targets-file rmsk.bed fully_annotated/maf_VEP_chrM.txt # in repeatmasker
bcftools view -O v -o Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_out.vcf --targets-file ^rmsk.bed fully_annotated/maf_VEP_chrM.txt # not in repeatmasker
echo "files are split by in and out of repeat masker regions!"

###################################################################################################################################################################################################################
###################################################################################################################################################################################################################
###################################################################################################################################################################################################################

# check for novel varaints
for i in {1..22}; do
    start_time=$(date +"%H:%M:%S")
    echo "chr$i start time: $start_time" >> "$status"
    file_name="maf_VEP_chr${i}.txt"
    
    #  NOVEL varinats not in  repeat maseker regions
    filter_vep -i Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_out.vcf -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr${i}_not_noRS.vcf --filter "not Existing_variation"

    #  NOVEL varinats in  repeat maseker regions
    filter_vep -i Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_In.vcf -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chr${i}_In_noRS.vcf --filter "not Existing_variation"
    echo "$file_name is done" 
    end_time=$(date +"%H:%M:%S")
    echo "chr$i end time: $end_time" >> "$status"
done 

# the same for chr X,Y,M
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_out.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrM_not_noRS.vcf -filter "not Existing_variation"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrM_In_noRS.vcf -filter "not Existing_variation"
echo "M done!"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_out.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrX_not_noRS.vcf -filter "not Existing_variation"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrX_In_noRS.vcf -filter "not Existing_variation"
echo "X done!"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_out.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrY_not_noRS.vcf -filter "not Existing_variation"
grep -E Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_novel/Table4_repeatmasker_chrY_In_noRS.vcf -filter "not Existing_variation"
echo "Y done!"

echo "files are split by NOVEL in and out of repeat masker regions!"


# check for known varaints
for i in {1..22}; do
    start_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "chr$i Known start time: $start_time" >> "$status"

    file_name="maf_VEP_chr${i}.txt"
    #  KNOWN varinats not in  repeat maseker regions
    filter_vep -i ../Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_out.vcf -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chr${i}_not_RS.vcf --filter "Existing_variation"
    endy_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "chr$i done not in RR: $endy_time" >> "$status"

    #  KNWON varinats in  repeat maseker regions
    filter_vep -i ../Table4_repeatmasker_byChr/Table4_repeatmasker_chr${i}_In.vcf -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chr${i}_In_RS.vcf --filter "Existing_variation"
    echo "$file_name is done" 
    end_time=$(TZ="Asia/Dubai" date +"%H:%M:%S")
    echo "chr$i Known end time: $end_time" >> "$status"
done 

# the same for chr X,Y,M
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_out.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrM_not_RS.vcf -filter "Existing_variation"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrM_In.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrM_In_RS.vcf -filter "Existing_variation"
echo "M done!"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_out.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrX_not_RS.vcf -filter "Existing_variation"../
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrX_In.vcf | filter_vep -o Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrX_In_RS.vcf -filter "Existing_variation"
echo "X done!"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_out.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrY_not_RS.vcf -filter "Existing_variation"
grep -E ../Table4_repeatmasker_byChr/Table4_repeatmasker_chrY_In.vcf | filter_vep -o ../Table4_repeatmasker_byChr/Table4_repeatmasker_Known/Table4_repeatmasker_chrY_In_RS.vcf -filter "Existing_variation"
echo "Y done!"

echo "files are split by KNOWN in and out of repeat masker regions!"

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
