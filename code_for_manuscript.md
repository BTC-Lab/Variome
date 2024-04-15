### Variant summary:
1. VCFs were first annotated with VEP
2. Bcftools plugin "split-vep" was used to read the VEP annotations
3. To filter out known variants, filter_vep was used 
- extract novel variants from each VCF - separated by chromosomes


```python
.filter_vep -i input_chr${i}.vcf -o novel_chr${i}.vcf --filter "not Existing_variation" 
```

4. prepare the data to create variant summary tables
- extract the unique variant id grouped by (chrom, pos, ref, alt), the AF, AC, novel variants, snps and indels from each VCF using the following shell script


```python
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    echo $i >> chroms 
    bcftools query -H -f '%CHROM:%POS:%REF:%ALT\t%AF\t%AC\n' VEP_chr${i}.vcf > AF_AC_chr${i}
    bcftools query -H -f '%CHROM:%POS:%REF:%ALT\n' novel_chr${i}.vcf | tail -n +2 > novel_variants_chr${i}
    bcftools view -v snps VEP_chr${i}.vcf | bcftools query -H -f '%CHROM:%POS:%REF:%ALT\t%AF\t%AC\n' > snps_VEP_chr${i}.txt 
    bcftools view -v indels VEP_chr${i}.vcf | bcftools query -H -f '%CHROM:%POS:%REF:%ALT\t%AF\t%AC\n'> indels_VEP_chr${i}.txt 
done
```

5. Create the tables

    a. All Variants


```python
for i in chroms:
    column1.append('chr'+str(i))
    
    df = pd.read_table(f'AF_AC_chr{i}')
    df.columns = ['CHROM:POS:REF:ALT', 'AF', 'AC']

    total = len(df)
    column2.append(total)
    all_total += total
    
    NOVEL = pd.read_table(f'novel_variants_chr{i}', header=None)
    NOVEL.columns = ['CHROM:POS:REF:ALT']

    df_novel = pd.merge(df, NOVEL, how='inner', on=['CHROM:POS:REF:ALT'])
    df_known = df[~df['CHROM:POS:REF:ALT'].isin(df_novel['CHROM:POS:REF:ALT'])]
    novel = len(df_novel)
    known = total - novel
    
    df_known = df_known[df_known['AC'] == 1]
    singleton = len(df_known)
    known = known - singleton

    column3.append(known)
    column4.append((known/total)*100)
    all_known += known
    column5.append(singleton)
    column6.append((singleton/total)*100)
    all_singletonK += singleton

    df_novel = df_novel[df_novel['AC'] == 1]
    singleton = len(df_novel)
    novel = novel - singleton
    column7.append(novel)
    column8.append((novel/total)*100)
    all_novel += novel
    column9.append(singleton)
    column10.append((singleton/total)*100)
    all_singletonN += singleton

    AF = df.loc[df['AF'] > 0.5 ]
    total_AF = len(AF)
    column11.append(total_AF)
    column12.append((total_AF/total)*100)
    all_maa += total_AF

    AF_novel = pd.merge(AF, NOVEL, how='inner', on=['CHROM:POS:REF:ALT'])
    column13.append(len(AF_novel))
    column14.append((len(AF_novel)/total_AF)*100)
    maa_novel += len(AF_novel)

```

    b. SNPs script and the same for INDELS


```python
for i in chroms:
    column1.append('chr'+str(i))
    
    df = pd.read_table(f'snps_VEP_chr{i}')
    df.columns = ['CHROM:POS:REF:ALT', 'AF', 'AC']
    total = len(df)
    column2.append(total)
    all_total += total

    df_common = df.loc[df['AF'] >= 0.05]
    common = len(df_common)
    column3.append(common)
    column4.append((common/total)*100)
    all_common += common

    df_rare = df.loc[(df['AC'] >= 2) & (df['AF'] < 0.05)]
    rare = len(df_rare)
    column5.append(rare)
    column6.append((rare/total)*100)
    all_rare += rare

    df_single = df.loc[df['AC'] == 1]
    singleton = len(df_single)
    column7.append(singleton)
    column8.append((singleton/total)*100)
    all_singleton += singleton

```

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
df = pd.read_csv("num.and.sum.of.regions.csv")
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

### Karyotype Plot:
1. Prepare the ROH data used to plot the ROH density
    
    a. merge the ROH files for each sample to one file "all.roh.bed"
    
    b. remove the centromeres from all.roh.bed
    
    c. extract columns (chromosome, start, end) and save to "roh.bed"
    
2. Extract the frequent and infrequent regions to plot 

    a. Create a bedgraph using the ROH regions for all samples


```python
## add all ROH files for all samples to one file
for sample in sampleList
do 
    cat ${sample}.roh.bed >> all.roh.bed
done

## remove the centromeres from all.roh.bed using bedtools subtract
bedtools subtract -a all.roh.bed -b centromeres.bed > roh.bed

## to get a file with the lengths of the chromosomes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom,size from hg38.chromInfo" | grep -v "_" | grep -v size | grep -v -i chrM 
| sort -Vf > hg38.genome

## be sure that your bed file is sorted based on the chrom and pos
sort -k 1,1 roh.bed > roh.sorted.bed

## produce a bedgraph output by dividing the regions based on the breakpoints in the samples and then counting how many times the region is found in the samples
bedtools genomecov -bga -i roh.sorted.bed -g hg38.genome > roh.merged.bed
```

    b. Get rare and frequent regions


```python
## written in R

filepath <- "roh.merged.bed"
merged <- read.table(filepath, header = TRUE, col.names = c("chrom", "start", "end", "count"))
merged <- merged[!(merged$chrom %in% c("chrX", "chrY")), ]

tmp <- quantile(merged$count, probs = c(0.05,0.95))
rare <- as.numeric(tmp[1]) # anything less than this is rare hom
freq <- as.numeric(tmp[2]) # anything more than this is freq hom

rare_regions <- merged[merged$count <= rare, ]
freq_regions <- merged[merged$count >= freq, ]

write.table(rare_regions, file = "rare_regions.bed", sep = "\t", row.names=FALSE, quote = FALSE)
write.table(freq_regions, file = "freq_regions.bed", sep = "\t", row.names=FALSE, quote = FALSE)
```

3. Extract deletion and duplicatin CNV variants
   - In order to get the regions we create two files ( one for DEL variants and one for DUP )
   - grep DEL first to get all DEL regions and then grep DEL & 1/1 to write the regions again for HOM variants
   - the same goes for DUP
   - the code for this is as follows:



```python
for SM in `cat sample-list`
do 

    grep -i DUP cnv/gff3/${SM}.cnv.gff3 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDUP.bed
    grep -i DUP cnv/gff3/${SM}.cnv.gff3 | grep -i 1/1 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDUP.bed
    
    grep -i DEL cnv/gff3/${SM}.cnv.gff3 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDEL.bed
    grep -i DEL cnv/gff3/${SM}.cnv.gff3 | grep -i 1/1 | sed -e 's/;/\t/g' | sed -e 's/=/\t/g' | cut -f 1,4,5,10,16 >> cnv-bedfiles/annotatedDEL.bed
done  
```


```python
cat cnv-bedfiles/annotatedDEL.bed | cut -f 1,2,3 > cnv-bedfiles/DEL.bed
cat cnv-bedfiles/annotatedDUP.bed | cut -f 1,2,3 > cnv-bedfiles/DUP.bed
```

4. Extract deletion and insertion SV variants
   - In order to get the regions we create two files ( one for DEL variants and one for INS )
   - grep DEL first to get all DEL regions and then grep DEL & 1/1 to write the regions again for HOM variants
   - the same goes for INS
   - the code for this is as follows:


```python
## INS
for SM in `cat sample-list`
do 
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | cut -f 1,2 > sv-bedfiles/INS1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/INS2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/INS3.bed
    paste -d '\t' INS1.bed INS2.bed INS3.bed >> annotatedINS.bed
    
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | grep -i 1/1 | cut -f 1,2 > sv-bedfiles/INS1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | grep -i 1/1 | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/INS2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<INS>" | grep -i 1/1 | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/INS3.bed
    paste -d '\t' INS1.bed INS2.bed INS3.bed >> annotatedINS.bed
done

## DEL
for SM in `cat sample-list`
do 
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | cut -f 1,2 > sv-bedfiles/DEL1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/DEL2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/DEL3.bed
    paste -d '\t' DEL1.bed DEL2.bed DEL3.bed >> annotatedDEL.bed
    
    echo $SM >> annotatedDEL.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | grep -i 1/1 | cut -f 1,2 > sv-bedfiles/DEL1.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | grep -i 1/1 | cut -f 8 | cut -d ';' -f 1 | sed -e 's/END=//g' > sv-bedfiles/DEL2.bed
    zcat sv/${SM}.sv.vcf.gz | grep -v "##" | grep -i "<DEL>" | grep -i 1/1 | cut -f 10 | cut -d ':' -f 1 > sv-bedfiles/DEL3.bed
    paste -d '\t' DEL1.bed DEL2.bed DEL3.bed >> annotatedDEL.bed
done
```


```python
cat sv-bedfiles/annotatedINS.bed | cut -f 1,2,3 > sv-bedfiles/INS.bed
cat sv-bedfiles/annotatedDEL.bed | cut -f 1,2,3 > sv-bedfiles/DEL.bed
```

5. Karyotype script in R


```python
suppressWarnings(library(karyoploteR, quietly = TRUE))

## ROH
filepath <- "roh.bed"
roh <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
roh <- toGRanges(roh)

## rare and freq regions
filepath <- "rare_regions.bed"
rare_regions <- read.table(filepath, header = TRUE)
roh_rare <- toGRanges(rare_regions)
roh_rare <- roh_rare[roh_rare$count > 0,]
roh_rare <- reduce(roh_rare)
roh_rare <- roh_rare[width(roh_rare)> 100000,]

filepath <- "freq_regions.bed" 
freq_regions <- read.table(filepath, header = TRUE)
roh_freq <- toGRanges(freq_regions)
roh_freq <- reduce(roh_freq)
roh_freq <- roh_freq[width(roh_freq)> 100000,]

### PLOT A - Autosomal ROH with rare and freq regions
pdf("../karyotype_plot_a.pdf", 10, 10) #width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = "autosomal")
kpPlotDensity(kp, roh, col = 'gray')
kpPlotRegions(kp, data=roh_rare, data.panel = 2, col = colors[[2]])
kpPlotRegions(kp, data=roh_freq, data.panel = 1, col = colors[[1]] )
legend(x = 0.7, y = 0.7, fill = c ('gray',colors[[2]], colors[[1]]), 
                        legend = c("ROH",'infrequent','frequent'), cex = 1)
dev.off()



## CNV Files
filepath <- "cnv-bedfiles/DEL.bed"
cnv_del <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
cnv_del <- toGRanges(cnv_del)

filepath <- "cnv-bedfiles/DUP.bed"
cnv_dup <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
cnv_dup <- toGRanges(cnv_dup)

### PLOT B - Autosomal CNV
pdf("../karyotype_plot_b.pdf", 10, 10)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = "autosomal")#, plot.params = pp)
kpPlotDensity(kp, cnv_dup, col = colors[[4]])
kpPlotDensity(kp, cnv_del,  data.panel = 2, col = colors[[5]])
legend(x = 0.7, y = 0.7, fill = c (colors[[4]],colors[[5]]), 
                        legend = c("CNV DUP","CNV DEL"), cex = 1)
dev.off()

### PLOT d - Sex CNV
pdf("../karyotype_plot_d.pdf", 7, 7)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = c("chrX", "chrY"))#, plot.params = pp)
kpPlotDensity(kp, cnv_dup, col = colors[[4]])
kpPlotDensity(kp, cnv_del,  data.panel = 2, col = colors[[5]])
legend(x = 0.7, y = 0.3, fill = c (colors[[4]],colors[[5]]), 
                        legend = c("CNV DUP","CNV DEL"), cex = 1)
dev.off()

## SV Files
filepath <- "sv-bedfiles/DEL.bed"
sv_del <- read.table(filepath, header = FALSE, fill = TRUE, col.names = c("chrom", "start", "end"))
sv_del <- toGRanges(sv_del)

## SV Files
filepath <- "INS.bed"
sv_ins <- read.table(filepath, header = FALSE, col.names = c("chrom", "start", "end"))
sv_ins <- toGRanges(sv_ins)

### PLOT C - Autosomal SV
pdf("../karyotype_plot_c.pdf", 10, 10)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = "autosomal")
kpPlotDensity(kp, sv_ins,  col = colors[[1]])
kpPlotDensity(kp, sv_del,  col = colors[[3]], data.panel=2)
legend(x = 0.7, y = 0.7, fill = c (colors[[1]],colors[[3]]), 
                        legend = c("SV INS","SV DEL"), cex = 1)
dev.off()

### PLOT e - Sex SV
pdf("../karyotype_plot_e.pdf", 7, 7)#width=500, height=1000)
kp <- plotKaryotype(plot.type=2, genome = "hg38", chromosomes = c("chrX", "chrY"))
kpPlotDensity(kp, sv_ins,  col = colors[[1]])
kpPlotDensity(kp, sv_del,  col = colors[[3]], data.panel=2)
legend(x = 0.7, y = 0.3, fill = c (colors[[1]],colors[[3]]), 
                        legend = c("SV INS","SV DEL"), cex = 1)
dev.off()
```

### Cluster regions to (Short, Medium, Long):
1. "mclust" package in R was used to categorize the regions lenghts to classess (short, medium, and long)


```python
## import libray
library(mclust)

## read ROH bed for all samples after removing the centromeres
df <- data.table::fread(file = file.path("all.roh.bed"))

## create length column
df <- df %>%
    dplyr::mutate(roh.length = End - Start) %>%
    dplyr::mutate(sampleID = Sample) %>%
    dplyr::left_join(meta, by = "sampleID") 

## cluster to three classes 
dummy <- Mclust(df$roh.length, G = 3)

## add classfication column to original dataframe
df$roh.class <- dummy$classification

## save dataframe to be used for creating plots
cluster.list[["full"]] <- df
saveRDS(cluster.list, file = file.path(file_dir, paste("UAE_roh_mclust",".RDS", sep="")), compress = TRUE)

```


```python

```
