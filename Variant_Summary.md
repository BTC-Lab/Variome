### Variant summary:
1. VCFs were first annotated with VEP
2. Bcftools plugin "split-vep" was used to read the VEP annotations
3. To filter out known variants, filter_vep was used 
- extract novel variants from each VCF - separated by chromosomes - {i} referes to chromosome name


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
