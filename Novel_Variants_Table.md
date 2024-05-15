```python
import os
import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 200)
```


```python
# filter annotated MAF table
!grep -v -E ";AF=0\.0[0-4]" maf_Illumina_VEP.txt | filter_vep -o maf_Illumina_VEP_filtered.txt -filter "(IMPACT in HIGH,MODERATE) and (gnomADe_AF < 0.01 or not gnomADe_AF)"

!bcftools view -O v -o noRepeatmasker.vcf --targets-file ^rmsk.bed  maf_Illumina_VEP_filtered.txt # filter out repeatmasker regions
!bcftools view -O v -o noRepeatmasker_noLCR.vcf --targets-file ^LCRFromHengHg38.bed noRepeatmasker.vcf # filter out low-complexity regions
!bcftools view -O v -o noRepeatmasker_noLCR_noSD.vcf --targets-file ^hglft_genome_eec9_942d10.bed noRepeatmasker_noLCR.vcf # filter out segmental duplication
!bcftools view -O v -o noRepeatmasker_noLCR_noSD_protien.vcf --targets-file prot_coding_variants.bed noRepeatmasker_noLCR_noSD.vcf # include protien coding variants

!filter_vep -i noRepeatmasker_noLCR_noSD_protien.vcf -o Table5_filtered.vcf -filter "(MAX_AF < 0.01 or not MAX_AF)" # filter out all variants that have AF > 0.01 in other populations
!grep -v "^#" Table5_filtered.vcf > Table5_filtered.txt #convert vcf to text file and remove the header
```

# Important message


```python
cols_to_drop = []
header = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS".split("|")
```


```python
# Read MGI and illumina files
illumina_variants = pd.read_csv("Table5_filtered.txt", sep="\t", header=None) #filtered illumina file
len(illumina_variants)
```




    74




```python
illumina_variants.head(1)
```





<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>7829912</td>
      <td>.</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>.</td>
      <td>PASS</td>
      <td>AC=31878;AF=0.365414154382265;RC=55360;RF=0.63...</td>
    </tr>
  </tbody>
</table>



# dataframe with the important column
In the VCF file some variants where associated with multiple genes in (((1 line))). 
these are split, in which each line is a row, then i removed the duplicate rows.

now we can see in the index the same number means they are of the same varinat.


```python
INFO_illumina = illumina_variants[7].str.split(";", expand=True)
INFO_illumina.head(3)
```





<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AC=31878</td>
      <td>AF=0.365414154382265</td>
      <td>RC=55360</td>
      <td>RF=0.634585845617735</td>
      <td>CSQ=-|inframe_deletion|MODERATE|PER3|ENSG00000...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AC=6117</td>
      <td>AF=0.0701185263302689</td>
      <td>RC=81121</td>
      <td>RF=0.929881473669731</td>
      <td>CSQ=T|missense_variant|MODERATE|PRAMEF6|ENSG00...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AC=8247</td>
      <td>AF=0.0945344918498819</td>
      <td>RC=78991</td>
      <td>RF=0.905465508150118</td>
      <td>CSQ=T|missense_variant|MODERATE|PRAMEF6|ENSG00...</td>
    </tr>
  </tbody>
</table>




```python
print(INFO_illumina.loc[0])
for i in (INFO_illumina[4].str.split(',')[0]):
    print('\n', i)
```

    0                                             AC=31878
    1                                 AF=0.365414154382265
    2                                             RC=55360
    3                                 RF=0.634585845617735
    4    CSQ=-|inframe_deletion|MODERATE|PER3|ENSG00000...
    Name: 0, dtype: object
    
     CSQ=-|inframe_deletion|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000361923|protein_coding|18/21||ENST00000361923.2:c.3019_3072del|ENSP00000355031.2:p.Ala1007_Ser1024del|3114-3167|2939-2992|980-998|ENPSHPTASALSTGSPPMK/E|gAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGAag/gag|||1||deletion|HGNC|HGNC:8847||||1|P4|CCDS89.1|ENSP00000355031|P56645.187|A2I2N5.99|UPI0000167B1D|P56645-1|1|||AlphaFold_DB_import:AF-P56645-F1&PANTHER:PTHR11269&MobiDB_lite:mobidb-lite||80|||||||||||||||||||||||||||||||||||||
    
     -|inframe_deletion|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000377532|protein_coding|19/22||ENST00000377532.8:c.3046_3099del|ENSP00000366755.3:p.Ala1016_Ser1033del|3276-3329|2966-3019|989-1007|ENPSHPTASALSTGSPPMK/E|gAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGAag/gag|||1||deletion|HGNC|HGNC:8847|YES|NM_001377275.1||1|A2|CCDS72695.1|ENSP00000366755|P56645.187||UPI00003664CA|P56645-2|1|||AlphaFold_DB_import:AF-P56645-F1&PANTHER:PTHR11269&MobiDB_lite:mobidb-lite&MobiDB_lite:mobidb-lite||80|||||||||||||||||||||||||||||||||||||
    
     -|upstream_gene_variant|MODIFIER||ENSG00000236266|Transcript|ENST00000451646|lncRNA|||||||||||2571|-1||deletion|||YES|||3||||||||||||||||||||||||||||||||||||||||||||||||||
    
     -|inframe_deletion|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000613533|protein_coding|19/22||ENST00000613533.4:c.3046_3099del|ENSP00000482093.1:p.Ala1016_Ser1033del|3230-3283|2966-3019|989-1007|ENPSHPTASALSTGSPPMK/E|gAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGAag/gag|||1||deletion|HGNC|HGNC:8847||||5|A2|CCDS72695.1|ENSP00000482093|P56645.187||UPI00003664CA|P56645-2|1|||AlphaFold_DB_import:AF-P56645-F1&PANTHER:PTHR11269&MobiDB_lite:mobidb-lite&MobiDB_lite:mobidb-lite||80|||||||||||||||||||||||||||||||||||||
    
     -|splice_donor_variant&splice_donor_5th_base_variant&coding_sequence_variant&intron_variant|HIGH|PER3|ENSG00000049246|Transcript|ENST00000614998|protein_coding|19/23|19/22|ENST00000614998.4:c.3002-13_3042del||3227-?|2963-?|988-?|||||1||deletion|HGNC|HGNC:8847||||1|A2|CCDS76097.1|ENSP00000479223||A0A087WV69.52|UPI000387D6D4||1|||||80|||||||||||||||||||||||||||||||||||||
    


```python
#Illumina files
# new dataframe only for the INFO column
#INFO_illumina = illumina_variants[7].str.split(";", expand=True)
INFO_illumina[4] = INFO_illumina[4].str.split(',')
INFO_illumina = INFO_illumina.explode(4)
INFO_illumina.insert(0, 'alt', illumina_variants[4])
INFO_illumina.insert(0, 'ref', illumina_variants[3])
INFO_illumina.insert(0, 'pos', illumina_variants[1])
INFO_illumina.insert(0, 'chrom', illumina_variants[0])
# illumina add more columns the annotations dataframe
gnomAD_illumina = INFO_illumina[4].str.split('|', expand=True)
gnomAD_illumina.columns = header
gnomAD_illumina = gnomAD_illumina.drop(columns = cols_to_drop)
gnomAD_illumina.insert(0, 'UAE_AF', INFO_illumina[1])
gnomAD_illumina.insert(0, 'alt', illumina_variants[4])
gnomAD_illumina.insert(0, 'ref', illumina_variants[3])
gnomAD_illumina.insert(0, 'pos', illumina_variants[1])
gnomAD_illumina.insert(0, 'chr', illumina_variants[0])
gnomAD_illumina = gnomAD_illumina.drop_duplicates()
```


```python
gnomAD_illumina.head(5)
```





<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chr</th>
      <th>pos</th>
      <th>ref</th>
      <th>alt</th>
      <th>UAE_AF</th>
      <th>Allele</th>
      <th>Consequence</th>
      <th>IMPACT</th>
      <th>SYMBOL</th>
      <th>Gene</th>
      <th>...</th>
      <th>MAX_AF_POPS</th>
      <th>CLIN_SIG</th>
      <th>SOMATIC</th>
      <th>PHENO</th>
      <th>PUBMED</th>
      <th>MOTIF_NAME</th>
      <th>MOTIF_POS</th>
      <th>HIGH_INF_POS</th>
      <th>MOTIF_SCORE_CHANGE</th>
      <th>TRANSCRIPTION_FACTORS</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>7829912</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>AF=0.365414154382265</td>
      <td>CSQ=-</td>
      <td>inframe_deletion</td>
      <td>MODERATE</td>
      <td>PER3</td>
      <td>ENSG00000049246</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>7829912</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>AF=0.365414154382265</td>
      <td>-</td>
      <td>inframe_deletion</td>
      <td>MODERATE</td>
      <td>PER3</td>
      <td>ENSG00000049246</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>7829912</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>AF=0.365414154382265</td>
      <td>-</td>
      <td>upstream_gene_variant</td>
      <td>MODIFIER</td>
      <td></td>
      <td>ENSG00000236266</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>7829912</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>AF=0.365414154382265</td>
      <td>-</td>
      <td>inframe_deletion</td>
      <td>MODERATE</td>
      <td>PER3</td>
      <td>ENSG00000049246</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>7829912</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>AF=0.365414154382265</td>
      <td>-</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>HIGH</td>
      <td>PER3</td>
      <td>ENSG00000049246</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
<p>5 rows × 83 columns</p>




```python
# genes in illumina
gnomAD_illumina["SYMBOL"].value_counts()
```




    JAKMIP3      46
                 46
    GSTT4        44
    PEX5         20
    NPRL3        20
    METTL22      13
    CPSF1        10
    MUC3A         9
    ZAN           8
    GOLGA6L2      8
    NBPF20        8
    DDT           8
    KRT10-AS1     7
    POM121C       7
    OPN3          7
    PDIA2         7
    ADCK5         7
    NBPF10        6
    PRDM9         6
    EPPK1         6
    ANXA8L1       6
    STOX1         5
    TXNDC2        5
    PRSS2         5
    SLC22A1       5
    ARHGDIG       4
    POTEB         4
    AXIN1         4
    PER3          4
    PRAMEF6       4
    MUC5AC        4
    POTEM         3
    GOLGA6L10     3
    TBC1D3I       3
    POTEG         3
    PON1          3
    PRAMEF13      2
    PCARE         2
    NPY4R         2
    RP1L1         2
    KRT10         2
    HBZ           2
    PLIN4         2
    FADS6         2
    GSTT3P        1
    PPIAL4F       1
    DDTL          1
    GOLGA6L1      1
    KRTAP4-3      1
    PRAMEF5       1
    GOLGA6L22     1
    MIR939        1
    HECA          1
    Name: SYMBOL, dtype: int64




```python
gnomAD_illumina.value_counts('IMPACT')
```




    IMPACT
    MODIFIER    201
    MODERATE    124
    HIGH         58
    dtype: int64



# Filter specific columns and merge illumina with mgi
 
### in this section only specific column names are chosen, and then the mgi and illumina tables are merged


```python
filter_cols=["chr","pos","ref","alt", "SYMBOL", "UAE_AF", "MAX_AF", "MAX_AF_POPS", "Existing_variation", "IMPACT", "Consequence"]
```


```python
#select sepcific columns only - FOR COMPARING
illumina_filter_cols = gnomAD_illumina#[filter_cols]
illumina_filter_cols = illumina_filter_cols.drop_duplicates()
illumina_filter_cols = illumina_filter_cols.sort_values(by = 'SYMBOL', ascending=False)
illumina_filter_cols["UAE_AF"] = illumina_filter_cols["UAE_AF"].str.split('=', expand=True)[1]
illumina_filter_cols["UAE_AF"] = pd.to_numeric(illumina_filter_cols["UAE_AF"])
illumina_filter_cols = illumina_filter_cols.reset_index().drop(columns = "index")
illumina_filter_cols['chr'] = illumina_filter_cols['chr'].astype(str)
```


```python
illumina_filter_cols.head(5)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chr</th>
      <th>pos</th>
      <th>ref</th>
      <th>alt</th>
      <th>UAE_AF</th>
      <th>Allele</th>
      <th>Consequence</th>
      <th>IMPACT</th>
      <th>SYMBOL</th>
      <th>Gene</th>
      <th>...</th>
      <th>MAX_AF_POPS</th>
      <th>CLIN_SIG</th>
      <th>SOMATIC</th>
      <th>PHENO</th>
      <th>PUBMED</th>
      <th>MOTIF_NAME</th>
      <th>MOTIF_POS</th>
      <th>HIGH_INF_POS</th>
      <th>MOTIF_SCORE_CHANGE</th>
      <th>TRANSCRIPTION_FACTORS</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7</td>
      <td>100773854</td>
      <td>G</td>
      <td>GTT</td>
      <td>0.499702</td>
      <td>TT</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>MODIFIER</td>
      <td>ZAN</td>
      <td>ENSG00000146839</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>1</th>
      <td>7</td>
      <td>100773854</td>
      <td>G</td>
      <td>GTT</td>
      <td>0.499702</td>
      <td>TT</td>
      <td>frameshift_variant</td>
      <td>HIGH</td>
      <td>ZAN</td>
      <td>ENSG00000146839</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>2</th>
      <td>7</td>
      <td>100773854</td>
      <td>G</td>
      <td>GTT</td>
      <td>0.499702</td>
      <td>TT</td>
      <td>frameshift_variant</td>
      <td>HIGH</td>
      <td>ZAN</td>
      <td>ENSG00000146839</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>3</th>
      <td>7</td>
      <td>100773854</td>
      <td>G</td>
      <td>GTT</td>
      <td>0.499702</td>
      <td>TT</td>
      <td>frameshift_variant</td>
      <td>HIGH</td>
      <td>ZAN</td>
      <td>ENSG00000146839</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>4</th>
      <td>7</td>
      <td>100773854</td>
      <td>G</td>
      <td>GTT</td>
      <td>0.499702</td>
      <td>TT</td>
      <td>frameshift_variant</td>
      <td>HIGH</td>
      <td>ZAN</td>
      <td>ENSG00000146839</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
<p>5 rows × 83 columns</p>
</div>




```python
merge = pd.merge(mgi_filter_cols,illumina_filter_cols, how='inner', on =['chr', 'ref', 'alt', 'pos'])
merge.rename(columns = {'chr_x': 'chr' , 'UAE_AF_x': 'UAE_AF_mgi', 'UAE_AF_y': 'UAE_AF_illumina', 'Existing_variation_x': 'Existing_variation_mgi', 'Existing_variation_y': 'Existing_variation_illumina'}, inplace=True)
merge = merge.sort_values(by = 'UAE_AF_illumina', ascending=False)
merge.head(1)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chr</th>
      <th>pos</th>
      <th>ref</th>
      <th>alt</th>
      <th>UAE_AF_mgi</th>
      <th>UAE_AF_illumina</th>
      <th>Allele</th>
      <th>Consequence</th>
      <th>IMPACT</th>
      <th>SYMBOL</th>
      <th>...</th>
      <th>MAX_AF_POPS</th>
      <th>CLIN_SIG</th>
      <th>SOMATIC</th>
      <th>PHENO</th>
      <th>PUBMED</th>
      <th>MOTIF_NAME</th>
      <th>MOTIF_POS</th>
      <th>HIGH_INF_POS</th>
      <th>MOTIF_SCORE_CHANGE</th>
      <th>TRANSCRIPTION_FACTORS</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>28</th>
      <td>6</td>
      <td>160139865</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>0.697853</td>
      <td>0.699615</td>
      <td>-</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>HIGH</td>
      <td>SLC22A1</td>
      <td>...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
<p>1 rows × 84 columns</p>
</div>




```python
final = illumina_filter_cols[['chr', 'ref', 'alt', 'pos', 'Allele', 'IMPACT', 'SYMBOL', 'UAE_AF', 'Consequence']].reset_index() #to delete
final['Allele'] = final['Allele'].str.replace('CSQ=', '')
final = final.drop(columns = [ 'index'])#, "IMPACT_x", "index", "Consequence_x", 'SYMBOL_x', 'Consequence_y']) #keep 1st 3 only
final = final.drop_duplicates()
#final = final.rename(columns ={'MAX_AF_x':'MAX_AF_illu', 'MAX_AF_POPS_x': 'MAX_AF_POPS_illu','IMPACT_y':'IMPACT', 'Consequence_y':'Consequence'})
final = final.sort_values(by = ['IMPACT', 'chr'])
final.head(2)
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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>pos</th>
      <th>Allele</th>
      <th>IMPACT</th>
      <th>SYMBOL</th>
      <th>UAE_AF</th>
      <th>Consequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>44</th>
      <td>1</td>
      <td>G</td>
      <td>C</td>
      <td>13196863</td>
      <td>C</td>
      <td>HIGH</td>
      <td>PRAMEF13</td>
      <td>0.050288</td>
      <td>stop_gained</td>
    </tr>
    <tr>
      <th>92</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>HIGH</td>
      <td>PER3</td>
      <td>0.365414</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
    </tr>
  </tbody>
</table>
</div>



# Table 5 - final results


```python
def countConseq(conseq):
    if '&' in conseq:
        return len(conseq.split('&'))
    else:
        return 1
```


```python
final['number of Consequences'] = final['Consequence'].apply(countConseq)
final = final.sort_values(by = 'number of Consequences', ascending=False)
```


```python
cols_to_consider = [col for col in final.columns if ((col != 'Consequence')  & (col != 'number of Consequences'))]
final = final.drop_duplicates(subset=cols_to_consider)
final 
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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>pos</th>
      <th>Allele</th>
      <th>IMPACT</th>
      <th>SYMBOL</th>
      <th>UAE_AF</th>
      <th>Consequence</th>
      <th>number of Consequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>22</th>
      <td>6</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>160139865</td>
      <td>-</td>
      <td>HIGH</td>
      <td>SLC22A1</td>
      <td>0.699615</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>5</td>
    </tr>
    <tr>
      <th>182</th>
      <td>10</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGGCTTGGA...</td>
      <td>A</td>
      <td>132140569</td>
      <td>-</td>
      <td>HIGH</td>
      <td>JAKMIP3</td>
      <td>0.157294</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>69</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>HIGH</td>
      <td>PEX5</td>
      <td>0.733041</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>9</th>
      <td>18</td>
      <td>GGAAGCCATCCAGCCCAAGGAGGGTGACATCCCCAAGTCCCCAGAA</td>
      <td>G</td>
      <td>9887391</td>
      <td>-</td>
      <td>HIGH</td>
      <td>TXNDC2</td>
      <td>0.234015</td>
      <td>splice_acceptor_variant&amp;coding_sequence_varian...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>92</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>HIGH</td>
      <td>PER3</td>
      <td>0.365414</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>128</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>HIGH</td>
      <td>NPRL3</td>
      <td>0.436186</td>
      <td>splice_acceptor_variant&amp;5_prime_UTR_variant&amp;NM...</td>
      <td>3</td>
    </tr>
    <tr>
      <th>273</th>
      <td>22</td>
      <td>G</td>
      <td>C</td>
      <td>23998567</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.360359</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>376</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODERATE</td>
      <td></td>
      <td>0.059504</td>
      <td>TFBS_ablation&amp;TF_binding_site_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>268</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998560</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.358685</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>55</th>
      <td>15</td>
      <td>G</td>
      <td>A</td>
      <td>21872128</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>POTEB</td>
      <td>0.051766</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>102</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>241604325</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>OPN3</td>
      <td>0.055962</td>
      <td>intron_variant&amp;non_coding_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>233</th>
      <td>22</td>
      <td>T</td>
      <td>C</td>
      <td>23998556</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.356817</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>50</th>
      <td>14</td>
      <td>T</td>
      <td>C</td>
      <td>19428910</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>POTEG</td>
      <td>0.061430</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>156</th>
      <td>7</td>
      <td>C</td>
      <td>CCACCACCGACTCGATCGT</td>
      <td>100958459</td>
      <td>CACCACCGACTCGATCGT</td>
      <td>MODERATE</td>
      <td>MUC3A</td>
      <td>0.250235</td>
      <td>inframe_insertion&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>83</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>PEX5</td>
      <td>0.733041</td>
      <td>intron_variant&amp;non_coding_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>260</th>
      <td>22</td>
      <td>TA</td>
      <td>T</td>
      <td>23998571</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.354903</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>150</th>
      <td>7</td>
      <td>G</td>
      <td>T</td>
      <td>100958301</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>MUC3A</td>
      <td>0.555618</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>256</th>
      <td>22</td>
      <td>TC</td>
      <td>T</td>
      <td>23998574</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.345331</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>25</th>
      <td>6</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>160139865</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>SLC22A1</td>
      <td>0.699615</td>
      <td>intron_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>254</th>
      <td>22</td>
      <td>T</td>
      <td>TCG</td>
      <td>23998579</td>
      <td>CG</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.353516</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>253</th>
      <td>22</td>
      <td>A</td>
      <td>C</td>
      <td>23998558</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.359224</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>59</th>
      <td>7</td>
      <td>A</td>
      <td>T</td>
      <td>95316772</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>PON1</td>
      <td>0.343211</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>249</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998546</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.356771</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>237</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998550</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.359327</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>245</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998547</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.363351</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>48</th>
      <td>14</td>
      <td>A</td>
      <td>G</td>
      <td>18972839</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>POTEM</td>
      <td>0.061132</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>240</th>
      <td>22</td>
      <td>C</td>
      <td>T</td>
      <td>23998549</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.353871</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>282</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23441010</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.056581</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>124</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>NPRL3</td>
      <td>0.059504</td>
      <td>5_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>167</th>
      <td>16</td>
      <td>C</td>
      <td>CGCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTG...</td>
      <td>8645986</td>
      <td>GCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTGC...</td>
      <td>HIGH</td>
      <td>METTL22</td>
      <td>0.149167</td>
      <td>stop_gained&amp;frameshift_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>335</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODIFIER</td>
      <td>ADCK5</td>
      <td>0.594695</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>114</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>NPRL3</td>
      <td>0.436186</td>
      <td>5_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>284</th>
      <td>15</td>
      <td>T</td>
      <td>A</td>
      <td>23441015</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.050666</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>278</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23440310</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.050471</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>280</th>
      <td>15</td>
      <td>A</td>
      <td>ATCT</td>
      <td>23440898</td>
      <td>TCT</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.093766</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>160</th>
      <td>16</td>
      <td>C</td>
      <td>CGCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTG...</td>
      <td>8645986</td>
      <td>GCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTGC...</td>
      <td>MODIFIER</td>
      <td>METTL22</td>
      <td>0.149167</td>
      <td>intron_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>5</th>
      <td>7</td>
      <td>G</td>
      <td>GTT</td>
      <td>100773854</td>
      <td>TT</td>
      <td>HIGH</td>
      <td>ZAN</td>
      <td>0.499702</td>
      <td>frameshift_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>229</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>HBZ</td>
      <td>0.436186</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>357</th>
      <td>2</td>
      <td>C</td>
      <td>CGCT</td>
      <td>29065060</td>
      <td>GCT</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.339852</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>230</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>HBZ</td>
      <td>0.059504</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>316</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>AXIN1</td>
      <td>0.061831</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>350</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.436186</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>378</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.059504</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>173</th>
      <td>17</td>
      <td>A</td>
      <td>AGCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTA...</td>
      <td>40818851</td>
      <td>GCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTAT...</td>
      <td>MODIFIER</td>
      <td>KRT10-AS1</td>
      <td>0.063917</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>353</th>
      <td>17</td>
      <td>T</td>
      <td>TGCAGCAGGTGGTCAG</td>
      <td>41167977</td>
      <td>GCAGCAGGTGGTCAG</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.178833</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>354</th>
      <td>17</td>
      <td>C</td>
      <td>CGGTTCCATGGGCTCCGTA</td>
      <td>74893497</td>
      <td>GGTTCCATGGGCTCCGTA</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.211869</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>363</th>
      <td>17</td>
      <td>A</td>
      <td>AGCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTA...</td>
      <td>40818851</td>
      <td>GCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTAT...</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.063917</td>
      <td>TF_binding_site_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>18</td>
      <td>GGAAGCCATCCAGCCCAAGGAGGGTGACATCCCCAAGTCCCCAGAA</td>
      <td>G</td>
      <td>9887391</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>TXNDC2</td>
      <td>0.234015</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>101</th>
      <td>2</td>
      <td>C</td>
      <td>CGCT</td>
      <td>29065060</td>
      <td>GCT</td>
      <td>MODIFIER</td>
      <td>PCARE</td>
      <td>0.339852</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>320</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>ARHGDIG</td>
      <td>0.061831</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>66</th>
      <td>7</td>
      <td>A</td>
      <td>G</td>
      <td>75422819</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>POM121C</td>
      <td>0.215342</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>37</th>
      <td>5</td>
      <td>C</td>
      <td>G</td>
      <td>23527226</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>PRDM9</td>
      <td>0.102776</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>356</th>
      <td>6</td>
      <td>A</td>
      <td>C</td>
      <td>139135620</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.071964</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>362</th>
      <td>6</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>160139865</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.699615</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>0</th>
      <td>7</td>
      <td>G</td>
      <td>GTT</td>
      <td>100773854</td>
      <td>TT</td>
      <td>MODIFIER</td>
      <td>ZAN</td>
      <td>0.499702</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>29</th>
      <td>7</td>
      <td>C</td>
      <td>G</td>
      <td>142771014</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>PRSS2</td>
      <td>0.051881</td>
      <td>intron_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>58</th>
      <td>7</td>
      <td>A</td>
      <td>T</td>
      <td>95316772</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>PON1</td>
      <td>0.343211</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>153</th>
      <td>7</td>
      <td>A</td>
      <td>AGT</td>
      <td>100953372</td>
      <td>GT</td>
      <td>MODIFIER</td>
      <td>MUC3A</td>
      <td>0.135893</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>348</th>
      <td>7</td>
      <td>C</td>
      <td>CCACCACCGACTCGATCGT</td>
      <td>100958459</td>
      <td>CACCACCGACTCGATCGT</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.250235</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>349</th>
      <td>7</td>
      <td>C</td>
      <td>G</td>
      <td>142771014</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.051881</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>355</th>
      <td>7</td>
      <td>G</td>
      <td>T</td>
      <td>100958301</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.555618</td>
      <td>TF_binding_site_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>365</th>
      <td>7</td>
      <td>A</td>
      <td>T</td>
      <td>95316772</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.343211</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>158</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODIFIER</td>
      <td>MIR939</td>
      <td>0.594695</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>306</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODIFIER</td>
      <td>CPSF1</td>
      <td>0.594695</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>351</th>
      <td>8</td>
      <td>TTCCTTC</td>
      <td>T</td>
      <td>10610170</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.051973</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>33</th>
      <td>5</td>
      <td>A</td>
      <td>C</td>
      <td>23527276</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>PRDM9</td>
      <td>0.058908</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>374</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.079277</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>302</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>DDT</td>
      <td>0.079277</td>
      <td>3_prime_UTR_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>94</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>PDIA2</td>
      <td>0.061831</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>275</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>GSTT3P</td>
      <td>0.079277</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>297</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>DDTL</td>
      <td>0.079277</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>44</th>
      <td>1</td>
      <td>G</td>
      <td>C</td>
      <td>13196863</td>
      <td>C</td>
      <td>HIGH</td>
      <td>PRAMEF13</td>
      <td>0.050288</td>
      <td>stop_gained</td>
      <td>1</td>
    </tr>
    <tr>
      <th>347</th>
      <td>15</td>
      <td>C</td>
      <td>T</td>
      <td>82345024</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.069775</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>279</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23440310</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.050471</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>145</th>
      <td>11</td>
      <td>C</td>
      <td>T</td>
      <td>1191713</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>MUC5AC</td>
      <td>0.085135</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>147</th>
      <td>11</td>
      <td>A</td>
      <td>AACAACCTCTGCTCCTACAACCAGCACAACTTCTGTCCCTACAACC...</td>
      <td>1184994</td>
      <td>ACAACCTCTGCTCCTACAACCAGCACAACTTCTGTCCCTACAACCA...</td>
      <td>MODERATE</td>
      <td>MUC5AC</td>
      <td>0.053394</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>148</th>
      <td>11</td>
      <td>A</td>
      <td>G</td>
      <td>1191700</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>MUC5AC</td>
      <td>0.085995</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>71</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>PEX5</td>
      <td>0.733041</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>276</th>
      <td>15</td>
      <td>A</td>
      <td>G</td>
      <td>22462442</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>GOLGA6L22</td>
      <td>0.056661</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>277</th>
      <td>15</td>
      <td>T</td>
      <td>A</td>
      <td>23441015</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.050666</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>281</th>
      <td>15</td>
      <td>A</td>
      <td>ATCT</td>
      <td>23440898</td>
      <td>TCT</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.093766</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>287</th>
      <td>15</td>
      <td>C</td>
      <td>T</td>
      <td>82345024</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>GOLGA6L10</td>
      <td>0.069775</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>283</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23441010</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.056581</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>286</th>
      <td>15</td>
      <td>C</td>
      <td>T</td>
      <td>82345024</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>GOLGA6L10</td>
      <td>0.069775</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>288</th>
      <td>15</td>
      <td>A</td>
      <td>T</td>
      <td>23130350</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>GOLGA6L1</td>
      <td>0.147344</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>93</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>PDIA2</td>
      <td>0.061831</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>13</th>
      <td>17</td>
      <td>C</td>
      <td>G</td>
      <td>36262937</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>TBC1D3I</td>
      <td>0.206940</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>172</th>
      <td>17</td>
      <td>T</td>
      <td>TGCAGCAGGTGGTCAG</td>
      <td>41167977</td>
      <td>GCAGCAGGTGGTCAG</td>
      <td>MODERATE</td>
      <td>KRTAP4-3</td>
      <td>0.178833</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>180</th>
      <td>17</td>
      <td>A</td>
      <td>AGCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTA...</td>
      <td>40818851</td>
      <td>GCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTAT...</td>
      <td>MODERATE</td>
      <td>KRT10</td>
      <td>0.063917</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>289</th>
      <td>17</td>
      <td>C</td>
      <td>CGGTTCCATGGGCTCCGTA</td>
      <td>74893497</td>
      <td>GGTTCCATGGGCTCCGTA</td>
      <td>MODERATE</td>
      <td>FADS6</td>
      <td>0.211869</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>18</td>
      <td>GGAAGCCATCCAGCCCAAGGAGGGTGACATCCCCAAGTCCCCAGAA</td>
      <td>G</td>
      <td>9887391</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>TXNDC2</td>
      <td>0.234015</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>324</th>
      <td>10</td>
      <td>G</td>
      <td>T</td>
      <td>46385403</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>ANXA8L1</td>
      <td>0.050207</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>189</th>
      <td>10</td>
      <td>A</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>132140569</td>
      <td>CTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>MODERATE</td>
      <td>JAKMIP3</td>
      <td>0.053761</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>109</th>
      <td>10</td>
      <td>A</td>
      <td>G</td>
      <td>46461918</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NPY4R</td>
      <td>0.345801</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>20</th>
      <td>10</td>
      <td>T</td>
      <td>G</td>
      <td>68827963</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>STOX1</td>
      <td>0.168665</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>146</th>
      <td>11</td>
      <td>A</td>
      <td>AAC</td>
      <td>1185062</td>
      <td>AC</td>
      <td>HIGH</td>
      <td>MUC5AC</td>
      <td>0.100277</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>252</th>
      <td>22</td>
      <td>TC</td>
      <td>T</td>
      <td>23998574</td>
      <td>-</td>
      <td>HIGH</td>
      <td>GSTT4</td>
      <td>0.345331</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>257</th>
      <td>22</td>
      <td>TA</td>
      <td>T</td>
      <td>23998571</td>
      <td>-</td>
      <td>HIGH</td>
      <td>GSTT4</td>
      <td>0.354903</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>274</th>
      <td>22</td>
      <td>T</td>
      <td>TCG</td>
      <td>23998579</td>
      <td>CG</td>
      <td>HIGH</td>
      <td>GSTT4</td>
      <td>0.353516</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>152</th>
      <td>7</td>
      <td>A</td>
      <td>AGT</td>
      <td>100953372</td>
      <td>GT</td>
      <td>HIGH</td>
      <td>MUC3A</td>
      <td>0.135893</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>39</th>
      <td>1</td>
      <td>C</td>
      <td>T</td>
      <td>12942514</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>PRAMEF6</td>
      <td>0.094534</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>40</th>
      <td>1</td>
      <td>C</td>
      <td>T</td>
      <td>12941488</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>PRAMEF6</td>
      <td>0.070119</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>43</th>
      <td>1</td>
      <td>G</td>
      <td>A</td>
      <td>13259315</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>PRAMEF5</td>
      <td>0.065338</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>46</th>
      <td>1</td>
      <td>T</td>
      <td>C</td>
      <td>144593367</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>PPIAL4F</td>
      <td>0.091268</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>89</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>PER3</td>
      <td>0.365414</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>104</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>241604325</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>OPN3</td>
      <td>0.055962</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>131</th>
      <td>1</td>
      <td>G</td>
      <td>A</td>
      <td>145338502</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.094317</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>132</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>145356828</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.121197</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>134</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>145371089</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.138013</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>137</th>
      <td>1</td>
      <td>T</td>
      <td>A</td>
      <td>145302789</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.074669</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>139</th>
      <td>1</td>
      <td>C</td>
      <td>G</td>
      <td>146076743</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NBPF10</td>
      <td>0.067345</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>140</th>
      <td>1</td>
      <td>T</td>
      <td>A</td>
      <td>146085406</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>NBPF10</td>
      <td>0.076767</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>142</th>
      <td>1</td>
      <td>G</td>
      <td>GGTCTTC</td>
      <td>146090900</td>
      <td>GTCTTC</td>
      <td>MODERATE</td>
      <td>NBPF10</td>
      <td>0.322050</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>67</th>
      <td>19</td>
      <td>T</td>
      <td>C</td>
      <td>4512914</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>PLIN4</td>
      <td>0.051893</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>100</th>
      <td>2</td>
      <td>C</td>
      <td>CGCT</td>
      <td>29065060</td>
      <td>GCT</td>
      <td>MODERATE</td>
      <td>PCARE</td>
      <td>0.339852</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>234</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998550</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.359327</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>359</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>241604325</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.055962</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>360</th>
      <td>1</td>
      <td>T</td>
      <td>C</td>
      <td>144593367</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.091268</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>361</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.365414</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>16</th>
      <td>10</td>
      <td>T</td>
      <td>G</td>
      <td>68827963</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>STOX1</td>
      <td>0.168665</td>
      <td>intron_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>194</th>
      <td>10</td>
      <td>A</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>132140569</td>
      <td>CTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>MODIFIER</td>
      <td>JAKMIP3</td>
      <td>0.053761</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>222</th>
      <td>10</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGGCTTGGA...</td>
      <td>A</td>
      <td>132140569</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>JAKMIP3</td>
      <td>0.157294</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>327</th>
      <td>10</td>
      <td>G</td>
      <td>T</td>
      <td>46385403</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>ANXA8L1</td>
      <td>0.050207</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>337</th>
      <td>10</td>
      <td>T</td>
      <td>G</td>
      <td>68827963</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.168665</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>338</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.733041</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>341</th>
      <td>14</td>
      <td>A</td>
      <td>G</td>
      <td>18972839</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.061132</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>342</th>
      <td>14</td>
      <td>T</td>
      <td>C</td>
      <td>19428910</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.061430</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>53</th>
      <td>15</td>
      <td>G</td>
      <td>A</td>
      <td>21872128</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>POTEB</td>
      <td>0.051766</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>143</th>
      <td>1</td>
      <td>G</td>
      <td>GGTCTTC</td>
      <td>146090900</td>
      <td>GTCTTC</td>
      <td>MODIFIER</td>
      <td>NBPF10</td>
      <td>0.322050</td>
      <td>5_prime_UTR_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>332</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODERATE</td>
      <td>ADCK5</td>
      <td>0.594695</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>238</th>
      <td>22</td>
      <td>C</td>
      <td>T</td>
      <td>23998549</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.353871</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>294</th>
      <td>8</td>
      <td>T</td>
      <td>G</td>
      <td>143865158</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>EPPK1</td>
      <td>0.070508</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>242</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998547</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.363351</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>246</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998546</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.356771</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>251</th>
      <td>22</td>
      <td>T</td>
      <td>C</td>
      <td>23998556</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.356817</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>261</th>
      <td>22</td>
      <td>G</td>
      <td>C</td>
      <td>23998567</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.360359</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>265</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998560</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.358685</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>270</th>
      <td>22</td>
      <td>A</td>
      <td>C</td>
      <td>23998558</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.359224</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>305</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>DDT</td>
      <td>0.079277</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>34</th>
      <td>5</td>
      <td>A</td>
      <td>C</td>
      <td>23527276</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>PRDM9</td>
      <td>0.058908</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>35</th>
      <td>5</td>
      <td>C</td>
      <td>G</td>
      <td>23527226</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>PRDM9</td>
      <td>0.102776</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>228</th>
      <td>6</td>
      <td>A</td>
      <td>C</td>
      <td>139135620</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>HECA</td>
      <td>0.071964</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>28</th>
      <td>7</td>
      <td>C</td>
      <td>G</td>
      <td>142771014</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>PRSS2</td>
      <td>0.051881</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>60</th>
      <td>7</td>
      <td>A</td>
      <td>G</td>
      <td>75422819</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>POM121C</td>
      <td>0.215342</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>26</th>
      <td>8</td>
      <td>TTCCTTC</td>
      <td>T</td>
      <td>10610170</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>RP1L1</td>
      <td>0.051973</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>27</th>
      <td>8</td>
      <td>TCTC</td>
      <td>T</td>
      <td>10610160</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>RP1L1</td>
      <td>0.058232</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>291</th>
      <td>8</td>
      <td>T</td>
      <td>C</td>
      <td>143864555</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>EPPK1</td>
      <td>0.082177</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>292</th>
      <td>8</td>
      <td>T</td>
      <td>C</td>
      <td>143858843</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>EPPK1</td>
      <td>0.214757</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>352</th>
      <td>8</td>
      <td>TCTC</td>
      <td>T</td>
      <td>10610160</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.058232</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
print("Number of varinats: " + str(len(final)))
final
```

    Number of varinats: 149
    




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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>pos</th>
      <th>Allele</th>
      <th>IMPACT</th>
      <th>SYMBOL</th>
      <th>UAE_AF</th>
      <th>Consequence</th>
      <th>number of Consequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>22</th>
      <td>6</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>160139865</td>
      <td>-</td>
      <td>HIGH</td>
      <td>SLC22A1</td>
      <td>0.699615</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>5</td>
    </tr>
    <tr>
      <th>182</th>
      <td>10</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGGCTTGGA...</td>
      <td>A</td>
      <td>132140569</td>
      <td>-</td>
      <td>HIGH</td>
      <td>JAKMIP3</td>
      <td>0.157294</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>69</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>HIGH</td>
      <td>PEX5</td>
      <td>0.733041</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>9</th>
      <td>18</td>
      <td>GGAAGCCATCCAGCCCAAGGAGGGTGACATCCCCAAGTCCCCAGAA</td>
      <td>G</td>
      <td>9887391</td>
      <td>-</td>
      <td>HIGH</td>
      <td>TXNDC2</td>
      <td>0.234015</td>
      <td>splice_acceptor_variant&amp;coding_sequence_varian...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>92</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>HIGH</td>
      <td>PER3</td>
      <td>0.365414</td>
      <td>splice_donor_variant&amp;splice_donor_5th_base_var...</td>
      <td>4</td>
    </tr>
    <tr>
      <th>128</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>HIGH</td>
      <td>NPRL3</td>
      <td>0.436186</td>
      <td>splice_acceptor_variant&amp;5_prime_UTR_variant&amp;NM...</td>
      <td>3</td>
    </tr>
    <tr>
      <th>273</th>
      <td>22</td>
      <td>G</td>
      <td>C</td>
      <td>23998567</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.360359</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>376</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODERATE</td>
      <td></td>
      <td>0.059504</td>
      <td>TFBS_ablation&amp;TF_binding_site_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>268</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998560</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.358685</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>55</th>
      <td>15</td>
      <td>G</td>
      <td>A</td>
      <td>21872128</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>POTEB</td>
      <td>0.051766</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>102</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>241604325</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>OPN3</td>
      <td>0.055962</td>
      <td>intron_variant&amp;non_coding_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>233</th>
      <td>22</td>
      <td>T</td>
      <td>C</td>
      <td>23998556</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.356817</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>50</th>
      <td>14</td>
      <td>T</td>
      <td>C</td>
      <td>19428910</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>POTEG</td>
      <td>0.061430</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>156</th>
      <td>7</td>
      <td>C</td>
      <td>CCACCACCGACTCGATCGT</td>
      <td>100958459</td>
      <td>CACCACCGACTCGATCGT</td>
      <td>MODERATE</td>
      <td>MUC3A</td>
      <td>0.250235</td>
      <td>inframe_insertion&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>83</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>PEX5</td>
      <td>0.733041</td>
      <td>intron_variant&amp;non_coding_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>260</th>
      <td>22</td>
      <td>TA</td>
      <td>T</td>
      <td>23998571</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.354903</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>150</th>
      <td>7</td>
      <td>G</td>
      <td>T</td>
      <td>100958301</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>MUC3A</td>
      <td>0.555618</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>256</th>
      <td>22</td>
      <td>TC</td>
      <td>T</td>
      <td>23998574</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.345331</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>25</th>
      <td>6</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>160139865</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>SLC22A1</td>
      <td>0.699615</td>
      <td>intron_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>254</th>
      <td>22</td>
      <td>T</td>
      <td>TCG</td>
      <td>23998579</td>
      <td>CG</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.353516</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>253</th>
      <td>22</td>
      <td>A</td>
      <td>C</td>
      <td>23998558</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.359224</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>59</th>
      <td>7</td>
      <td>A</td>
      <td>T</td>
      <td>95316772</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>PON1</td>
      <td>0.343211</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>249</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998546</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.356771</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>237</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998550</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.359327</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>245</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998547</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.363351</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>48</th>
      <td>14</td>
      <td>A</td>
      <td>G</td>
      <td>18972839</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>POTEM</td>
      <td>0.061132</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>240</th>
      <td>22</td>
      <td>C</td>
      <td>T</td>
      <td>23998549</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>GSTT4</td>
      <td>0.353871</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>282</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23441010</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.056581</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>124</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>NPRL3</td>
      <td>0.059504</td>
      <td>5_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>167</th>
      <td>16</td>
      <td>C</td>
      <td>CGCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTG...</td>
      <td>8645986</td>
      <td>GCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTGC...</td>
      <td>HIGH</td>
      <td>METTL22</td>
      <td>0.149167</td>
      <td>stop_gained&amp;frameshift_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>335</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODIFIER</td>
      <td>ADCK5</td>
      <td>0.594695</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>114</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>NPRL3</td>
      <td>0.436186</td>
      <td>5_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>284</th>
      <td>15</td>
      <td>T</td>
      <td>A</td>
      <td>23441015</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.050666</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>278</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23440310</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.050471</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>280</th>
      <td>15</td>
      <td>A</td>
      <td>ATCT</td>
      <td>23440898</td>
      <td>TCT</td>
      <td>MODIFIER</td>
      <td>GOLGA6L2</td>
      <td>0.093766</td>
      <td>3_prime_UTR_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>160</th>
      <td>16</td>
      <td>C</td>
      <td>CGCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTG...</td>
      <td>8645986</td>
      <td>GCCTGCTCCGTGGCTGACATGTTGTCTAACTACTCTCCTTGCCTGC...</td>
      <td>MODIFIER</td>
      <td>METTL22</td>
      <td>0.149167</td>
      <td>intron_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>5</th>
      <td>7</td>
      <td>G</td>
      <td>GTT</td>
      <td>100773854</td>
      <td>TT</td>
      <td>HIGH</td>
      <td>ZAN</td>
      <td>0.499702</td>
      <td>frameshift_variant&amp;NMD_transcript_variant</td>
      <td>2</td>
    </tr>
    <tr>
      <th>229</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>HBZ</td>
      <td>0.436186</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>357</th>
      <td>2</td>
      <td>C</td>
      <td>CGCT</td>
      <td>29065060</td>
      <td>GCT</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.339852</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>230</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>HBZ</td>
      <td>0.059504</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>316</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>AXIN1</td>
      <td>0.061831</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>350</th>
      <td>16</td>
      <td>AGAGGAGGACGGAGCCGGAGGCGGAGGGGGCCT</td>
      <td>A</td>
      <td>138304</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.436186</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>378</th>
      <td>16</td>
      <td>GAGGGCAGGGGTGGAGGCGGAGGGGGCCTGAGGAGGGCAGGGGTGG...</td>
      <td>G</td>
      <td>138467</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.059504</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>173</th>
      <td>17</td>
      <td>A</td>
      <td>AGCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTA...</td>
      <td>40818851</td>
      <td>GCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTAT...</td>
      <td>MODIFIER</td>
      <td>KRT10-AS1</td>
      <td>0.063917</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>353</th>
      <td>17</td>
      <td>T</td>
      <td>TGCAGCAGGTGGTCAG</td>
      <td>41167977</td>
      <td>GCAGCAGGTGGTCAG</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.178833</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>354</th>
      <td>17</td>
      <td>C</td>
      <td>CGGTTCCATGGGCTCCGTA</td>
      <td>74893497</td>
      <td>GGTTCCATGGGCTCCGTA</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.211869</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>363</th>
      <td>17</td>
      <td>A</td>
      <td>AGCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTA...</td>
      <td>40818851</td>
      <td>GCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTAT...</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.063917</td>
      <td>TF_binding_site_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>18</td>
      <td>GGAAGCCATCCAGCCCAAGGAGGGTGACATCCCCAAGTCCCCAGAA</td>
      <td>G</td>
      <td>9887391</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>TXNDC2</td>
      <td>0.234015</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>101</th>
      <td>2</td>
      <td>C</td>
      <td>CGCT</td>
      <td>29065060</td>
      <td>GCT</td>
      <td>MODIFIER</td>
      <td>PCARE</td>
      <td>0.339852</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>320</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>ARHGDIG</td>
      <td>0.061831</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>66</th>
      <td>7</td>
      <td>A</td>
      <td>G</td>
      <td>75422819</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>POM121C</td>
      <td>0.215342</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>37</th>
      <td>5</td>
      <td>C</td>
      <td>G</td>
      <td>23527226</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>PRDM9</td>
      <td>0.102776</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>356</th>
      <td>6</td>
      <td>A</td>
      <td>C</td>
      <td>139135620</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.071964</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>362</th>
      <td>6</td>
      <td>CTGGTAAGT</td>
      <td>C</td>
      <td>160139865</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.699615</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>0</th>
      <td>7</td>
      <td>G</td>
      <td>GTT</td>
      <td>100773854</td>
      <td>TT</td>
      <td>MODIFIER</td>
      <td>ZAN</td>
      <td>0.499702</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>29</th>
      <td>7</td>
      <td>C</td>
      <td>G</td>
      <td>142771014</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>PRSS2</td>
      <td>0.051881</td>
      <td>intron_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>58</th>
      <td>7</td>
      <td>A</td>
      <td>T</td>
      <td>95316772</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>PON1</td>
      <td>0.343211</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>153</th>
      <td>7</td>
      <td>A</td>
      <td>AGT</td>
      <td>100953372</td>
      <td>GT</td>
      <td>MODIFIER</td>
      <td>MUC3A</td>
      <td>0.135893</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>348</th>
      <td>7</td>
      <td>C</td>
      <td>CCACCACCGACTCGATCGT</td>
      <td>100958459</td>
      <td>CACCACCGACTCGATCGT</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.250235</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>349</th>
      <td>7</td>
      <td>C</td>
      <td>G</td>
      <td>142771014</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.051881</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>355</th>
      <td>7</td>
      <td>G</td>
      <td>T</td>
      <td>100958301</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.555618</td>
      <td>TF_binding_site_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>365</th>
      <td>7</td>
      <td>A</td>
      <td>T</td>
      <td>95316772</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.343211</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>158</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODIFIER</td>
      <td>MIR939</td>
      <td>0.594695</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>306</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODIFIER</td>
      <td>CPSF1</td>
      <td>0.594695</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>351</th>
      <td>8</td>
      <td>TTCCTTC</td>
      <td>T</td>
      <td>10610170</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.051973</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>33</th>
      <td>5</td>
      <td>A</td>
      <td>C</td>
      <td>23527276</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td>PRDM9</td>
      <td>0.058908</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>374</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.079277</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>302</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>DDT</td>
      <td>0.079277</td>
      <td>3_prime_UTR_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>94</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>PDIA2</td>
      <td>0.061831</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>275</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>GSTT3P</td>
      <td>0.079277</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>297</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>DDTL</td>
      <td>0.079277</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>44</th>
      <td>1</td>
      <td>G</td>
      <td>C</td>
      <td>13196863</td>
      <td>C</td>
      <td>HIGH</td>
      <td>PRAMEF13</td>
      <td>0.050288</td>
      <td>stop_gained</td>
      <td>1</td>
    </tr>
    <tr>
      <th>347</th>
      <td>15</td>
      <td>C</td>
      <td>T</td>
      <td>82345024</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.069775</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>279</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23440310</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.050471</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>145</th>
      <td>11</td>
      <td>C</td>
      <td>T</td>
      <td>1191713</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>MUC5AC</td>
      <td>0.085135</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>147</th>
      <td>11</td>
      <td>A</td>
      <td>AACAACCTCTGCTCCTACAACCAGCACAACTTCTGTCCCTACAACC...</td>
      <td>1184994</td>
      <td>ACAACCTCTGCTCCTACAACCAGCACAACTTCTGTCCCTACAACCA...</td>
      <td>MODERATE</td>
      <td>MUC5AC</td>
      <td>0.053394</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>148</th>
      <td>11</td>
      <td>A</td>
      <td>G</td>
      <td>1191700</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>MUC5AC</td>
      <td>0.085995</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>71</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>PEX5</td>
      <td>0.733041</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>276</th>
      <td>15</td>
      <td>A</td>
      <td>G</td>
      <td>22462442</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>GOLGA6L22</td>
      <td>0.056661</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>277</th>
      <td>15</td>
      <td>T</td>
      <td>A</td>
      <td>23441015</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.050666</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>281</th>
      <td>15</td>
      <td>A</td>
      <td>ATCT</td>
      <td>23440898</td>
      <td>TCT</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.093766</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>287</th>
      <td>15</td>
      <td>C</td>
      <td>T</td>
      <td>82345024</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>GOLGA6L10</td>
      <td>0.069775</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>283</th>
      <td>15</td>
      <td>G</td>
      <td>C</td>
      <td>23441010</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GOLGA6L2</td>
      <td>0.056581</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>286</th>
      <td>15</td>
      <td>C</td>
      <td>T</td>
      <td>82345024</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>GOLGA6L10</td>
      <td>0.069775</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>288</th>
      <td>15</td>
      <td>A</td>
      <td>T</td>
      <td>23130350</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>GOLGA6L1</td>
      <td>0.147344</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>93</th>
      <td>16</td>
      <td>G</td>
      <td>A</td>
      <td>286898</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>PDIA2</td>
      <td>0.061831</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>13</th>
      <td>17</td>
      <td>C</td>
      <td>G</td>
      <td>36262937</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>TBC1D3I</td>
      <td>0.206940</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>172</th>
      <td>17</td>
      <td>T</td>
      <td>TGCAGCAGGTGGTCAG</td>
      <td>41167977</td>
      <td>GCAGCAGGTGGTCAG</td>
      <td>MODERATE</td>
      <td>KRTAP4-3</td>
      <td>0.178833</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>180</th>
      <td>17</td>
      <td>A</td>
      <td>AGCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTA...</td>
      <td>40818851</td>
      <td>GCTGCCGCCGCCGTATCCGCCGCCGGAGCTGCTGCCGCCGCCGTAT...</td>
      <td>MODERATE</td>
      <td>KRT10</td>
      <td>0.063917</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>289</th>
      <td>17</td>
      <td>C</td>
      <td>CGGTTCCATGGGCTCCGTA</td>
      <td>74893497</td>
      <td>GGTTCCATGGGCTCCGTA</td>
      <td>MODERATE</td>
      <td>FADS6</td>
      <td>0.211869</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>18</td>
      <td>GGAAGCCATCCAGCCCAAGGAGGGTGACATCCCCAAGTCCCCAGAA</td>
      <td>G</td>
      <td>9887391</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>TXNDC2</td>
      <td>0.234015</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>324</th>
      <td>10</td>
      <td>G</td>
      <td>T</td>
      <td>46385403</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>ANXA8L1</td>
      <td>0.050207</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>189</th>
      <td>10</td>
      <td>A</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>132140569</td>
      <td>CTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>MODERATE</td>
      <td>JAKMIP3</td>
      <td>0.053761</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>109</th>
      <td>10</td>
      <td>A</td>
      <td>G</td>
      <td>46461918</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NPY4R</td>
      <td>0.345801</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>20</th>
      <td>10</td>
      <td>T</td>
      <td>G</td>
      <td>68827963</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>STOX1</td>
      <td>0.168665</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>146</th>
      <td>11</td>
      <td>A</td>
      <td>AAC</td>
      <td>1185062</td>
      <td>AC</td>
      <td>HIGH</td>
      <td>MUC5AC</td>
      <td>0.100277</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>252</th>
      <td>22</td>
      <td>TC</td>
      <td>T</td>
      <td>23998574</td>
      <td>-</td>
      <td>HIGH</td>
      <td>GSTT4</td>
      <td>0.345331</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>257</th>
      <td>22</td>
      <td>TA</td>
      <td>T</td>
      <td>23998571</td>
      <td>-</td>
      <td>HIGH</td>
      <td>GSTT4</td>
      <td>0.354903</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>274</th>
      <td>22</td>
      <td>T</td>
      <td>TCG</td>
      <td>23998579</td>
      <td>CG</td>
      <td>HIGH</td>
      <td>GSTT4</td>
      <td>0.353516</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>152</th>
      <td>7</td>
      <td>A</td>
      <td>AGT</td>
      <td>100953372</td>
      <td>GT</td>
      <td>HIGH</td>
      <td>MUC3A</td>
      <td>0.135893</td>
      <td>frameshift_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>39</th>
      <td>1</td>
      <td>C</td>
      <td>T</td>
      <td>12942514</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>PRAMEF6</td>
      <td>0.094534</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>40</th>
      <td>1</td>
      <td>C</td>
      <td>T</td>
      <td>12941488</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>PRAMEF6</td>
      <td>0.070119</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>43</th>
      <td>1</td>
      <td>G</td>
      <td>A</td>
      <td>13259315</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>PRAMEF5</td>
      <td>0.065338</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>46</th>
      <td>1</td>
      <td>T</td>
      <td>C</td>
      <td>144593367</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>PPIAL4F</td>
      <td>0.091268</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>89</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>PER3</td>
      <td>0.365414</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>104</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>241604325</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>OPN3</td>
      <td>0.055962</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>131</th>
      <td>1</td>
      <td>G</td>
      <td>A</td>
      <td>145338502</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.094317</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>132</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>145356828</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.121197</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>134</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>145371089</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.138013</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>137</th>
      <td>1</td>
      <td>T</td>
      <td>A</td>
      <td>145302789</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>NBPF20</td>
      <td>0.074669</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>139</th>
      <td>1</td>
      <td>C</td>
      <td>G</td>
      <td>146076743</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>NBPF10</td>
      <td>0.067345</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>140</th>
      <td>1</td>
      <td>T</td>
      <td>A</td>
      <td>146085406</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>NBPF10</td>
      <td>0.076767</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>142</th>
      <td>1</td>
      <td>G</td>
      <td>GGTCTTC</td>
      <td>146090900</td>
      <td>GTCTTC</td>
      <td>MODERATE</td>
      <td>NBPF10</td>
      <td>0.322050</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>67</th>
      <td>19</td>
      <td>T</td>
      <td>C</td>
      <td>4512914</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>PLIN4</td>
      <td>0.051893</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>100</th>
      <td>2</td>
      <td>C</td>
      <td>CGCT</td>
      <td>29065060</td>
      <td>GCT</td>
      <td>MODERATE</td>
      <td>PCARE</td>
      <td>0.339852</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>234</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998550</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.359327</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>359</th>
      <td>1</td>
      <td>A</td>
      <td>G</td>
      <td>241604325</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.055962</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>360</th>
      <td>1</td>
      <td>T</td>
      <td>C</td>
      <td>144593367</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.091268</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>361</th>
      <td>1</td>
      <td>GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGC...</td>
      <td>G</td>
      <td>7829912</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.365414</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>16</th>
      <td>10</td>
      <td>T</td>
      <td>G</td>
      <td>68827963</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td>STOX1</td>
      <td>0.168665</td>
      <td>intron_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>194</th>
      <td>10</td>
      <td>A</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>132140569</td>
      <td>CTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGG</td>
      <td>MODIFIER</td>
      <td>JAKMIP3</td>
      <td>0.053761</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>222</th>
      <td>10</td>
      <td>ACTTGGAGGAGGTAACGAGGGTCTCCTGCCGGGTCCTGGGCTTGGA...</td>
      <td>A</td>
      <td>132140569</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td>JAKMIP3</td>
      <td>0.157294</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>327</th>
      <td>10</td>
      <td>G</td>
      <td>T</td>
      <td>46385403</td>
      <td>T</td>
      <td>MODIFIER</td>
      <td>ANXA8L1</td>
      <td>0.050207</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>337</th>
      <td>10</td>
      <td>T</td>
      <td>G</td>
      <td>68827963</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.168665</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>338</th>
      <td>12</td>
      <td>GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA</td>
      <td>G</td>
      <td>7190512</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.733041</td>
      <td>upstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>341</th>
      <td>14</td>
      <td>A</td>
      <td>G</td>
      <td>18972839</td>
      <td>G</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.061132</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>342</th>
      <td>14</td>
      <td>T</td>
      <td>C</td>
      <td>19428910</td>
      <td>C</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.061430</td>
      <td>downstream_gene_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>53</th>
      <td>15</td>
      <td>G</td>
      <td>A</td>
      <td>21872128</td>
      <td>A</td>
      <td>MODIFIER</td>
      <td>POTEB</td>
      <td>0.051766</td>
      <td>non_coding_transcript_exon_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>143</th>
      <td>1</td>
      <td>G</td>
      <td>GGTCTTC</td>
      <td>146090900</td>
      <td>GTCTTC</td>
      <td>MODIFIER</td>
      <td>NBPF10</td>
      <td>0.322050</td>
      <td>5_prime_UTR_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>332</th>
      <td>8</td>
      <td>T</td>
      <td>TGGGGGTGCAAGGTGA</td>
      <td>144392334</td>
      <td>GGGGGTGCAAGGTGA</td>
      <td>MODERATE</td>
      <td>ADCK5</td>
      <td>0.594695</td>
      <td>inframe_insertion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>238</th>
      <td>22</td>
      <td>C</td>
      <td>T</td>
      <td>23998549</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.353871</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>294</th>
      <td>8</td>
      <td>T</td>
      <td>G</td>
      <td>143865158</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>EPPK1</td>
      <td>0.070508</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>242</th>
      <td>22</td>
      <td>T</td>
      <td>A</td>
      <td>23998547</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.363351</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>246</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998546</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.356771</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>251</th>
      <td>22</td>
      <td>T</td>
      <td>C</td>
      <td>23998556</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.356817</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>261</th>
      <td>22</td>
      <td>G</td>
      <td>C</td>
      <td>23998567</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.360359</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>265</th>
      <td>22</td>
      <td>C</td>
      <td>A</td>
      <td>23998560</td>
      <td>A</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.358685</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>270</th>
      <td>22</td>
      <td>A</td>
      <td>C</td>
      <td>23998558</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>GSTT4</td>
      <td>0.359224</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>305</th>
      <td>22</td>
      <td>A</td>
      <td>T</td>
      <td>23973190</td>
      <td>T</td>
      <td>MODERATE</td>
      <td>DDT</td>
      <td>0.079277</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>34</th>
      <td>5</td>
      <td>A</td>
      <td>C</td>
      <td>23527276</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>PRDM9</td>
      <td>0.058908</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>35</th>
      <td>5</td>
      <td>C</td>
      <td>G</td>
      <td>23527226</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>PRDM9</td>
      <td>0.102776</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>228</th>
      <td>6</td>
      <td>A</td>
      <td>C</td>
      <td>139135620</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>HECA</td>
      <td>0.071964</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>28</th>
      <td>7</td>
      <td>C</td>
      <td>G</td>
      <td>142771014</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>PRSS2</td>
      <td>0.051881</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>60</th>
      <td>7</td>
      <td>A</td>
      <td>G</td>
      <td>75422819</td>
      <td>G</td>
      <td>MODERATE</td>
      <td>POM121C</td>
      <td>0.215342</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>26</th>
      <td>8</td>
      <td>TTCCTTC</td>
      <td>T</td>
      <td>10610170</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>RP1L1</td>
      <td>0.051973</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>27</th>
      <td>8</td>
      <td>TCTC</td>
      <td>T</td>
      <td>10610160</td>
      <td>-</td>
      <td>MODERATE</td>
      <td>RP1L1</td>
      <td>0.058232</td>
      <td>inframe_deletion</td>
      <td>1</td>
    </tr>
    <tr>
      <th>291</th>
      <td>8</td>
      <td>T</td>
      <td>C</td>
      <td>143864555</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>EPPK1</td>
      <td>0.082177</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>292</th>
      <td>8</td>
      <td>T</td>
      <td>C</td>
      <td>143858843</td>
      <td>C</td>
      <td>MODERATE</td>
      <td>EPPK1</td>
      <td>0.214757</td>
      <td>missense_variant</td>
      <td>1</td>
    </tr>
    <tr>
      <th>352</th>
      <td>8</td>
      <td>TCTC</td>
      <td>T</td>
      <td>10610160</td>
      <td>-</td>
      <td>MODIFIER</td>
      <td></td>
      <td>0.058232</td>
      <td>regulatory_region_variant</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
final.to_csv("/home/jovyan/aisha/Paper/results/Table5_partA_77variants_alldata.csv", index = False)
```


```python
def readFile(fileName, chrCol):
    file = pd.read_csv(fileName, sep="\t", header=None)
    # new dataframe only for the INFO column
    file_INFO = file[7].str.split(";", expand=True)
    file_INFO[4] = file_INFO[4].str.split(',')
    file_INFO = file_INFO.explode(4)
    file_final = file_INFO[4].str.split('|', expand=True)
    file_final.columns = header
    file_final = file_final.drop(columns = cols_to_drop)
    file_final.insert(0, 'UAE_AF', file_INFO[1])
    file_final.insert(0, 'alt', file[4])
    file_final.insert(0, 'ref', file[3])
    file_final.insert(0, 'pos', file[1])
    file_final.insert(0, 'chr', file[chrCol])
    file_final = file_final.drop_duplicates()
    
    return file_final
```


```python
# only the new start positions of the variants (subset of Table5_variants_liftove_all.txt)
liftover_pos = pd.read_csv("positions/Table5_variants_liftove_all.txt", sep="\t", header=None, names=['chr','ref','alt','gene','pos','end'])
liftover_pos.head(1)
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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>gene</th>
      <th>pos</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr2</td>
      <td>C</td>
      <td>G</td>
      <td>ACADL</td>
      <td>210205677.0</td>
      <td>210205678.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
#  filter the illumina MAFtable using the new positions
illuMAF = pd.read_csv("positions/Table5_allVariants_pos_MAFtable.txt", sep=",", header=None, names=['filename','pos','end','ref','alt','referenceCounts','alternateCounts','ReferenceFreq','AlternateFreq'])
illuMAF = illuMAF.sort_values(by='pos')
illuMAF.insert(0, 'chr', illuMAF['filename'].str.split(":", expand=True)[1])
illuMAF = illuMAF.drop(columns=["filename", "referenceCounts", "alternateCounts","ReferenceFreq"])
illuMAF = illuMAF.rename(columns={'AlternateFreq': 'UAE_AF'})
illuMAF.head(1)
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
      <th>chr</th>
      <th>pos</th>
      <th>end</th>
      <th>ref</th>
      <th>alt</th>
      <th>UAE_AF</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>546</th>
      <td>chr16</td>
      <td>1018086</td>
      <td>1018087</td>
      <td>G</td>
      <td>A</td>
      <td>0.001364</td>
    </tr>
  </tbody>
</table>
</div>




```python
maf_Illumina_VEP_variants = readFile("testing_file_table5.txt", 0)
maf_Illumina_VEP_variants = maf_Illumina_VEP_variants[["chr","pos","ref","alt", "SYMBOL", "UAE_AF", "MAX_AF", "MAX_AF_POPS", "Existing_variation", "IMPACT", "Consequence"]]
print("number of entries found from the new positions of the variants in the paper: "+ str(len(maf_Illumina_VEP_variants)))
maf_Illumina_VEP_variants.sort_values(by='pos').head(1)
```

    number of entries found from the new positions of the variants in the paper: 5637
    




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
      <th>chr</th>
      <th>pos</th>
      <th>ref</th>
      <th>alt</th>
      <th>SYMBOL</th>
      <th>UAE_AF</th>
      <th>MAX_AF</th>
      <th>MAX_AF_POPS</th>
      <th>Existing_variation</th>
      <th>IMPACT</th>
      <th>Consequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>68</th>
      <td>1</td>
      <td>88674</td>
      <td>AACCCAGAGTGCTCGTGTTGGCCGGCAGAGCCGGCCCCCATCTCCT...</td>
      <td>A</td>
      <td></td>
      <td>AF=0.00163919392925101</td>
      <td></td>
      <td></td>
      <td></td>
      <td>MODIFIER</td>
      <td>downstream_gene_variant</td>
    </tr>
  </tbody>
</table>
</div>




```python
#  filter the mgi MAFtable using the new positions
mgiMAF = pd.read_csv("positions/Table5_mgi.txt", sep=",", header=None, names=['filename','start','end','referenceAllele','alternateAlleles','referenceCounts','alternateCounts','ReferenceFreq','AlternateFreq'])
mgiMAF = mgiMAF.sort_values(by='start')
mgiMAF = mgiMAF.rename(columns = {'filename': 'chr'})
mgiMAF.head(1)
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
      <th>chr</th>
      <th>start</th>
      <th>end</th>
      <th>referenceAllele</th>
      <th>alternateAlleles</th>
      <th>referenceCounts</th>
      <th>alternateCounts</th>
      <th>ReferenceFreq</th>
      <th>AlternateFreq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>83</th>
      <td>chr11</td>
      <td>1018086</td>
      <td>1018087</td>
      <td>G</td>
      <td>A</td>
      <td>12761</td>
      <td>1</td>
      <td>0.999922</td>
      <td>0.000078</td>
    </tr>
  </tbody>
</table>
</div>




```python
# make sure the tables have the same headers and values formats
illuMAF['chr'] = illuMAF['chr'].str.replace('chr', '')
liftover_pos['chr'] = liftover_pos['chr'].str.replace('chr', '')
liftover_pos = liftover_pos[:-2]
liftover_pos['pos'] = liftover_pos['pos'].astype(int)
liftover_pos['end'] = liftover_pos['end'].astype(int)

maf_Illumina_VEP_variants['UAE_AF'] = maf_Illumina_VEP_variants['UAE_AF'].str.replace('AF=', '')

illuMAF['pos'] = illuMAF['pos'].astype(int)+1
liftover_pos['pos'] = liftover_pos['pos'].astype(int)+1
```

    /tmp/ipykernel_59/2303359823.py:5: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      liftover_pos['pos'] = liftover_pos['pos'].astype(int)
    /tmp/ipykernel_59/2303359823.py:6: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      liftover_pos['end'] = liftover_pos['end'].astype(int)
    /tmp/ipykernel_59/2303359823.py:11: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      liftover_pos['pos'] = liftover_pos['pos'].astype(int)+1
    


```python
print("Total number of rows of each dataframe")

print("[illuMAF]: " + str(len(illuMAF)))
print("[liftover_pos]: " + str(len(liftover_pos)))
print("[maf_Illumina_VEP_variants]: " + str(len(maf_Illumina_VEP_variants)))
```

    Total number of rows of each dataframe
    [illuMAF]: 759
    [liftover_pos]: 44
    [maf_Illumina_VEP_variants]: 5637
    

## Merging


```python
# merge liftover_pos & maf_Illumina_VEP_variants
# in other terms, merge the dataframe of the new positions with the annotated dataframe
merge_VEP_liftoverPos = pd.merge(liftover_pos,maf_Illumina_VEP_variants, how='inner', on =['chr', 'pos', 'ref','alt'])
print("Total number of rows: " + str(len(merge_VEP_liftoverPos)))
merge_VEP_liftoverPos.head(1)
```

    Total number of rows: 73
    




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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>gene</th>
      <th>pos</th>
      <th>end</th>
      <th>SYMBOL</th>
      <th>UAE_AF</th>
      <th>MAX_AF</th>
      <th>MAX_AF_POPS</th>
      <th>Existing_variation</th>
      <th>IMPACT</th>
      <th>Consequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2</td>
      <td>C</td>
      <td>G</td>
      <td>ACADL</td>
      <td>210205678</td>
      <td>210205678</td>
      <td>ACADL</td>
      <td>0.0279923886379789</td>
      <td>0.008239</td>
      <td>gnomADe_ASJ</td>
      <td>rs146511220&amp;COSV106086073</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
  </tbody>
</table>
</div>




```python
# merge liftover_pos & maf_Illumina_VEP_variants
# in other terms, merge the dataframe of the new positions with the MAF dataframe
# why? because the annoatted dataframe have filtered out the variants with very low AF
merge_illuminaMAF_liftoverPos = pd.merge(liftover_pos,illuMAF, how='inner', on =['chr', 'pos', 'ref','alt'])
print("Length before dropping duplicates: " + str(len(merge_illuminaMAF_liftoverPos)))
merge_illuminaMAF_liftoverPos = merge_illuminaMAF_liftoverPos.drop_duplicates()
print("Length after dropping duplicates: " + str(len(merge_illuminaMAF_liftoverPos)))
merge_illuminaMAF_liftoverPos = merge_illuminaMAF_liftoverPos.sort_values(by='pos')
merge_illuminaMAF_liftoverPos.head(1)
```

    Length before dropping duplicates: 75
    Length after dropping duplicates: 28
    




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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>gene</th>
      <th>pos</th>
      <th>end_x</th>
      <th>end_y</th>
      <th>UAE_AF</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>73</th>
      <td>11</td>
      <td>TG</td>
      <td>T</td>
      <td>MUC6</td>
      <td>1018088</td>
      <td>1018089</td>
      <td>1018089</td>
      <td>0.000768</td>
    </tr>
  </tbody>
</table>
</div>




```python
# merge both daatframes together to get the final updated table
update_oldTable5 = pd.merge(merge_illuminaMAF_liftoverPos,merge_VEP_liftoverPos, how='outer', on =['chr', 'pos', 'ref', 'alt'])
print("Total number of rows: " + str(len(update_oldTable5)))
update_oldTable5.head(5)
```

    Total number of rows: 82
    




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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>gene_x</th>
      <th>pos</th>
      <th>end_x</th>
      <th>end_y</th>
      <th>UAE_AF_x</th>
      <th>gene_y</th>
      <th>end</th>
      <th>SYMBOL</th>
      <th>UAE_AF_y</th>
      <th>MAX_AF</th>
      <th>MAX_AF_POPS</th>
      <th>Existing_variation</th>
      <th>IMPACT</th>
      <th>Consequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>11</td>
      <td>TG</td>
      <td>T</td>
      <td>MUC6</td>
      <td>1018088</td>
      <td>1018089</td>
      <td>1018089</td>
      <td>0.000768</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11</td>
      <td>G</td>
      <td>GCA</td>
      <td>MUC6</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>0.330968</td>
      <td>MUC6</td>
      <td>1018215.0</td>
      <td>MUC6</td>
      <td>0.330968156078773</td>
      <td>0.2744</td>
      <td>gnomADg_MID</td>
      <td>rs769713098&amp;COSV70144397</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
    <tr>
      <th>2</th>
      <td>11</td>
      <td>G</td>
      <td>GCA</td>
      <td>MUC6</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>0.330968</td>
      <td>MUC6</td>
      <td>1018215.0</td>
      <td>MUC6</td>
      <td>0.330968156078773</td>
      <td>0.2744</td>
      <td>gnomADg_MID</td>
      <td>rs769713098&amp;COSV70144397</td>
      <td>MODIFIER</td>
      <td>downstream_gene_variant</td>
    </tr>
    <tr>
      <th>3</th>
      <td>11</td>
      <td>G</td>
      <td>GCA</td>
      <td>MUC6</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>0.330968</td>
      <td>MUC6</td>
      <td>1018215.0</td>
      <td>MUC6</td>
      <td>0.330968156078773</td>
      <td>0.2744</td>
      <td>gnomADg_MID</td>
      <td>rs769713098&amp;COSV70144397</td>
      <td>MODIFIER</td>
      <td>upstream_gene_variant</td>
    </tr>
    <tr>
      <th>4</th>
      <td>11</td>
      <td>AAT</td>
      <td>A</td>
      <td>MUC6</td>
      <td>1018222</td>
      <td>1018224</td>
      <td>1018224</td>
      <td>0.314232</td>
      <td>MUC6</td>
      <td>1018224.0</td>
      <td>MUC6</td>
      <td>0.314232329947958</td>
      <td>0.04129</td>
      <td>gnomADg_FIN</td>
      <td>rs780061827</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
  </tbody>
</table>
</div>




```python
update_oldTable5_noMod = update_oldTable5.loc[update_oldTable5['IMPACT']!="MODIFIER"]
update_oldTable5_noMod
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
      <th>chr</th>
      <th>ref</th>
      <th>alt</th>
      <th>gene_x</th>
      <th>pos</th>
      <th>end_x</th>
      <th>end_y</th>
      <th>UAE_AF_x</th>
      <th>gene_y</th>
      <th>end</th>
      <th>SYMBOL</th>
      <th>UAE_AF_y</th>
      <th>MAX_AF</th>
      <th>MAX_AF_POPS</th>
      <th>Existing_variation</th>
      <th>IMPACT</th>
      <th>Consequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>11</td>
      <td>TG</td>
      <td>T</td>
      <td>MUC6</td>
      <td>1018088</td>
      <td>1018089</td>
      <td>1018089</td>
      <td>0.000768</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11</td>
      <td>G</td>
      <td>GCA</td>
      <td>MUC6</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>1018215</td>
      <td>0.330968</td>
      <td>MUC6</td>
      <td>1018215.0</td>
      <td>MUC6</td>
      <td>0.330968156078773</td>
      <td>0.2744</td>
      <td>gnomADg_MID</td>
      <td>rs769713098&amp;COSV70144397</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
    <tr>
      <th>4</th>
      <td>11</td>
      <td>AAT</td>
      <td>A</td>
      <td>MUC6</td>
      <td>1018222</td>
      <td>1018224</td>
      <td>1018224</td>
      <td>0.314232</td>
      <td>MUC6</td>
      <td>1018224.0</td>
      <td>MUC6</td>
      <td>0.314232329947958</td>
      <td>0.04129</td>
      <td>gnomADg_FIN</td>
      <td>rs780061827</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
    <tr>
      <th>7</th>
      <td>19</td>
      <td>GGGC</td>
      <td>G</td>
      <td>BTBD2</td>
      <td>2015541</td>
      <td>2015559</td>
      <td>2015544</td>
      <td>0.457198</td>
      <td>BTBD2</td>
      <td>2015559.0</td>
      <td>BTBD2</td>
      <td>0.457197551525711</td>
      <td>0.5076</td>
      <td>gnomADg_NFE</td>
      <td>rs533915623</td>
      <td>MODERATE</td>
      <td>inframe_deletion</td>
    </tr>
    <tr>
      <th>12</th>
      <td>X</td>
      <td>A</td>
      <td>T</td>
      <td>ARSD</td>
      <td>2918000</td>
      <td>2918000</td>
      <td>2918000</td>
      <td>0.000034</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>13</th>
      <td>X</td>
      <td>C</td>
      <td>T</td>
      <td>ARSD</td>
      <td>2918006</td>
      <td>2918006</td>
      <td>2918006</td>
      <td>0.000034</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>14</th>
      <td>X</td>
      <td>A</td>
      <td>T</td>
      <td>ARSD</td>
      <td>2918170</td>
      <td>2918170</td>
      <td>2918170</td>
      <td>0.000011</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>15</th>
      <td>X</td>
      <td>G</td>
      <td>A</td>
      <td>ARSD</td>
      <td>2918197</td>
      <td>2918197</td>
      <td>2918197</td>
      <td>0.000103</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>16</th>
      <td>7</td>
      <td>G</td>
      <td>A</td>
      <td>FAM220A</td>
      <td>6330718</td>
      <td>6330718</td>
      <td>6330718</td>
      <td>0.038343</td>
      <td>FAM220A</td>
      <td>6330718.0</td>
      <td>FAM220A</td>
      <td>0.0383433824709416</td>
      <td>0.006329</td>
      <td>gnomADg_MID</td>
      <td>rs75910050&amp;COSV57627270</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>22</th>
      <td>1</td>
      <td>CAGCTT</td>
      <td>C</td>
      <td>ESPN</td>
      <td>6445776</td>
      <td>6445781</td>
      <td>6445781</td>
      <td>0.003783</td>
      <td>ESPN</td>
      <td>6445781.0</td>
      <td>ESPN</td>
      <td>0.00378275522134849</td>
      <td>0.01669</td>
      <td>gnomADg_ASJ</td>
      <td>rs753994746</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
    <tr>
      <th>24</th>
      <td>17</td>
      <td>CTGT</td>
      <td>C</td>
      <td>KDM6B</td>
      <td>7847704</td>
      <td>7847707</td>
      <td>7847707</td>
      <td>0.000390</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>25</th>
      <td>1</td>
      <td>CTG</td>
      <td>C</td>
      <td>SPEN</td>
      <td>15935980</td>
      <td>15935982</td>
      <td>15935982</td>
      <td>0.000149</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>26</th>
      <td>2</td>
      <td>A</td>
      <td>G</td>
      <td>WDR35</td>
      <td>19973675</td>
      <td>19973675</td>
      <td>19973675</td>
      <td>0.050047</td>
      <td>WDR35</td>
      <td>19973675.0</td>
      <td>WDR35</td>
      <td>0.0500469978679016</td>
      <td>0.04747</td>
      <td>gnomADg_MID</td>
      <td>rs142955097</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>27</th>
      <td>2</td>
      <td>A</td>
      <td>G</td>
      <td>WDR35</td>
      <td>19973675</td>
      <td>19973675</td>
      <td>19973675</td>
      <td>0.050047</td>
      <td>WDR35</td>
      <td>19973675.0</td>
      <td>WDR35</td>
      <td>0.0500469978679016</td>
      <td>0.04747</td>
      <td>gnomADg_MID</td>
      <td>rs142955097</td>
      <td>MODERATE</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
    </tr>
    <tr>
      <th>28</th>
      <td>2</td>
      <td>GGGC</td>
      <td>G</td>
      <td>GDF7</td>
      <td>20667362</td>
      <td>20667374</td>
      <td>20667365</td>
      <td>0.031683</td>
      <td>GDF7</td>
      <td>20667374.0</td>
      <td>GDF7</td>
      <td>0.0316834407024462</td>
      <td>0.03488</td>
      <td>gnomADg_AFR</td>
      <td>rs563301908</td>
      <td>MODERATE</td>
      <td>inframe_deletion</td>
    </tr>
    <tr>
      <th>30</th>
      <td>2</td>
      <td>G</td>
      <td>A</td>
      <td>LTBP1</td>
      <td>33347499</td>
      <td>33347499</td>
      <td>33347499</td>
      <td>0.042161</td>
      <td>LTBP1</td>
      <td>33347499.0</td>
      <td>LTBP1</td>
      <td>0.0421605263761205</td>
      <td>0.02532</td>
      <td>gnomADg_MID</td>
      <td>rs141080282&amp;COSV106083540</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>32</th>
      <td>17</td>
      <td>C</td>
      <td>T</td>
      <td>CASC3</td>
      <td>40164128</td>
      <td>40164128</td>
      <td>40164128</td>
      <td>0.045657</td>
      <td>CASC3</td>
      <td>40164128.0</td>
      <td>CASC3</td>
      <td>0.0456567092322153</td>
      <td>0.03481</td>
      <td>gnomADg_MID</td>
      <td>rs140375987</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>36</th>
      <td>17</td>
      <td>C</td>
      <td>T</td>
      <td>CASC3</td>
      <td>40164128</td>
      <td>40164128</td>
      <td>40164128</td>
      <td>0.045657</td>
      <td>CASC3</td>
      <td>40164128.0</td>
      <td>CASC3</td>
      <td>0.0456567092322153</td>
      <td>0.03481</td>
      <td>gnomADg_MID</td>
      <td>rs140375987</td>
      <td>MODERATE</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
    </tr>
    <tr>
      <th>39</th>
      <td>19</td>
      <td>A</td>
      <td>T</td>
      <td>ZNF28</td>
      <td>52800274</td>
      <td>52800274</td>
      <td>52800274</td>
      <td>0.000298</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>40</th>
      <td>16</td>
      <td>C</td>
      <td>T</td>
      <td>ACD</td>
      <td>67660199</td>
      <td>67660199</td>
      <td>67660199</td>
      <td>0.047044</td>
      <td>ACD</td>
      <td>67660199.0</td>
      <td>ACD</td>
      <td>0.0470437194800431</td>
      <td>0.03481</td>
      <td>gnomADg_MID</td>
      <td>rs149365469</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>46</th>
      <td>16</td>
      <td>C</td>
      <td>T</td>
      <td>ACD</td>
      <td>67660199</td>
      <td>67660199</td>
      <td>67660199</td>
      <td>0.047044</td>
      <td>ACD</td>
      <td>67660199.0</td>
      <td>ACD</td>
      <td>0.0470437194800431</td>
      <td>0.03481</td>
      <td>gnomADg_MID</td>
      <td>rs149365469</td>
      <td>MODERATE</td>
      <td>missense_variant&amp;NMD_transcript_variant</td>
    </tr>
    <tr>
      <th>50</th>
      <td>3</td>
      <td>G</td>
      <td>GTT</td>
      <td>ZNF717</td>
      <td>75739225</td>
      <td>75739225</td>
      <td>75739225</td>
      <td>0.133680</td>
      <td>ZNF717</td>
      <td>75739225.0</td>
      <td>ZNF717</td>
      <td>0.133680276943534</td>
      <td>0.3291</td>
      <td>gnomADg_AMI</td>
      <td>rs749453662</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
    <tr>
      <th>54</th>
      <td>2</td>
      <td>C</td>
      <td>T</td>
      <td>ANKRD23</td>
      <td>96842402</td>
      <td>96842402</td>
      <td>96842402</td>
      <td>0.041633</td>
      <td>ANKRD23</td>
      <td>96842402.0</td>
      <td>ANKRD23</td>
      <td>0.0416332332240537</td>
      <td>0.04114</td>
      <td>gnomADg_MID</td>
      <td>rs143372458&amp;COSV58413586</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>58</th>
      <td>11</td>
      <td>T</td>
      <td>C</td>
      <td>APOA4</td>
      <td>116820959</td>
      <td>116820959</td>
      <td>116820959</td>
      <td>0.044808</td>
      <td>APOA4</td>
      <td>116820959.0</td>
      <td>APOA4</td>
      <td>0.0448084550310644</td>
      <td>0.02423</td>
      <td>gnomADg_AMI</td>
      <td>rs675&amp;CM031126&amp;COSV63360269</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>66</th>
      <td>2</td>
      <td>A</td>
      <td>G</td>
      <td>MCM6</td>
      <td>135840873</td>
      <td>135840873</td>
      <td>135840873</td>
      <td>0.200291</td>
      <td>MCM6</td>
      <td>135840873.0</td>
      <td>MCM6</td>
      <td>0.200291157523098</td>
      <td>0.1139</td>
      <td>gnomADg_MID</td>
      <td>rs55660827&amp;CR080768</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>71</th>
      <td>5</td>
      <td>A</td>
      <td>C</td>
      <td>GPR151</td>
      <td>146516045</td>
      <td>146516045</td>
      <td>146516045</td>
      <td>0.028187</td>
      <td>GPR151</td>
      <td>146516045.0</td>
      <td>GPR151</td>
      <td>0.0281872578463514</td>
      <td>0.01899</td>
      <td>gnomADg_MID</td>
      <td>rs144066680</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
    <tr>
      <th>73</th>
      <td>4</td>
      <td>A</td>
      <td>C</td>
      <td>DCLK2</td>
      <td>150256161</td>
      <td>150256161</td>
      <td>150256161</td>
      <td>0.000046</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>74</th>
      <td>2</td>
      <td>TCGCA</td>
      <td>T</td>
      <td>NRP2</td>
      <td>205776515</td>
      <td>205776519</td>
      <td>205776519</td>
      <td>0.073512</td>
      <td>NRP2</td>
      <td>205776519.0</td>
      <td>NRP2</td>
      <td>0.0735115431348724</td>
      <td>0.1143</td>
      <td>gnomADg_AMI</td>
      <td>rs527478913</td>
      <td>HIGH</td>
      <td>frameshift_variant</td>
    </tr>
    <tr>
      <th>78</th>
      <td>2</td>
      <td>C</td>
      <td>G</td>
      <td>ACADL</td>
      <td>210205678</td>
      <td>210205678</td>
      <td>210205678</td>
      <td>0.027992</td>
      <td>ACADL</td>
      <td>210205678.0</td>
      <td>ACADL</td>
      <td>0.0279923886379789</td>
      <td>0.008239</td>
      <td>gnomADe_ASJ</td>
      <td>rs146511220&amp;COSV106086073</td>
      <td>MODERATE</td>
      <td>missense_variant</td>
    </tr>
  </tbody>
</table>
</div>




```python
#update_oldTable5_noMod.to_csv("/home/jovyan/aisha/Final_outputs/Table5_update_partb.csv", index = False)
```


```python
# merge with MGI data with new positions
merge_mgiMAF_liftoverPos = pd.merge(liftover_pos,mgiMAF, how='inner', on =['chr', 'start', 'end', 'referenceAllele','alternateAlleles'])
print("Total number of rows: " + str(len(merge_mgiMAF_liftoverPos)))
merge_mgiMAF_liftoverPos.head(3)
```


    ---------------------------------------------------------------------------

    KeyError                                  Traceback (most recent call last)

    Input In [140], in <cell line: 2>()
          1 # merge with MGI data with new positions
    ----> 2 merge_mgiMAF_liftoverPos = pd.merge(liftover_pos,mgiMAF, how='inner', on =['chr', 'start', 'end', 'referenceAllele','alternateAlleles'])
          3 print("Total number of rows: " + str(len(merge_mgiMAF_liftoverPos)))
          4 merge_mgiMAF_liftoverPos.head(3)
    

    File /opt/conda/lib/python3.10/site-packages/pandas/core/reshape/merge.py:111, in merge(left, right, how, on, left_on, right_on, left_index, right_index, sort, suffixes, copy, indicator, validate)
         94 @Substitution("\nleft : DataFrame or named Series")
         95 @Appender(_merge_doc, indents=0)
         96 def merge(
       (...)
        109     validate: str | None = None,
        110 ) -> DataFrame:
    --> 111     op = _MergeOperation(
        112         left,
        113         right,
        114         how=how,
        115         on=on,
        116         left_on=left_on,
        117         right_on=right_on,
        118         left_index=left_index,
        119         right_index=right_index,
        120         sort=sort,
        121         suffixes=suffixes,
        122         indicator=indicator,
        123         validate=validate,
        124     )
        125     return op.get_result(copy=copy)
    

    File /opt/conda/lib/python3.10/site-packages/pandas/core/reshape/merge.py:706, in _MergeOperation.__init__(self, left, right, how, on, left_on, right_on, axis, left_index, right_index, sort, suffixes, indicator, validate)
        699 self._cross = cross_col
        701 # note this function has side effects
        702 (
        703     self.left_join_keys,
        704     self.right_join_keys,
        705     self.join_names,
    --> 706 ) = self._get_merge_keys()
        708 # validate the merge keys dtypes. We may need to coerce
        709 # to avoid incompatible dtypes
        710 self._maybe_coerce_merge_keys()
    

    File /opt/conda/lib/python3.10/site-packages/pandas/core/reshape/merge.py:1182, in _MergeOperation._get_merge_keys(self)
       1178 if lk is not None:
       1179     # Then we're either Hashable or a wrong-length arraylike,
       1180     #  the latter of which will raise
       1181     lk = cast(Hashable, lk)
    -> 1182     left_keys.append(left._get_label_or_level_values(lk))
       1183     join_names.append(lk)
       1184 else:
       1185     # work-around for merge_asof(left_index=True)
    

    File /opt/conda/lib/python3.10/site-packages/pandas/core/generic.py:1849, in NDFrame._get_label_or_level_values(self, key, axis)
       1843     values = (
       1844         self.axes[axis]
       1845         .get_level_values(key)  # type: ignore[assignment]
       1846         ._values
       1847     )
       1848 else:
    -> 1849     raise KeyError(key)
       1851 # Check for duplicates
       1852 if values.ndim > 1:
    

    KeyError: 'start'



```python
maf_Illumina_VEP_variants_final_noMod = maf_Illumina_VEP_variants_final.loc[maf_Illumina_VEP_variants_final['IMPACT']!="MODIFIER"]
```


```python
merge_illuminaMAF_liftoverPos['chr'] = merge_illuminaMAF_liftoverPos['chr'].str.replace('chr', '')
merge_illuminaMAF_liftoverPos['start'] = merge_illuminaMAF_liftoverPos['start'].astype(int)+1
merge_illuminaMAF_liftoverPos['end'] = merge_illuminaMAF_liftoverPos['end'].astype(int)
merge_illuminaMAF_liftoverPos

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
      <th>chr</th>
      <th>referenceAllele</th>
      <th>alternateAlleles</th>
      <th>gene</th>
      <th>start</th>
      <th>end</th>
      <th>AlternateFreq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>47</th>
      <td>1</td>
      <td>CAGCTT</td>
      <td>C</td>
      <td>ESPN</td>
      <td>6445776</td>
      <td>6445781</td>
      <td>0.003783</td>
    </tr>
  </tbody>
</table>
</div>


