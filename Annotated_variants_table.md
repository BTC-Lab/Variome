```python
import os
import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 200)
```


```python
# function to capitalize first letter only
def cap_first(input):
    return input[0].upper() + input[1:]
```


```python
# step 2 All variants
all_varinats = pd.read_csv("../output_files/step1_allVariants_noFiler.csv", delimiter=',')
all_varinats = all_varinats.drop(['chr12.txt_warnings.txt', 'chr17.txt_warnings.txt', 'transcript_ablation', 'chrrange21.txt', 'chrrange22.txt'])
all_varinats = all_varinats.astype(int)
all_varinats.head(2)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr10.txt</th>
      <td>0</td>
      <td>1489</td>
      <td>2129</td>
      <td>2345</td>
      <td>3254</td>
      <td>206</td>
      <td>435</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>2142178</td>
      <td>2235880</td>
      <td>2031</td>
      <td>0</td>
      <td>219427</td>
      <td>6</td>
      <td>0</td>
      <td>2593746</td>
      <td>8495226</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr11.txt</th>
      <td>1</td>
      <td>1938</td>
      <td>2639</td>
      <td>3440</td>
      <td>5583</td>
      <td>380</td>
      <td>709</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>2893846</td>
      <td>2933779</td>
      <td>1501</td>
      <td>0</td>
      <td>227946</td>
      <td>13</td>
      <td>0</td>
      <td>2608430</td>
      <td>12436500</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 41 columns</p>
</div>




```python
#step 3 Novel not in repeat masker region
novel_notRepeatmasker = pd.read_csv("../output_files/step3_novelVariants_notRepeat.csv", delimiter=',', index_col=0)
novel_notRepeatmasker
novel_notRepeatmasker.head(2)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr10_not_noRS.txt</th>
      <td>0</td>
      <td>738</td>
      <td>1123</td>
      <td>1092</td>
      <td>1976</td>
      <td>107</td>
      <td>247</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>583628</td>
      <td>592481</td>
      <td>1080</td>
      <td>0</td>
      <td>81758</td>
      <td>2</td>
      <td>0</td>
      <td>821964</td>
      <td>1004107</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr11_not_noRS.txt</th>
      <td>0</td>
      <td>887</td>
      <td>1213</td>
      <td>1324</td>
      <td>3376</td>
      <td>187</td>
      <td>402</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>794332</td>
      <td>773803</td>
      <td>764</td>
      <td>0</td>
      <td>87159</td>
      <td>0</td>
      <td>0</td>
      <td>812178</td>
      <td>705432</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 41 columns</p>
</div>




```python
#step 4 Novel in repeat masker region
novel_inRepeatmasker = pd.read_csv("../output_files/step4_novelVariants_inRepeat.csv", delimiter=',')
novel_inRepeatmasker
novel_inRepeatmasker.head(2)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr10_in_noRS.txt</th>
      <td>0</td>
      <td>191</td>
      <td>133</td>
      <td>26</td>
      <td>135</td>
      <td>9</td>
      <td>14</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>545223</td>
      <td>575100</td>
      <td>500</td>
      <td>0</td>
      <td>33009</td>
      <td>4</td>
      <td>0</td>
      <td>524616</td>
      <td>4793002</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr11_in_noRS.txt</th>
      <td>0</td>
      <td>257</td>
      <td>241</td>
      <td>54</td>
      <td>192</td>
      <td>13</td>
      <td>13</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>730965</td>
      <td>764373</td>
      <td>361</td>
      <td>0</td>
      <td>32759</td>
      <td>7</td>
      <td>0</td>
      <td>548292</td>
      <td>9052558</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 41 columns</p>
</div>




```python
# all variants in repeat
all_inRepeatmasker = pd.read_csv("../output_files/step5_inRepeat.txt", delimiter=',')
all_inRepeatmasker.head(2)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr10_In.vcf</th>
      <td>0</td>
      <td>332</td>
      <td>244</td>
      <td>66</td>
      <td>229</td>
      <td>19</td>
      <td>25</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1037992</td>
      <td>1095761</td>
      <td>682</td>
      <td>0</td>
      <td>63041</td>
      <td>4</td>
      <td>0</td>
      <td>1002513</td>
      <td>6546667</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr10_out.vcf</th>
      <td>0</td>
      <td>1157</td>
      <td>1885</td>
      <td>2279</td>
      <td>3025</td>
      <td>187</td>
      <td>410</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1104186</td>
      <td>1140119</td>
      <td>1349</td>
      <td>0</td>
      <td>156386</td>
      <td>2</td>
      <td>0</td>
      <td>1591233</td>
      <td>1948559</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 41 columns</p>
</div>




```python
# all variants in repeat
notRepeatmasker_sub = pd.read_csv("../output_files/step5_notRepeat.csv", delimiter=',', names = list(novel_inRepeatmasker.columns))
notRepeatmasker_sub.head(1)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr20_out.vcf</th>
      <td>5</td>
      <td>684</td>
      <td>901</td>
      <td>1163</td>
      <td>1773</td>
      <td>108</td>
      <td>304</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>627011</td>
      <td>646895</td>
      <td>1009</td>
      <td>0</td>
      <td>106054</td>
      <td>0</td>
      <td>0</td>
      <td>855938</td>
      <td>817495</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 41 columns</p>
</div>




```python
#step 5 All variants in repeat masker region
inRepeatmasker = all_inRepeatmasker[~all_inRepeatmasker.index.str.endswith("_out.vcf")]
inRepeatmasker.head(1)#.sum(axis=0)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr10_In.vcf</th>
      <td>0</td>
      <td>332</td>
      <td>244</td>
      <td>66</td>
      <td>229</td>
      <td>19</td>
      <td>25</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1037992</td>
      <td>1095761</td>
      <td>682</td>
      <td>0</td>
      <td>63041</td>
      <td>4</td>
      <td>0</td>
      <td>1002513</td>
      <td>6546667</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 41 columns</p>
</div>




```python
# All variants not in repeat masker region
notRepeatmasker = all_inRepeatmasker[all_inRepeatmasker.index.str.endswith("_out.vcf")]
notRepeatmasker.head(4)#.sum(axis=0)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr10_out.vcf</th>
      <td>0</td>
      <td>1157</td>
      <td>1885</td>
      <td>2279</td>
      <td>3025</td>
      <td>187</td>
      <td>410</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1104186</td>
      <td>1140119</td>
      <td>1349</td>
      <td>0</td>
      <td>156386</td>
      <td>2</td>
      <td>0</td>
      <td>1591233</td>
      <td>1948559</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr11_out.vcf</th>
      <td>1</td>
      <td>1537</td>
      <td>2189</td>
      <td>3319</td>
      <td>5272</td>
      <td>349</td>
      <td>691</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1499863</td>
      <td>1484056</td>
      <td>978</td>
      <td>0</td>
      <td>166051</td>
      <td>0</td>
      <td>0</td>
      <td>1568634</td>
      <td>1358073</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr12_out.vcf</th>
      <td>1</td>
      <td>1317</td>
      <td>2055</td>
      <td>2798</td>
      <td>4618</td>
      <td>433</td>
      <td>535</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1315042</td>
      <td>1336343</td>
      <td>822</td>
      <td>0</td>
      <td>149364</td>
      <td>1</td>
      <td>0</td>
      <td>1391400</td>
      <td>1305352</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr13_out.vcf</th>
      <td>0</td>
      <td>502</td>
      <td>763</td>
      <td>824</td>
      <td>1328</td>
      <td>67</td>
      <td>186</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>601442</td>
      <td>620100</td>
      <td>946</td>
      <td>0</td>
      <td>62046</td>
      <td>3</td>
      <td>0</td>
      <td>840344</td>
      <td>1822901</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>4 rows × 41 columns</p>
</div>




```python
# Known variants  in repeat masker region
inRepeatmasker_Known = pd.read_csv("../output_files/step6_notRepeat_Known.csv", delimiter=',')
inRepeatmasker_Known = inRepeatmasker_Known.drop(['transcript_ablation'])
inRepeatmasker_Known.head(4)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>novel/Table4_repeatmasker_chr1_not_RS.vcf</th>
      <td>0</td>
      <td>1009</td>
      <td>1692</td>
      <td>3578</td>
      <td>2966</td>
      <td>204</td>
      <td>432</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1108160</td>
      <td>1121346</td>
      <td>413</td>
      <td>0</td>
      <td>152109</td>
      <td>0</td>
      <td>0</td>
      <td>1355222</td>
      <td>1145045</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>novel/Table4_repeatmasker_chr2_not_RS.vcf</th>
      <td>1</td>
      <td>853</td>
      <td>1401</td>
      <td>2136</td>
      <td>1822</td>
      <td>145</td>
      <td>279</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>969767</td>
      <td>1028508</td>
      <td>408</td>
      <td>0</td>
      <td>117849</td>
      <td>1</td>
      <td>0</td>
      <td>1268010</td>
      <td>1490990</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>novel/Table4_repeatmasker_chr3_not_RS.vcf</th>
      <td>0</td>
      <td>655</td>
      <td>974</td>
      <td>1710</td>
      <td>2086</td>
      <td>136</td>
      <td>220</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>706644</td>
      <td>766199</td>
      <td>239</td>
      <td>0</td>
      <td>88086</td>
      <td>0</td>
      <td>0</td>
      <td>987484</td>
      <td>1071389</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>novel/Table4_repeatmasker_chr4_not_RS.vcf</th>
      <td>0</td>
      <td>533</td>
      <td>725</td>
      <td>1337</td>
      <td>998</td>
      <td>84</td>
      <td>174</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>570809</td>
      <td>619281</td>
      <td>156</td>
      <td>0</td>
      <td>59172</td>
      <td>1</td>
      <td>0</td>
      <td>713264</td>
      <td>1411726</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>4 rows × 41 columns</p>
</div>




```python
# Known variants  in repeat masker region
inRepeatmasker_Known = pd.read_csv("../output_files/step7_inRepeat_Known.csv", delimiter=',')
inRepeatmasker_Known.head(4)
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
      <th>transcript_ablation</th>
      <th>splice_acceptor_variant</th>
      <th>splice_donor_variant</th>
      <th>stop_gained</th>
      <th>frameshift_variant</th>
      <th>stop_lost</th>
      <th>start_lost</th>
      <th>transcript_amplification</th>
      <th>feature_elongation</th>
      <th>feature_truncation</th>
      <th>...</th>
      <th>upstream_gene_variant</th>
      <th>downstream_gene_variant</th>
      <th>TFBS_ablation</th>
      <th>TFBS_amplification</th>
      <th>TF_binding_site_variant</th>
      <th>regulatory_region_ablation</th>
      <th>regulatory_region_amplification</th>
      <th>regulatory_region_variant</th>
      <th>intergenic_variant</th>
      <th>sequence_variant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>novel/Table4_repeatmasker_chr1_In_RS.vcf</th>
      <td>0</td>
      <td>271</td>
      <td>279</td>
      <td>127</td>
      <td>239</td>
      <td>23</td>
      <td>22</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1135689</td>
      <td>1164746</td>
      <td>375</td>
      <td>0</td>
      <td>62860</td>
      <td>0</td>
      <td>0</td>
      <td>932216</td>
      <td>3081709</td>
      <td>0</td>
    </tr>
    <tr>
      <th>novel/Table4_repeatmasker_chr2_In_RS.vcf</th>
      <td>0</td>
      <td>238</td>
      <td>245</td>
      <td>70</td>
      <td>146</td>
      <td>17</td>
      <td>16</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>878199</td>
      <td>942827</td>
      <td>280</td>
      <td>0</td>
      <td>46559</td>
      <td>2</td>
      <td>0</td>
      <td>807973</td>
      <td>2953257</td>
      <td>0</td>
    </tr>
    <tr>
      <th>novel/Table4_repeatmasker_chr3_In_RS.vcf</th>
      <td>0</td>
      <td>182</td>
      <td>249</td>
      <td>76</td>
      <td>105</td>
      <td>21</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>698275</td>
      <td>772803</td>
      <td>195</td>
      <td>0</td>
      <td>37845</td>
      <td>1</td>
      <td>0</td>
      <td>682002</td>
      <td>2255839</td>
      <td>0</td>
    </tr>
    <tr>
      <th>novel/Table4_repeatmasker_chr4_In_RS.vcf</th>
      <td>0</td>
      <td>141</td>
      <td>188</td>
      <td>43</td>
      <td>152</td>
      <td>16</td>
      <td>7</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>571638</td>
      <td>619046</td>
      <td>137</td>
      <td>0</td>
      <td>31369</td>
      <td>3</td>
      <td>0</td>
      <td>522468</td>
      <td>2513784</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>4 rows × 41 columns</p>
</div>




```python
imapact_sev= pd.read_csv("Impact_sev.txt", delimiter='\t', header=None, index_col=0)
imapact_sev =imapact_sev.rename_axis('impact',index=True)
imapact_sev= imapact_sev.sort_index(ascending=True)
imapact_sev = imapact_sev.rename(columns={1:'severity'})
imapact_sev.index = imapact_sev.index.to_series().apply(cap_first).str.replace('_', ' ' )
```


```python
# table of combined sums of the number of variants
combinedsums = pd.DataFrame({'All': all_varinats.sum(axis=0),'In repeat regions':inRepeatmasker.sum(axis=0), 'In repeat regions "Novel"': novel_inRepeatmasker.sum(axis=0), 'Not in repeat regions "Novel"':novel_notRepeatmasker.sum(axis=0)}) 
combinedsums.index = combinedsums.index.to_series().apply(cap_first).str.replace('_', ' ' )
combinedsums['Not in repeat regions "Novel"'] = combinedsums['Not in repeat regions "Novel"'].astype("int")
combinedsums['Not in repeat regions'] = combinedsums['All'] - combinedsums['In repeat regions']
combinedsums =combinedsums.rename_axis('impact',index=True)
combinedsums = combinedsums.sort_index(ascending=True)
```


```python
Table4 = combinedsums.join(imapact_sev)
Table4['Not in repeat regions'] = Table4['All'] - Table4['In repeat regions']
impact_order =  ["HIGH", "MODERATE", "LOW", "MODIFIER"]
Table4['severity'] = pd.Categorical(Table4['severity'], categories=impact_order)
Table4=Table4.sort_values(by=['severity', 'All', 'In repeat regions','In repeat regions "Novel"','Not in repeat regions "Novel"'] , ascending=[True, False, False, False, False])
Table4
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
      <th>All</th>
      <th>In repeat regions</th>
      <th>In repeat regions "Novel"</th>
      <th>Not in repeat regions "Novel"</th>
      <th>Not in repeat regions</th>
      <th>severity</th>
    </tr>
    <tr>
      <th>impact</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Frameshift variant</th>
      <td>100058</td>
      <td>9578</td>
      <td>6433</td>
      <td>60749</td>
      <td>90480</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Stop gained</th>
      <td>67541</td>
      <td>2595</td>
      <td>1310</td>
      <td>32768</td>
      <td>64946</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Splice donor variant</th>
      <td>53057</td>
      <td>8308</td>
      <td>4658</td>
      <td>26663</td>
      <td>44749</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Splice acceptor variant</th>
      <td>41162</td>
      <td>8751</td>
      <td>5507</td>
      <td>20815</td>
      <td>32411</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Start lost</th>
      <td>11908</td>
      <td>555</td>
      <td>375</td>
      <td>7058</td>
      <td>11353</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Stop lost</th>
      <td>6649</td>
      <td>780</td>
      <td>441</td>
      <td>3493</td>
      <td>5869</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Transcript ablation</th>
      <td>40</td>
      <td>4</td>
      <td>2</td>
      <td>30</td>
      <td>36</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Feature elongation</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Feature truncation</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Transcript amplification</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Missense variant</th>
      <td>2145817</td>
      <td>66084</td>
      <td>29995</td>
      <td>839368</td>
      <td>2079733</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Inframe deletion</th>
      <td>41380</td>
      <td>6841</td>
      <td>1662</td>
      <td>17610</td>
      <td>34539</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Inframe insertion</th>
      <td>31242</td>
      <td>9486</td>
      <td>4926</td>
      <td>14910</td>
      <td>21756</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Protein altering variant</th>
      <td>2711</td>
      <td>349</td>
      <td>302</td>
      <td>2040</td>
      <td>2362</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Synonymous variant</th>
      <td>1082712</td>
      <td>31936</td>
      <td>11984</td>
      <td>321638</td>
      <td>1050776</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice polypyrimidine tract variant</th>
      <td>494512</td>
      <td>88600</td>
      <td>45217</td>
      <td>191553</td>
      <td>405912</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice region variant</th>
      <td>436389</td>
      <td>67837</td>
      <td>34526</td>
      <td>178253</td>
      <td>368552</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice donor region variant</th>
      <td>87462</td>
      <td>11074</td>
      <td>6195</td>
      <td>43792</td>
      <td>76388</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice donor 5th base variant</th>
      <td>33389</td>
      <td>4985</td>
      <td>2540</td>
      <td>14866</td>
      <td>28404</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Stop retained variant</th>
      <td>2331</td>
      <td>261</td>
      <td>120</td>
      <td>831</td>
      <td>2070</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Incomplete terminal codon variant</th>
      <td>1497</td>
      <td>45</td>
      <td>23</td>
      <td>503</td>
      <td>1452</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Start retained variant</th>
      <td>118</td>
      <td>3</td>
      <td>1</td>
      <td>49</td>
      <td>115</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Intron variant</th>
      <td>191104957</td>
      <td>102665745</td>
      <td>54905588</td>
      <td>47256251</td>
      <td>88439212</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Intergenic variant</th>
      <td>190823822</td>
      <td>152297552</td>
      <td>114829872</td>
      <td>21245508</td>
      <td>38526270</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Non coding transcript variant</th>
      <td>126272380</td>
      <td>68951066</td>
      <td>37254622</td>
      <td>30885929</td>
      <td>57321314</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Downstream gene variant</th>
      <td>59587422</td>
      <td>29943700</td>
      <td>16548264</td>
      <td>16207237</td>
      <td>29643722</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Upstream gene variant</th>
      <td>57520146</td>
      <td>28820980</td>
      <td>16088250</td>
      <td>15896133</td>
      <td>28699166</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Regulatory region variant</th>
      <td>55339899</td>
      <td>22740521</td>
      <td>12447354</td>
      <td>17393322</td>
      <td>32599378</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>NMD transcript variant</th>
      <td>49121558</td>
      <td>25556115</td>
      <td>13321535</td>
      <td>12079318</td>
      <td>23565443</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Non coding transcript exon variant</th>
      <td>10403124</td>
      <td>2729075</td>
      <td>1425573</td>
      <td>3923469</td>
      <td>7674049</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>TF binding site variant</th>
      <td>4920885</td>
      <td>1529165</td>
      <td>863109</td>
      <td>1830014</td>
      <td>3391720</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>3 prime UTR variant</th>
      <td>4918843</td>
      <td>1160803</td>
      <td>589683</td>
      <td>1814098</td>
      <td>3758040</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>5 prime UTR variant</th>
      <td>1402780</td>
      <td>235774</td>
      <td>130752</td>
      <td>633907</td>
      <td>1167006</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>TFBS ablation</th>
      <td>38811</td>
      <td>14615</td>
      <td>11369</td>
      <td>18695</td>
      <td>24196</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Mature miRNA variant</th>
      <td>7797</td>
      <td>1770</td>
      <td>1100</td>
      <td>3216</td>
      <td>6027</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Coding sequence variant</th>
      <td>5770</td>
      <td>529</td>
      <td>324</td>
      <td>2700</td>
      <td>5241</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Regulatory region ablation</th>
      <td>135</td>
      <td>90</td>
      <td>73</td>
      <td>38</td>
      <td>45</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Coding transcript variant</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Regulatory region amplification</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Sequence variant</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>TFBS amplification</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>MODIFIER</td>
    </tr>
  </tbody>
</table>
</div>




```python
Table4_cleaned = combinedsums[combinedsums.values.sum(axis=1)!=0]

Table4_cleaned = Table4_cleaned.join(imapact_sev)
impact_order =  ["HIGH", "MODERATE", "LOW", "MODIFIER"]
Table4_cleaned['severity'] = pd.Categorical(Table4_cleaned['severity'], categories=impact_order)
Table4_cleaned = Table4_cleaned.sort_values(by=['severity', 'All', 'In repeat regions','In repeat regions "Novel"','Not in repeat regions "Novel"'] , ascending=[True, False, False, False, False])
```


```python
z = Table4_cleaned['In repeat regions']/Table4_cleaned['All']*100
Table4_cleaned.insert(2, 'In RR', z.round(0))

z = Table4_cleaned['In repeat regions "Novel"']/Table4_cleaned['All']*100
Table4_cleaned.insert(4, 'In RR Novel', z.round(0))

z = Table4_cleaned['Not in repeat regions "Novel"']/Table4_cleaned['All']*100
Table4_cleaned.insert(6, 'Not RR Novel', z.round(0))

z = Table4_cleaned['Not in repeat regions']/Table4_cleaned['All']*100
Table4_cleaned.insert(8, 'Not RR', z.round(0))
```


```python
Table4_cleaned
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
      <th>All</th>
      <th>In repeat regions</th>
      <th>In RR</th>
      <th>In repeat regions "Novel"</th>
      <th>In RR Novel</th>
      <th>Not in repeat regions "Novel"</th>
      <th>Not RR Novel</th>
      <th>Not in repeat regions</th>
      <th>Not RR</th>
      <th>severity</th>
    </tr>
    <tr>
      <th>impact</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Frameshift variant</th>
      <td>100058</td>
      <td>9578</td>
      <td>10.0</td>
      <td>6433</td>
      <td>6.0</td>
      <td>60749</td>
      <td>61.0</td>
      <td>90480</td>
      <td>90.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Stop gained</th>
      <td>67541</td>
      <td>2595</td>
      <td>4.0</td>
      <td>1310</td>
      <td>2.0</td>
      <td>32768</td>
      <td>49.0</td>
      <td>64946</td>
      <td>96.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Splice donor variant</th>
      <td>53057</td>
      <td>8308</td>
      <td>16.0</td>
      <td>4658</td>
      <td>9.0</td>
      <td>26663</td>
      <td>50.0</td>
      <td>44749</td>
      <td>84.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Splice acceptor variant</th>
      <td>41162</td>
      <td>8751</td>
      <td>21.0</td>
      <td>5507</td>
      <td>13.0</td>
      <td>20815</td>
      <td>51.0</td>
      <td>32411</td>
      <td>79.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Start lost</th>
      <td>11908</td>
      <td>555</td>
      <td>5.0</td>
      <td>375</td>
      <td>3.0</td>
      <td>7058</td>
      <td>59.0</td>
      <td>11353</td>
      <td>95.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Stop lost</th>
      <td>6649</td>
      <td>780</td>
      <td>12.0</td>
      <td>441</td>
      <td>7.0</td>
      <td>3493</td>
      <td>53.0</td>
      <td>5869</td>
      <td>88.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Transcript ablation</th>
      <td>40</td>
      <td>4</td>
      <td>10.0</td>
      <td>2</td>
      <td>5.0</td>
      <td>30</td>
      <td>75.0</td>
      <td>36</td>
      <td>90.0</td>
      <td>HIGH</td>
    </tr>
    <tr>
      <th>Missense variant</th>
      <td>2145817</td>
      <td>66084</td>
      <td>3.0</td>
      <td>29995</td>
      <td>1.0</td>
      <td>839368</td>
      <td>39.0</td>
      <td>2079733</td>
      <td>97.0</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Inframe deletion</th>
      <td>41380</td>
      <td>6841</td>
      <td>17.0</td>
      <td>1662</td>
      <td>4.0</td>
      <td>17610</td>
      <td>43.0</td>
      <td>34539</td>
      <td>83.0</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Inframe insertion</th>
      <td>31242</td>
      <td>9486</td>
      <td>30.0</td>
      <td>4926</td>
      <td>16.0</td>
      <td>14910</td>
      <td>48.0</td>
      <td>21756</td>
      <td>70.0</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Protein altering variant</th>
      <td>2711</td>
      <td>349</td>
      <td>13.0</td>
      <td>302</td>
      <td>11.0</td>
      <td>2040</td>
      <td>75.0</td>
      <td>2362</td>
      <td>87.0</td>
      <td>MODERATE</td>
    </tr>
    <tr>
      <th>Synonymous variant</th>
      <td>1082712</td>
      <td>31936</td>
      <td>3.0</td>
      <td>11984</td>
      <td>1.0</td>
      <td>321638</td>
      <td>30.0</td>
      <td>1050776</td>
      <td>97.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice polypyrimidine tract variant</th>
      <td>494512</td>
      <td>88600</td>
      <td>18.0</td>
      <td>45217</td>
      <td>9.0</td>
      <td>191553</td>
      <td>39.0</td>
      <td>405912</td>
      <td>82.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice region variant</th>
      <td>436389</td>
      <td>67837</td>
      <td>16.0</td>
      <td>34526</td>
      <td>8.0</td>
      <td>178253</td>
      <td>41.0</td>
      <td>368552</td>
      <td>84.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice donor region variant</th>
      <td>87462</td>
      <td>11074</td>
      <td>13.0</td>
      <td>6195</td>
      <td>7.0</td>
      <td>43792</td>
      <td>50.0</td>
      <td>76388</td>
      <td>87.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Splice donor 5th base variant</th>
      <td>33389</td>
      <td>4985</td>
      <td>15.0</td>
      <td>2540</td>
      <td>8.0</td>
      <td>14866</td>
      <td>45.0</td>
      <td>28404</td>
      <td>85.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Stop retained variant</th>
      <td>2331</td>
      <td>261</td>
      <td>11.0</td>
      <td>120</td>
      <td>5.0</td>
      <td>831</td>
      <td>36.0</td>
      <td>2070</td>
      <td>89.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Incomplete terminal codon variant</th>
      <td>1497</td>
      <td>45</td>
      <td>3.0</td>
      <td>23</td>
      <td>2.0</td>
      <td>503</td>
      <td>34.0</td>
      <td>1452</td>
      <td>97.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Start retained variant</th>
      <td>118</td>
      <td>3</td>
      <td>3.0</td>
      <td>1</td>
      <td>1.0</td>
      <td>49</td>
      <td>42.0</td>
      <td>115</td>
      <td>97.0</td>
      <td>LOW</td>
    </tr>
    <tr>
      <th>Intron variant</th>
      <td>191104957</td>
      <td>102665745</td>
      <td>54.0</td>
      <td>54905588</td>
      <td>29.0</td>
      <td>47256251</td>
      <td>25.0</td>
      <td>88439212</td>
      <td>46.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Intergenic variant</th>
      <td>190823822</td>
      <td>152297552</td>
      <td>80.0</td>
      <td>114829872</td>
      <td>60.0</td>
      <td>21245508</td>
      <td>11.0</td>
      <td>38526270</td>
      <td>20.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Non coding transcript variant</th>
      <td>126272380</td>
      <td>68951066</td>
      <td>55.0</td>
      <td>37254622</td>
      <td>30.0</td>
      <td>30885929</td>
      <td>24.0</td>
      <td>57321314</td>
      <td>45.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Downstream gene variant</th>
      <td>59587422</td>
      <td>29943700</td>
      <td>50.0</td>
      <td>16548264</td>
      <td>28.0</td>
      <td>16207237</td>
      <td>27.0</td>
      <td>29643722</td>
      <td>50.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Upstream gene variant</th>
      <td>57520146</td>
      <td>28820980</td>
      <td>50.0</td>
      <td>16088250</td>
      <td>28.0</td>
      <td>15896133</td>
      <td>28.0</td>
      <td>28699166</td>
      <td>50.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Regulatory region variant</th>
      <td>55339899</td>
      <td>22740521</td>
      <td>41.0</td>
      <td>12447354</td>
      <td>22.0</td>
      <td>17393322</td>
      <td>31.0</td>
      <td>32599378</td>
      <td>59.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>NMD transcript variant</th>
      <td>49121558</td>
      <td>25556115</td>
      <td>52.0</td>
      <td>13321535</td>
      <td>27.0</td>
      <td>12079318</td>
      <td>25.0</td>
      <td>23565443</td>
      <td>48.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Non coding transcript exon variant</th>
      <td>10403124</td>
      <td>2729075</td>
      <td>26.0</td>
      <td>1425573</td>
      <td>14.0</td>
      <td>3923469</td>
      <td>38.0</td>
      <td>7674049</td>
      <td>74.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>TF binding site variant</th>
      <td>4920885</td>
      <td>1529165</td>
      <td>31.0</td>
      <td>863109</td>
      <td>18.0</td>
      <td>1830014</td>
      <td>37.0</td>
      <td>3391720</td>
      <td>69.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>3 prime UTR variant</th>
      <td>4918843</td>
      <td>1160803</td>
      <td>24.0</td>
      <td>589683</td>
      <td>12.0</td>
      <td>1814098</td>
      <td>37.0</td>
      <td>3758040</td>
      <td>76.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>5 prime UTR variant</th>
      <td>1402780</td>
      <td>235774</td>
      <td>17.0</td>
      <td>130752</td>
      <td>9.0</td>
      <td>633907</td>
      <td>45.0</td>
      <td>1167006</td>
      <td>83.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>TFBS ablation</th>
      <td>38811</td>
      <td>14615</td>
      <td>38.0</td>
      <td>11369</td>
      <td>29.0</td>
      <td>18695</td>
      <td>48.0</td>
      <td>24196</td>
      <td>62.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Mature miRNA variant</th>
      <td>7797</td>
      <td>1770</td>
      <td>23.0</td>
      <td>1100</td>
      <td>14.0</td>
      <td>3216</td>
      <td>41.0</td>
      <td>6027</td>
      <td>77.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Coding sequence variant</th>
      <td>5770</td>
      <td>529</td>
      <td>9.0</td>
      <td>324</td>
      <td>6.0</td>
      <td>2700</td>
      <td>47.0</td>
      <td>5241</td>
      <td>91.0</td>
      <td>MODIFIER</td>
    </tr>
    <tr>
      <th>Regulatory region ablation</th>
      <td>135</td>
      <td>90</td>
      <td>67.0</td>
      <td>73</td>
      <td>54.0</td>
      <td>38</td>
      <td>28.0</td>
      <td>45</td>
      <td>33.0</td>
      <td>MODIFIER</td>
    </tr>
  </tbody>
</table>
</div>




```python

```
